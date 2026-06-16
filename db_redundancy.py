#!/usr/bin/env python3
import os
import sys
import subprocess
import argparse
import logging
import shutil
from pathlib import Path
from Bio import SeqIO
from datetime import datetime
from collections import defaultdict
import re
import glob

"""
A Program desgined to reduce the redundancy in metaproteomics protein database
Run after downloading protein sequences with UniprotDownloader_v3.sh
- It uses CD-HIT to remove the redundancy
- cluster-id - for removing over-all duplicate sequences
- species-id - for removing species-level duplicates - extra layer
"""

# -----------------------------
# Logging setup
# -----------------------------
def setup_logger(verbose=False):
    logger = logging.getLogger("DBPipeline")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    if not logger.handlers:
        fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        ch = logging.StreamHandler(sys.stdout)
        ch.setFormatter(fmt)
        logger.addHandler(ch)
    return logger

logger = setup_logger()

# -----------------------------
# Utility functions
# -----------------------------
def run_cmd(cmd, timeout=3600):
    logger.info(f"Running: {' '.join(map(str, cmd))}")
    try:
        res = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=timeout)
        if res.stdout:
            logger.debug(res.stdout)
        return res
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed ({e.returncode}): {e.stderr}")
        raise
    except subprocess.TimeoutExpired:
        logger.error(f"Command timed out after {timeout}s")
        raise

def check_dependencies():
    for dep in ['awk', 'bash']:
        if not shutil.which(dep):
            raise EnvironmentError(f"Missing dependency: {dep}")
        
def species_from_header(header: str) -> str:
    """Extract species name from '[species=...]' tag in the header."""
    m = re.search(r"\[species\s*=\s*([^\]]+)\]", header)
    if m:
        return m.group(1).strip()
    return "Unknown"

def combine_species_fastas(input_dir: Path, combined_out: Path):
    """
    Combine all FASTA files from input_dir into one file.
    Adds species name (from filename) as a tag in each header.
    """
    input_dir = Path(input_dir)
    #fasta_files = sorted(input_dir.glob("*.fa*"))
    fasta_files = [
        Path(p) for p in glob.glob(os.path.join(str(input_dir), "*.fa*"))]

    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found in {input_dir}")

    logger.info(f"Combining {len(fasta_files)} species FASTA files from {input_dir}")

    with open(combined_out, "w") as out:
        for fasta_file in fasta_files:
            species_name = fasta_file.stem.replace(" ", "_").replace("-", "_")
            for record in SeqIO.parse(fasta_file, "fasta"):
                # Append species info in header
                record.description = f"{record.description} [species={species_name}]"
                out.write(f">{record.description}\n{record.seq}\n")

    logger.info(f"Combined input written to {combined_out}")

def cleanup_file(fasta_file, keep_intermediate):
    if fasta_file and Path(fasta_file).exists() and not keep_intermediate:
        Path(fasta_file).unlink(missing_ok=True)
        logger.debug(f"Removed intermediate file: {fasta_file}")

# -----------------------------
# Species collapse helpers
# -----------------------------
def safe_mkdir(path):
    path.mkdir(parents=True, exist_ok=True)

def species_from_header(header):
    """Extract species tag from header like '[species=Allistipes_sp]'."""
    if "[species=" in header:
        sp = header.split("[species=")[-1].split("]")[0]
        return sp.strip()
    return "Unknown"

def split_by_species(infile: str, outdir: Path):
    """Split combined FASTA into per-species FASTA files, preserving headers."""
    safe_mkdir(outdir)
    species_records = defaultdict(list)
    header_maps = {}

    for record in SeqIO.parse(infile, "fasta"):
        sp = species_from_header(record.description)
        if sp == "Unknown":
            logger.warning(f"No species tag found in header: {record.description[:60]}")
        species_records[sp].append(record)

    for sp, records in species_records.items():
        outfile = outdir / f"{sp}.fasta"
        header_map = {}
        with open(outfile, "w") as f:
            for i, rec in enumerate(records, start=1):
                temp_id = f"{sp}_seq{i}"
                header_map[temp_id] = rec.description
                # Replace only the ID, not the full description
                f.write(f">{temp_id}\n{rec.seq}\n")
        header_maps[sp] = header_map

    map_dir = outdir / "maps"
    safe_mkdir(map_dir)
    for sp, hmap in header_maps.items():
        with open(map_dir / f"{sp}.map", "w") as m:
            for tid, header in hmap.items():
                m.write(f"{tid}\t{header}\n")

    logger.info(f"Split {len(species_records)} species bins into {outdir}")
    return list(species_records.keys())

# Pipeline class
# -----------------------------
class DatabaseCurationPipeline:
    def __init__(self, cfg):
        self.cfg = cfg

    def s1_dereplicate(self):
        logger.info("[Step 1] Dereplication...")
        input_fasta = self.cfg.input_fasta
        output_fasta = self.cfg.derep_fasta
        seq_dict = {}
        for rec in SeqIO.parse(input_fasta, "fasta"):
            seq_str = str(rec.seq)
            if seq_str not in seq_dict:
                seq_dict[seq_str] = rec.description
        
        count = len(seq_dict)
        with open(output_fasta, "w") as out:
            for seq, header in seq_dict.items():
                out.write(f">{header}\n{seq}\n")
        
        logger.info(f"Step 1 complete: {count} unique sequences")

    def s2_cluster99(self):
        logger.info(f"[Step 2] Clustering at {self.cfg.cluster_identity*100:.0f}% identity...")
        input_fasta = self.cfg.derep_fasta
        output_fasta = self.cfg.cluster99_fasta

        run_cmd([
            self.cfg.cdhit_bin,
            "-i", input_fasta,
            "-c", str(self.cfg.cluster_identity),
            "-o", output_fasta,
            "-T", str(self.cfg.threads),
            "-M", "0",
            "-d", "0"
        ])

        # Count clusters
        count = sum(1 for _ in SeqIO.parse(output_fasta, "fasta"))
        logger.info(f"Step 2 complete: {count} clusters")
        
        cleanup_file(input_fasta, self.cfg.keep_intermediate)

    def run_cdhit(self, infile: Path, outprefix: Path, identity: float, threads: int) -> Path:
        """
        Run CD-HIT clustering for protein sequences.
        Returns:
        Path: Path to the CD-HIT .clstr file.
        """
        cmd = [
            str(self.cfg.cdhit_bin),
            "-i", str(infile),
            "-o", str(outprefix),
            "-c", str(identity),
            "-n", "5",        # word size for high-identity protein clustering
            "-T", str(threads),
            "-M", "0",        # use all available memory
            "-d", "0"         # keep full headers
        ]
        run_cmd(cmd)  # Assumes run_cmd handles logging, errors, and execution
        return Path(f"{outprefix}.fasta")

    def s3_species_collapse(self):
        """Step 3: Species-level collapse with CD-HIT per species bins."""
        logger.info(f"[Step 3] Species-level collapse with CD-HIT at {self.cfg.species_identity*100:.0f}% identity...")
        input_fasta = Path(self.cfg.cluster99_fasta)
        output_fasta = Path(self.cfg.species_collapsed_fasta)
        workdir = Path("species_bins")
        safe_mkdir(workdir)

        species_list = split_by_species(input_fasta, workdir)
        merged_count = 0

        with open(output_fasta, "w") as merged:
            for sp in species_list:
                sp_fa = workdir / f"{sp}.fasta"
                seq_count = sum(1 for _ in SeqIO.parse(sp_fa, "fasta"))

                if seq_count == 0:
                    logger.warning(f"No sequences for species {sp}, skipping")
                    continue
                elif seq_count == 1:
                    logger.info(f"Species {sp} has only 1 sequence, skipping CD-HIT")
                    with open(sp_fa) as f:
                        merged.write(f.read())
                    merged_count += 1
                    continue

                outprefix = workdir / f"{sp}_cdhit"
                # Avoid identity=1.0 to prevent CD-HIT from failing
                identity = min(self.cfg.species_identity, 0.999)

                # Run CD-HIT
                self.run_cdhit(sp_fa, outprefix, identity=identity, threads=self.cfg.threads)

                # Clustered sequences are written to outprefix (FASTA), not .clstr
                cluster_fasta = outprefix
                if not cluster_fasta.exists():
                    logger.warning(f"CD-HIT failed or produced no output for {sp}")
                    continue
                
                map_file = workdir / "maps" / f"{sp}.map"
                header_map = {}
                with open(map_file) as m:
                    for line in m:
                        tid, header = line.strip().split("\t", 1)
                        header_map[tid] = header

                for rec in SeqIO.parse(cluster_fasta, "fasta"):
                    full_header = header_map.get(rec.id, rec.id)
                    # Remove [species=...] pattern from header
                    full_header = re.sub(r"\[species=[^\]]+\]", "", full_header).strip()
                    merged.write(f">{full_header}\n{rec.seq}\n")
                    merged_count += 1

        logger.info(f"Step 3 complete: {merged_count} sequences after species collapse")

        cleanup_file(input_fasta, self.cfg.keep_intermediate)
        if not self.cfg.keep_intermediate:
            shutil.rmtree(workdir, ignore_errors=True)

    def finalize(self):
        """Write the final curated DB, removing sequences with MULTISPECIES in the header."""
        final_source = self.cfg.species_collapsed_fasta
        with open(final_source) as inp, open(self.cfg.output_fasta, "w") as out:
            kept_count = 0
            for rec in SeqIO.parse(inp, "fasta"):
                if "MULTISPECIES:" in rec.description:
                    continue  # skip sequences with MULTISPECIES
                elif "uncultured Clostridium sp" in rec.description:
                    continue
                out.write(f">{rec.description}\n{rec.seq}\n")
                kept_count += 1

        logger.info(f"Final curated DB written to: {self.cfg.output_fasta} ({kept_count} sequences kept)")
        cleanup_file(final_source, self.cfg.keep_intermediate)

    def run_pipeline(self):
        start = datetime.now()

        # If input is a directory, combine all species FASTA files
        input_path = Path(self.cfg.input_fasta)
        if input_path.is_dir():
            combined = Path("combined_input.fasta")
            combine_species_fastas(input_path, combined)
            self.cfg.input_fasta = str(combined)
            logger.info(f"Using combined input: {combined}")
        elif not input_path.exists():
            raise FileNotFoundError(f"Input not found: {input_path}")

        if Path(self.cfg.cluster99_fasta).exists() and not self.cfg.keep_intermediate:
            Path(self.cfg.cluster99_fasta).unlink()

        # Count sequences
        input_count = sum(1 for _ in SeqIO.parse(self.cfg.input_fasta, "fasta"))
        logger.info(f"Starting pipeline with {input_count} input sequences")

        self.s1_dereplicate() 
        self.s2_cluster99() 
        self.s3_species_collapse() 
        self.finalize() 
        elapsed = datetime.now() - start 
        logger.info(f"Pipeline completed in {elapsed}")

# -----------------------------
# Configuration
# -----------------------------
class PipelineConfig:
    def __init__(self, args):
        self.input_fasta = args.input
        self.output_fasta = args.output
        self.derep_fasta = "step1_derep.fasta"
        self.cluster99_fasta = "step2_cluster99.fasta"
        self.species_collapsed_fasta = "step3_species_collapsed.fasta"
        self.threads = args.threads
        self.cdhit_bin = args.cdhit_bin
        self.keep_intermediate = args.keep_intermediate
        self.cluster_identity = args.cluster_identity
        self.species_identity = args.species_identity

# Argument parser
# -----------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Database curation pipeline for protein sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-i", "--input", required=True,
                   help="Input FASTA file or directory containing per-species FASTA files")
    parser.add_argument("-o", "--output", default="final_curated_db.fasta",
                       help="Output FASTA file (default: final_curated_db.fasta)")
    parser.add_argument("-t", "--threads", type=int, default=os.cpu_count() or 4,
                       help="Number of threads (default: all available cores)")
    parser.add_argument("--cluster-id", dest="cluster_identity", type=float, default=0.99,
                       help="Identity threshold for step 2 clustering (default: 0.99)")
    parser.add_argument("--species-id", dest="species_identity", type=float, default=1.0,
                       help="Identity threshold for species-level collapse (default: 1.0)")
    parser.add_argument("--cdhit-bin", default="cd-hit",
                       help="Path to cdhit binary (default: cd-hit)")
    parser.add_argument("--keep-intermediate", action="store_true",
                       help="Keep intermediate files")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Verbose output")
    return parser.parse_args()

# Main #
# -----------------------------
def main():
    args = parse_args()
    global logger
    logger = setup_logger(verbose=args.verbose)

    if Path(args.output).exists():
        logger.warning(f"Overwriting output file: {args.output}")

    check_dependencies()
    config = PipelineConfig(args)
    pipeline = DatabaseCurationPipeline(config)
    pipeline.run_pipeline()

if __name__ == "__main__":
    main()