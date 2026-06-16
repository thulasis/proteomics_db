
#!/bin/bash
#
# A Program to download the bacterial proteomes (FASTA Sequences) from NCBI/Uniprot
#To construct protein database for MS-GF+ search engine
# Description: Exhaustive search in different Databases NCBI/and different Uniprot databases
# Method: Uses Webscraping and API integration
# Exhaustive search and provides less redundant protein sequences
# Author: Tulasi Rao Relangi, PhD

# Usage:
#   ./Uniprot_downloader_v3.sh -q queries.txt [-o output_dir] [-m min_count]
set -o pipefail

# -------- Defaults --------
QUERY_FILE=""
OUTPUT_DIR="fasta_downloads"
MIN_COUNT=1000

# -------- Logging --------
log()   { printf "[%s] %s\n" "$(date +%H:%M:%S)" "$*" >&2; }
die()   { log "ERROR: $*"; exit 1; }

# -------- Dependency check --------
need_bins=(curl jq awk grep sed unzip gunzip file)
for b in "${need_bins[@]}"; do
  command -v "$b" >/dev/null 2>&1 || die "Missing dependency: $b"
done

# -------- Usage --------
usage() {
  cat <<EOF
Usage: $0 -q <query_file> [-o <output_dir>] [-m <min_count>]

Options:
  -q   File containing query names (required)
  -o   Output directory (default: ${OUTPUT_DIR})
  -m   Minimum required sequences per query to avoid warning (default: ${MIN_COUNT})

The script tries NCBI (taxonomy → assembly → proteins) and, if needed, falls back
to UniProtKB, then UniRef, then UniParc. Outputs one FASTA per query.
EOF
  exit 1
}

# -------- Normalize spaces (convert Unicode NBSP/U+00A0 and NNBSP/U+202F to ASCII) --------
normalize_spaces() {
  perl -CS -pe 's/\x{00A0}/ /g; s/\x{202F}/ /g' 2>/dev/null || sed 's/\xc2\xa0/ /g; s/\xe2\x80\xaf/ /g'
}

# -------- URL encoding (simple) --------
urlenc() {
  # Encode spaces and () minimally for UniProt queries; adjust if needed
  printf "%s" "$1" | sed 's/ /%20/g; s/(/%28/g; s/)/%29/g'
}

# -------- FASTA formatting: one sequence line per header --------
format_fasta() {
  # format_fasta <input> <output>
  local input="$1" output="$2"
  [[ -s "$input" ]] || return 1
  awk '
    BEGIN { have=0; seq="" }
    /^>/ {
      if (have) {
        gsub(/[ \t\r]/,"",seq);
        print prev "\n" seq
      }
      prev=$0; have=1; seq=""
      next
    }
    { if (NF) seq = seq $0 }
    END {
      if (have) {
        gsub(/[ \t\r]/,"",seq);
        print prev "\n" seq
      }
    }
  ' "$input" >"${output}.tmp" && mv "${output}.tmp" "$output"
}

# -------- Merge + canonicalize --------
merge_fastas() {
  # merge_fastas <output> <input1> [input2 ...]
  local out="$1"; shift
  : >"$out"
  local f
  for f in "$@"; do
    [[ -s "$f" ]] || continue
    format_fasta "$f" "${f}.fmt" || continue
    cat "${f}.fmt" >>"$out"
    rm -f "${f}.fmt"
  done
  [[ -s "$out" ]] && format_fasta "$out" "${out}.tmp" && mv "${out}.tmp" "$out"
}

# -------- Deduplicate by header (keep first) --------
dedup_by_header() {
  # dedup_by_header <input> <output>
  awk '
    /^>/ { h=$0; if (!seen[h]++) { print h; print_next=1 } else { print_next=0 } next }
    { if (print_next) print }
  ' "$1" >"$2.tmp" && mv "$2.tmp" "$2"
}

# -------- Count sequences (headers) --------
count_seqs() {
  # count_seqs <file>
  [[ -s "$1" ]] || { echo 0; return; }
  grep -c '^>' "$1" 2>/dev/null || echo 0
}

# -------- NCBI: get taxonomy id from name --------
get_taxon_id() {
  local query="$1"
  local q
  q="$(printf "%s" "$query" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//;s/[[:space:]]\+/ /g' | sed 's/ /%20/g')"
  local url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=${q}&retmode=json"
  curl -sSL -A "Mozilla/5.0" "$url" | jq -r '.esearchresult.idlist[0] // empty'
}

# -------- NCBI: get assembly accession from taxid (Datasets v2) --------
get_assembly_accession() {
  local taxid="$1"
  local url="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/${taxid}/dataset_report"
  local acc
  acc="$(curl -sSL -H "Accept: application/json" "$url" \
        | jq -r '(.reports[0].accession // .reports[0].current_accession
                 // .assemblies.refseq_assemblies[0].assembly_accession
                 // .assemblies.genbank_assemblies[0].assembly_accession
                 // .assemblies[0].assembly_accession) // empty')"
  printf "%s" "$acc"
}

# -------- NCBI: download proteins from assembly --------
download_proteins_from_assembly() {
  # download_proteins_from_assembly <assembly> <out_fasta> <safe_name>
  local assembly="$1" out="$2" safe="$3"
  local tmpd zip hydrated_url ftp_base listing asm_dir proturl

  tmpd="$(mktemp -d "dl_${safe}_XXXXXX")" || return 1
  zip="${tmpd}/${safe}.zip"

  # Try FULLY_HYDRATED datasets zip
  hydrated_url="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${assembly}/download?include_annotation_type=PROT_FASTA,GENOME_FASTA&hydrated=FULLY_HYDRATED"
  if curl -sSLf "$hydrated_url" -o "$zip"; then
    if file "$zip" 2>/dev/null | grep -q "Zip archive"; then
      unzip -q -o "$zip" -d "$tmpd" 2>/dev/null || true
      # pick *protein*.faa|fa candidates
      mapfile -t protfiles < <(find "$tmpd" -type f \( -iname "*protein*.faa" -o -iname "*protein*.fa" -o -iname "*.faa" -o -iname "*prot*.fa" \) 2>/dev/null)
      if [[ ${#protfiles[@]} -gt 0 ]]; then
        merge_fastas "$out" "${protfiles[@]}"
        [[ -s "$out" ]] && rm -rf "$tmpd" && return 0
      fi
    else
      # maybe plain fasta
      if [[ -s "$zip" ]]; then
        format_fasta "$zip" "$out" && [[ -s "$out" ]] && rm -rf "$tmpd" && return 0
      fi
    fi
  fi

  # Fallback: FTP layout
  if [[ "$assembly" =~ ^(GC[AF])_[0-9]{3}[0-9]{3}[0-9]{3}\.[0-9]+$ ]]; then
    local prefix="${BASH_REMATCH[1]}"
    local digits="${assembly#${prefix}_}"           # e.g., 000/000/000 parts from accession
    local d1="${digits:0:3}" d2="${digits:3:3}" d3="${digits:6:3}"
    ftp_base="https://ftp.ncbi.nlm.nih.gov/genomes/all/${prefix}/${d1}/${d2}/${d3}"
    listing="$(curl -sSL "${ftp_base}/" || true)"
    asm_dir="$(printf "%s" "$listing" | grep -oP "${assembly}_[^/\"<>\s]+" | head -n1 || true)"
    if [[ -n "$asm_dir" ]]; then
      proturl="${ftp_base}/${asm_dir}/${asm_dir}_protein.faa.gz"
      if curl -sSLf "$proturl" -o "${tmpd}/prot.gz"; then
        gunzip -c "${tmpd}/prot.gz" > "${tmpd}/prot.faa" 2>/dev/null || true
        if [[ -s "${tmpd}/prot.faa" ]]; then
          format_fasta "${tmpd}/prot.faa" "$out" && [[ -s "$out" ]] && rm -rf "$tmpd" && return 0
        fi
      fi
    fi
  fi

  rm -rf "$tmpd"
  return 1
}

# -------- Uni* generic fetcher (KB / UniRef / UniParc) --------
fetch_uniprot_stream() {
  # fetch_uniprot_stream <endpoint> <query> <outfile>
  # endpoint: uniprotkb | uniref | uniparc
  local ep="$1" q="$2" out="$3"
  local qparam enc url
  enc="$(urlenc "$q")"
  qparam="(%28${enc}%29)"   # keep parentheses to group; tweak if you use fielded queries
  url="https://rest.uniprot.org/${ep}/stream?compressed=false&download=true&format=fasta&query=${qparam}"
  curl -sSL "$url" -o "${out}.raw" || return 1
  [[ -s "${out}.raw" ]] || return 1
  format_fasta "${out}.raw" "${out}" || return 1
  grep -q '^>' "${out}" || return 1
  return 0
}

# -------- CLI --------
while getopts ":q:o:m:" opt; do
  case "${opt}" in
    q) QUERY_FILE="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    m) MIN_COUNT="$OPTARG" ;;
    *) usage ;;
  esac
done

[[ -z "$QUERY_FILE" ]] && usage
[[ -f "$QUERY_FILE" ]] || die "Query file not found: $QUERY_FILE"
mkdir -p "$OUTPUT_DIR" || die "Cannot create output directory: $OUTPUT_DIR"

MISSING_LOG="missing.log"
{
  echo "=== Missing Queries Log ==="
  echo "Date: $(date --iso-8601=seconds 2>/dev/null || date)"
  echo "Query file: $QUERY_FILE"
  echo "Output directory: $OUTPUT_DIR"
  echo "Min count: $MIN_COUNT"
  echo "==========================="
  echo
} >"$MISSING_LOG"

# -------- Main loop --------
total=0; success=0; failed=0

while IFS= read -r raw || [[ -n "$raw" ]]; do
  # trim & skip comments/blanks
  query="$(printf "%s" "$raw" | normalize_spaces | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
  [[ -z "$query" || "$query" =~ ^# ]] && continue

  ((total++))
  safe_name="$(printf "%s" "$query" | sed 's/[^a-zA-Z0-9]/_/g')"
  out="${OUTPUT_DIR}/${safe_name}.fasta"
  log "=== Query: ${query} ==="

  qtmp="$(mktemp -d "${OUTPUT_DIR}/tmp_${safe_name}_XXXXXX")" || { log "ERROR: cannot create temp dir in $OUTPUT_DIR"; echo "$query" >>"$MISSING_LOG"; ((failed++)); continue; }
  gotfile=""

  # --- 1) Try NCBI route ---
  taxid="$(get_taxon_id "$query" || true)"
  if [[ -n "$taxid" ]]; then
    assembly="$(get_assembly_accession "$taxid" || true)"
    if [[ -n "$assembly" ]]; then
      if download_proteins_from_assembly "$assembly" "${qtmp}/ncbi.fasta" "$safe_name"; then
        if [[ -s "${qtmp}/ncbi.fasta" ]]; then
          format_fasta "${qtmp}/ncbi.fasta" "${qtmp}/ncbi.formatted" || true
          gotfile="${qtmp}/ncbi.formatted"
          log "NCBI: got proteins for ${query}"
        else
          log "NCBI: empty proteins for ${query}"
        fi
      else
        log "NCBI: failed to retrieve proteins via assembly ${assembly}"
      fi
    else
      log "NCBI: no assembly for taxid ${taxid}"
    fi
  else
    log "NCBI: taxid not found for ${query}"
  fi

  # Decide if Uni* fallback needed
  need_uniprot=true
  if [[ -n "$gotfile" && -s "$gotfile" ]]; then
    cnt="$(count_seqs "$gotfile")"
    [[ "$cnt" -ge 10 ]] && need_uniprot=false
  fi

  # --- 2) UniProtKB / UniRef / UniParc cascade ---
  uniprot_ok=false
  uniref_ok=false
  uniparc_ok=false

  if $need_uniprot; then
    log "Trying UniProtKB..."
    if fetch_uniprot_stream "uniprotkb" "$query" "${qtmp}/uniprot.fasta"; then
      if [[ -n "$gotfile" ]]; then
        merge_fastas "${qtmp}/merged.fasta" "$gotfile" "${qtmp}/uniprot.fasta"
        gotfile="${qtmp}/merged.fasta"
      else
        gotfile="${qtmp}/uniprot.fasta"
      fi
      uniprot_ok=true
      log "UniProtKB: retrieved for ${query}"
    else
      log "UniProtKB: no results for ${query}"
    fi

    # UniRef if UniProtKB failed or still thin
    if ! $uniprot_ok || { [[ -n "$gotfile" ]] && [[ "$(count_seqs "$gotfile")" -lt 10 ]]; }; then
      log "Trying UniRef..."
      if fetch_uniprot_stream "uniref" "$query" "${qtmp}/uniref.fasta"; then
        merge_fastas "${qtmp}/merged.fasta" ${gotfile:+$gotfile} "${qtmp}/uniref.fasta"
        gotfile="${qtmp}/merged.fasta"
        uniref_ok=true
        log "UniRef: merged results"
      else
        log "UniRef: no results"
      fi
    fi

    # UniParc if still empty/thin
    if ! $uniprot_ok && ! $uniref_ok || { [[ -n "$gotfile" ]] && [[ "$(count_seqs "$gotfile")" -lt 10 ]]; }; then
      log "Trying UniParc..."
      if fetch_uniprot_stream "uniparc" "$query" "${qtmp}/uniparc.fasta"; then
        merge_fastas "${qtmp}/merged.fasta" ${gotfile:+$gotfile} "${qtmp}/uniparc.fasta"
        gotfile="${qtmp}/merged.fasta"
        uniparc_ok=true
        log "UniParc: merged results"
      else
        log "UniParc: no results"
      fi
    fi
  else
    log "Skipping Uniprot (NCBI sufficient)"
  fi

  # --- 3) Finalize: format + dedup + move ---
  if [[ -n "$gotfile" && -s "$gotfile" && "$(count_seqs "$gotfile")" -gt 0 ]]; then
    format_fasta "$gotfile" "${qtmp}/final.fasta" || true
    dedup_by_header "${qtmp}/final.fasta" "${qtmp}/final.dedup"
    if [[ -s "${qtmp}/final.dedup" ]]; then
      mv "${qtmp}/final.dedup" "$out"
      cnt="$(count_seqs "$out")"
      ((success++))
      log "Saved: ${out} (${cnt} sequences)"
      if [[ "$cnt" -lt "$MIN_COUNT" ]]; then
        msg="WARNING: ${out} contains only ${cnt} proteins (minimum required: ${MIN_COUNT})"
        log "$msg"
        echo "$msg" >>"$MISSING_LOG"
      fi
    else
      log "ERROR: final FASTA empty after processing for ${query}"
      echo "$query" >>"$MISSING_LOG"
      ((failed++))
    fi
  else
    log "ERROR: No sequences found for ${query}"
    echo "$query" >>"$MISSING_LOG"
    ((failed++))
  fi

  rm -rf "$qtmp" || true
  log ""
done <"$QUERY_FILE"

# -------- Summary --------
{
  echo
  echo "=== Summary ==="
  echo "Total queries processed: $total"
  echo "Successful queries:      $success"
  echo "Failed queries:          $failed"
  if [[ $total -gt 0 ]]; then
    echo "Success rate: $((success * 100 / total))%"
  else
    echo "Success rate:            0%"
  fi
} >>"$MISSING_LOG"

echo "Total: $total, Success: $success, Failed: $failed"
echo "Results: $OUTPUT_DIR"
echo "Missing logged: $MISSING_LOG"

