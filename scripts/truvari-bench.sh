#!/bin/bash
# truvari-bench.sh
# This script benchmarks CVCF (calls) against BVCF (truth) by SV type.
# It computes sensitivity (recall) and precision for deletions, duplications, and insertions.
#
# Usage:
#   ./truvari-bench.sh -b <BVCF> -c <CVCF> -w <WORKDIR> -r <REFERENCE> [-d <B_DUP_IDS>] [-D <C_DUP_IDS>]
#
#   -b : Path to the truth VCF (BVCF)
#   -c : Path to the calls VCF (CVCF)
#   -w : Working directory for intermediate files and results
#   -r : Reference FASTA file (e.g. GRCh38)
#   -p : Path to truvari executable
#   -d : (Optional) File with truth IDs treated as DUP (BVCF-side reclassification)
#   -D : (Optional) File with call IDs treated as DUP (CVCF-side reclassification)
#
# Example:
#   ./truvari-bench.sh -b truth.vcf.gz -c calls.vcf.gz \
#      -w truvari-bench -r ~/references/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa \
#      -d dup_ids.truth.txt -D dup_ids.calls.txt

set -euo pipefail

sv_sensitivity() {
  local tp=$(jq -r '."TP-base"' "$1")
  local fn=$(jq -r '.FN' "$1")
  local N=$tp M=$((tp+fn))
  echo "${2:-DEL} SENSITIVITY: $N/$M = $(awk -v n=$N -v m=$M 'BEGIN{printf "%.2f",n/m}')"
}

sv_precision() {
  local tp=$(jq -r '."TP-comp"' "$1")
  local fp=$(jq -r '.FP' "$1")
  local N=$tp M=$((tp+fp))
  echo "${2:-DEL} PRECISION: $N/$M = $(awk -v n=$N -v m=$M 'BEGIN{printf "%.2f",n/m}')"
}

# Parse command-line arguments.
B_DUP_IDS=""
C_DUP_IDS=""
while getopts "b:c:w:r:p:d:D:" opt; do
    case $opt in
        b) BVCF="$OPTARG" ;;
        c) CVCF="$OPTARG" ;;
        w) WORKDIR="$OPTARG" ;;
        r) REFERENCE="$OPTARG" ;;
        p) TRUVARI="$OPTARG" ;;
        d) B_DUP_IDS="$OPTARG" ;;
        D) C_DUP_IDS="$OPTARG" ;;
        *) echo "Usage: $0 -b <BVCF> -c <CVCF> -w <WORKDIR> -r <REFERENCE> [-d <B_DUP_IDS>]" && exit 1 ;;
    esac
done

if [[ -z "${BVCF:-}" || -z "${CVCF:-}" || -z "${WORKDIR:-}" || -z "${REFERENCE:-}" || -z "${TRUVARI:-}" ]]; then
    echo "Missing required arguments."
    echo "Usage: $0 -b <BVCF> -c <CVCF> -w <WORKDIR> -r <REFERENCE> -p <TRUVARI_PATH> [-d <B_DUP_IDS>]"
    exit 1
fi

echo "Benchmarking:"
echo "Truth VCF: $BVCF"
echo "Call VCF:  $CVCF"
echo "Workdir:   $WORKDIR"
echo "Reference: $REFERENCE"
echo "Truvari path:   $TRUVARI"
if [[ -n "$B_DUP_IDS" ]]; then
    echo "Duplication IDs file: $B_DUP_IDS"
fi
if [[ -n "$C_DUP_IDS" ]]; then
    echo "Calls DUP ID file (C_DUP_IDS): $C_DUP_IDS"
fi
echo ""

# Create directories for intermediate VCFs and results.
mkdir -p "$WORKDIR/vcfs"
mkdir -p "$WORKDIR/results"

# Helper function: filter, compress (bgzip), and index a VCF.
# Arguments:
#   $1: Filter expression (for bcftools view -i)
#   $2: Input VCF file
#   $3: Output file (should end with .vcf.gz)
filter_and_index() {
    bcftools view -i "$1" "$2" -O z -o "$3"
    tabix -p vcf "$3"
}

#######################################
# 1. Benchmark Deletions
#######################################
echo "Filtering deletions..."
filter_and_index 'INFO/SVTYPE=="DEL"' "$BVCF" "$WORKDIR/vcfs/BVCF_DEL.vcf.gz"
filter_and_index 'INFO/SVTYPE=="DEL"' "$CVCF" "$WORKDIR/vcfs/CVCF_DEL.vcf.gz"

# Remove preexisting output directory if it exists.
if [ -d "$WORKDIR/results/del" ]; then rm -rf "$WORKDIR/results/del"; fi

echo "Running truvari for deletions..."
"$TRUVARI" bench \
    -b "$WORKDIR/vcfs/BVCF_DEL.vcf.gz" \
    -c "$WORKDIR/vcfs/CVCF_DEL.vcf.gz" \
    -o "$WORKDIR/results/del" \
    --reference "$REFERENCE" \
    --no-ref a --pick multi

#######################################
# 2. Prepare INS/DUP VCFs for Truth and Calls
#######################################
if [[ -n "$B_DUP_IDS" ]]; then
    echo "Reclassifying truth VCF records based on duplication ID file: $B_DUP_IDS"
    
    echo "Creating truth duplications (BVCF_DUP): records with SVTYPE==\"DUP\" OR ID in file"
    bcftools view -i 'INFO/SVTYPE=="DUP" || ID==@'"$B_DUP_IDS" "$BVCF" -O z -o "$WORKDIR/vcfs/BVCF_DUP.vcf.gz"
    tabix -p vcf "$WORKDIR/vcfs/BVCF_DUP.vcf.gz"
    
    echo "Creating truth insertions (BVCF_INS): records with SVTYPE==\"INS\" AND ID not in file"
    bcftools view -i 'INFO/SVTYPE=="INS" && ID!=@'"$B_DUP_IDS" "$BVCF" -O z -o "$WORKDIR/vcfs/BVCF_INS.temp.vcf.gz"
    python3 ~/surveyor_paper/scripts/add_svlen.py "$WORKDIR/vcfs/BVCF_INS.temp.vcf.gz" "$WORKDIR/vcfs/BVCF_INS.vcf.gz"
    tabix -p vcf "$WORKDIR/vcfs/BVCF_INS.vcf.gz"
    
    echo "Creating merged truth INS/DUP VCF (BVCF_INS_DUP)..."
    bcftools concat -a "$WORKDIR/vcfs/BVCF_DUP.vcf.gz" "$WORKDIR/vcfs/BVCF_INS.vcf.gz" -O z -o "$WORKDIR/vcfs/BVCF_INS_DUP.vcf.gz"
    tabix -p vcf "$WORKDIR/vcfs/BVCF_INS_DUP.vcf.gz"
else
    echo "No duplication IDs file provided. Using standard filtering for truth INS/DUP VCFs."
    filter_and_index 'INFO/SVTYPE=="DUP"' "$BVCF" "$WORKDIR/vcfs/BVCF_DUP.vcf.gz"
    filter_and_index 'INFO/SVTYPE=="INS"' "$BVCF" "$WORKDIR/vcfs/BVCF_INS.temp.vcf.gz"
    python3 ~/surveyor_paper/scripts/add_svlen.py "$WORKDIR/vcfs/BVCF_INS.temp.vcf.gz" "$WORKDIR/vcfs/BVCF_INS.vcf.gz"
    tabix -p vcf "$WORKDIR/vcfs/BVCF_INS.vcf.gz"

    echo "Creating merged truth INS/DUP VCF (BVCF_INS_DUP)..."
    bcftools concat -a "$WORKDIR/vcfs/BVCF_DUP.vcf.gz" "$WORKDIR/vcfs/BVCF_INS.vcf.gz" -O z -o "$WORKDIR/vcfs/BVCF_INS_DUP.vcf.gz"
    tabix -p vcf "$WORKDIR/vcfs/BVCF_INS_DUP.vcf.gz"
fi

if [[ -n "$C_DUP_IDS" ]]; then
    echo "Reclassifying calls VCF records based on duplication ID file (C_DUP_IDS): $C_DUP_IDS"
    echo "Creating calls duplications (CVCF_DUP): SVTYPE==DUP OR ID in file"
    bcftools view -i 'INFO/SVTYPE=="DUP" || ID==@'"$C_DUP_IDS" "$CVCF" -O z -o "$WORKDIR/vcfs/CVCF_DUP.vcf.gz"
    tabix -p vcf "$WORKDIR/vcfs/CVCF_DUP.vcf.gz"

    echo "Creating calls insertions (CVCF_INS): SVTYPE==INS AND ID not in file"
    bcftools view -i 'INFO/SVTYPE=="INS" && ID!=@'"$C_DUP_IDS" "$CVCF" -O z -o "$WORKDIR/vcfs/CVCF_INS.temp.vcf.gz"
    python3 ~/surveyor_paper/scripts/add_svlen.py "$WORKDIR/vcfs/CVCF_INS.temp.vcf.gz" "$WORKDIR/vcfs/CVCF_INS.vcf.gz"
    tabix -p vcf "$WORKDIR/vcfs/CVCF_INS.vcf.gz"

    echo "Creating merged calls INS/DUP VCF (CVCF_INS_DUP)..."
    bcftools concat -a "$WORKDIR/vcfs/CVCF_DUP.vcf.gz" "$WORKDIR/vcfs/CVCF_INS.vcf.gz" -O z -o "$WORKDIR/vcfs/CVCF_INS_DUP.vcf.gz"
    tabix -p vcf "$WORKDIR/vcfs/CVCF_INS_DUP.vcf.gz"
else
    echo "No C_DUP_IDS provided. Using standard filtering for calls INS/DUP VCFs."
    filter_and_index 'INFO/SVTYPE=="DUP"' "$CVCF" "$WORKDIR/vcfs/CVCF_DUP.vcf.gz"
    filter_and_index 'INFO/SVTYPE=="INS"' "$CVCF" "$WORKDIR/vcfs/CVCF_INS.temp.vcf.gz"
    python3 ~/surveyor_paper/scripts/add_svlen.py "$WORKDIR/vcfs/CVCF_INS.temp.vcf.gz" "$WORKDIR/vcfs/CVCF_INS.vcf.gz"
    tabix -p vcf "$WORKDIR/vcfs/CVCF_INS.vcf.gz"

    echo "Creating merged calls INS/DUP VCF (CVCF_INS_DUP)..."
    bcftools concat -a "$WORKDIR/vcfs/CVCF_DUP.vcf.gz" "$WORKDIR/vcfs/CVCF_INS.vcf.gz" -O z -o "$WORKDIR/vcfs/CVCF_INS_DUP.vcf.gz"
    tabix -p vcf "$WORKDIR/vcfs/CVCF_INS_DUP.vcf.gz"
fi

#######################################
# 3. Benchmark Duplications
#######################################
# Duplication Sensitivity:
#   Truth: Only duplications from BVCF.
#   Call: Merged INS+DUP from CVCF.
echo "Benchmarking duplication sensitivity..."
if [ -d "$WORKDIR/results/dup_sens" ]; then rm -rf "$WORKDIR/results/dup_sens"; fi
"$TRUVARI" bench \
    -b "$WORKDIR/vcfs/BVCF_DUP.vcf.gz" \
    -c "$WORKDIR/vcfs/CVCF_INS_DUP.vcf.gz" \
    -o "$WORKDIR/results/dup_sens" \
    --reference "$REFERENCE" \
    --dup-to-ins --no-ref a --pick multi -p 0 -P 0

# Duplication Precision:
#   Truth: Merged INS/DUP from BVCF.
#   Call: Only duplications from CVCF.
echo "Benchmarking duplication precision..."
if [ -d "$WORKDIR/results/dup_prec" ]; then rm -rf "$WORKDIR/results/dup_prec"; fi
"$TRUVARI" bench \
    -b "$WORKDIR/vcfs/BVCF_INS_DUP.vcf.gz" \
    -c "$WORKDIR/vcfs/CVCF_DUP.vcf.gz" \
    -o "$WORKDIR/results/dup_prec" \
    --reference "$REFERENCE" \
    --dup-to-ins --no-ref a --pick multi -p 0 -P 0

#######################################
# 4. Benchmark Insertions
#######################################
# Insertion Sensitivity:
#   Truth: Only insertions from BVCF.
#   Call: Merged INS+DUP from CVCF.
echo "Benchmarking insertion sensitivity..."
if [ -d "$WORKDIR/results/ins_sens" ]; then rm -rf "$WORKDIR/results/ins_sens"; fi
"$TRUVARI" bench \
    -b "$WORKDIR/vcfs/BVCF_INS.vcf.gz" \
    -c "$WORKDIR/vcfs/CVCF_INS_DUP.vcf.gz" \
    -o "$WORKDIR/results/ins_sens" \
    --reference "$REFERENCE" \
    --dup-to-ins --no-ref a --pick multi -p 0 -P 0

# Insertion Precision:
#   Truth: Merged INS/DUP from BVCF.
#   Call: Only insertions from CVCF.
echo "Benchmarking insertion precision..."
if [ -d "$WORKDIR/results/ins_prec" ]; then rm -rf "$WORKDIR/results/ins_prec"; fi
"$TRUVARI" bench \
    -b "$WORKDIR/vcfs/BVCF_INS_DUP.vcf.gz" \
    -c "$WORKDIR/vcfs/CVCF_INS.vcf.gz" \
    -o "$WORKDIR/results/ins_prec" \
    --reference "$REFERENCE" \
    --dup-to-ins --no-ref a --pick multi -p 0 -P 0

#######################################
# 5. Summarize Results
#######################################
echo ""
echo "Benchmarking complete. Summary metrics (recall and precision) are reported below:"
echo "--------------------------------------------------------------------------------------------"
{
    sv_sensitivity "$WORKDIR/results/del/summary.json" DEL
    sv_precision "$WORKDIR/results/del/summary.json" DEL
    sv_sensitivity "$WORKDIR/results/dup_sens/summary.json" DUP
    sv_precision "$WORKDIR/results/dup_prec/summary.json" DUP
    sv_sensitivity "$WORKDIR/results/ins_sens/summary.json" INS
    sv_precision "$WORKDIR/results/ins_prec/summary.json" INS
} > "$WORKDIR/results/summary.txt"
cat "$WORKDIR/results/summary.txt"
echo "All intermediate VCFs and results can be found in $WORKDIR."
