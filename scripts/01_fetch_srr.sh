#   1. CREATE a text file with 11 SRR IDs in accessions.txt
#   2. bash fetch_srr.sh [accessions.txt] [output_dir] [threads]

set -euo pipefail

ACC_FILE="${1:-accessions.txt}"
OUTDIR="${2:-fastq}"
THREADS="${3:-8}"
SRA_DIR="${OUTDIR}/sra_cache"
LOG="${OUTDIR}/download.log"

#checks
for tool in prefetch fasterq-dump; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: '$tool' not found. Install sra-tools:"
        echo "  conda install -c bioconda sra-tools"
        exit 1
    }
done

[[ -f "$ACC_FILE" ]] || {
    echo "ERROR: $ACC_FILE not found. Create it with one SRR ID per line."
    exit 1
}

mkdir -p "$OUTDIR" "$SRA_DIR"

#read accessions, skipping blanks/comments, trimming whitespace
mapfile -t ACCS < <(grep -v '^\s*#' "$ACC_FILE" | grep -v '^\s*$' | awk '{print $1}')
echo "Found ${#ACCS[@]} accessions in $ACC_FILE"
[[ ${#ACCS[@]} -eq 11 ]] || echo "WARNING: expected 11 accessions, got ${#ACCS[@]}"

#download loop
FAILED=()
for srr in "${ACCS[@]}"; do
    echo "=== $srr ===" | tee -a "$LOG"

    #skip if already done
    if ls "${OUTDIR}/${srr}"*.fastq* >/dev/null 2>&1; then
        echo "  already downloaded, skipping" | tee -a "$LOG"
        continue
    fi

    #download .sra file + retry
    if ! prefetch "$srr" --output-directory "$SRA_DIR" --max-size 100G >>"$LOG" 2>&1; then
        echo "  prefetch FAILED" | tee -a "$LOG"
        FAILED+=("$srr"); continue
    fi

    #convert to FASTQ; --split-files gives _1/_2 for paired-end
    if ! fasterq-dump "${SRA_DIR}/${srr}/${srr}.sra" \
            --split-files --threads "$THREADS" \
            --outdir "$OUTDIR" --temp "$OUTDIR" >>"$LOG" 2>&1; then
        echo "  fasterq-dump FAILED" | tee -a "$LOG"
        FAILED+=("$srr"); continue
    fi

    #compress, we use pigs because it is faster than gzip and supports multithreading
    if command -v pigz >/dev/null 2>&1; then
        pigz -p "$THREADS" "${OUTDIR}/${srr}"*.fastq
    else
        gzip "${OUTDIR}/${srr}"*.fastq
    fi

    #free up space by removing the .sra file
    rm -rf "${SRA_DIR:?}/${srr}"
    echo "  done" | tee -a "$LOG"
done

#summary
echo
echo "Finished. FASTQ files in: $OUTDIR"
ls -lh "$OUTDIR"/*.fastq.gz 2>/dev/null || true

if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo "FAILED accessions (rerun the script to retry): ${FAILED[*]}"
    exit 1
fi
