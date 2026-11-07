#   bash 03_trim_quant.sh [fastq_dir] [salmon_index] [threads]
#   DELETE_RAW=1 bash 03_trim_quant.sh ...   # also delete raw FASTQ after a sample's quant succeeds
#
# Paired-end only (_1/_2.fastq.gz). Resumable: samples with quant.sf are skipped.

set -euo pipefail

INDIR="${1:-fastq}"
INDEX="${2:-ref/salmon_index_mm10}"
THREADS="${3:-8}"
DELETE_RAW="${DELETE_RAW:-0}"

TMPTRIM="trimmed_tmp"
QCDIR="qc"
QUANT="results/salmon_quant"
LOG="$QCDIR/trim_quant.log"

for tool in fastqc fastp salmon; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: '$tool' not found"; exit 1; }
done
[[ -d "$INDIR" ]]  || { echo "ERROR: $INDIR not found"; exit 1; }
[[ -f "$INDEX/versionInfo.json" || -f "$INDEX/info.json" ]] || { echo "ERROR: Salmon index not found at $INDEX"; exit 1; }

mkdir -p "$TMPTRIM" "$QCDIR/raw" "$QCDIR/trimmed" "$QCDIR/fastp" "$QUANT"

mapfile -t SAMPLES < <(find "$INDIR" -maxdepth 1 -name "*_1.fastq.gz" -printf "%f\n" \
                       | sed 's/_1\.fastq\.gz$//' | sort)
echo "Found ${#SAMPLES[@]} paired-end samples (DELETE_RAW=$DELETE_RAW)" | tee "$LOG"

FAILED=()
for s in "${SAMPLES[@]}"; do
    R1="$INDIR/${s}_1.fastq.gz"; R2="$INDIR/${s}_2.fastq.gz"
    T1="$TMPTRIM/${s}_trimmed_1.fastq.gz"; T2="$TMPTRIM/${s}_trimmed_2.fastq.gz"

    if [[ -f "$QUANT/$s/quant.sf" ]]; then
        echo "=== $s: quant.sf exists, skipping ===" | tee -a "$LOG"; continue
    fi
    [[ -f "$R2" ]] || { echo "=== $s: missing _2 file ===" | tee -a "$LOG"; FAILED+=("$s"); continue; }

    echo "=== $s ($(date +%H:%M)) ===" | tee -a "$LOG"
    {
        fastqc -t "$THREADS" -o "$QCDIR/raw" "$R1" "$R2" &&
        fastp -i "$R1" -I "$R2" -o "$T1" -O "$T2" --detect_adapter_for_pe \
              --thread "$((THREADS > 16 ? 16 : THREADS))" \
              --qualified_quality_phred 20 --length_required 36 \
              --html "$QCDIR/fastp/${s}.html" --json "$QCDIR/fastp/${s}.json" &&
        fastqc -t "$THREADS" -o "$QCDIR/trimmed" "$T1" "$T2" &&
        salmon quant -i "$INDEX" -l A -1 "$T1" -2 "$T2" \
              --validateMappings --gcBias --seqBias \
              -p "$THREADS" -o "$QUANT/$s"
    } >>"$LOG" 2>&1 || { echo "  FAILED" | tee -a "$LOG"; FAILED+=("$s"); rm -f "$T1" "$T2"; continue; }

    rm -f "$T1" "$T2"           # free trimmed copy right away
    if [[ "$DELETE_RAW" == "1" && -f "$QUANT/$s/quant.sf" ]]; then
        rm -f "$R1" "$R2"
        echo "  done (raw FASTQ deleted)" | tee -a "$LOG"
    else
        echo "  done" | tee -a "$LOG"
    fi
done

rmdir "$TMPTRIM" 2>/dev/null || true

if command -v multiqc >/dev/null 2>&1; then
    multiqc "$QCDIR" "$QUANT" -o "$QCDIR/multiqc" -f >>"$LOG" 2>&1 || true
fi

echo; echo "Quant results: $QUANT   QC: $QCDIR/multiqc/multiqc_report.html"
[[ ${#FAILED[@]} -eq 0 ]] || { echo "FAILED: ${FAILED[*]}"; exit 1; }
