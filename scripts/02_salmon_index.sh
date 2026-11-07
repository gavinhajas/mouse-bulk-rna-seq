#   bash 02_salmon_index.sh [ref_dir] [threads]
#
# expects in ref_dir:
#   Mus_musculus.GRCm38.dna.primary_assembly.fa
#   Mus_musculus.GRCm38.102.gtf
# builds (skipping anything that already exists):
#   Mus_musculus.GRCm38.transcripts.fa, decoys.txt, gentrome.fa, salmon_index_mm10/, tx2gene.tsv

set -euo pipefail

REF="${1:-ref}"
THREADS="${2:-8}"

GENOME="$REF/Mus_musculus.GRCm38.dna.primary_assembly.fa"
GTF="$REF/Mus_musculus.GRCm38.102.gtf"
TRANSCRIPTS="$REF/Mus_musculus.GRCm38.transcripts.fa"
DECOYS="$REF/decoys.txt"
GENTROME="$REF/gentrome.fa"
INDEX="$REF/salmon_index_mm10"
TX2GENE="$REF/tx2gene.tsv"

command -v salmon >/dev/null 2>&1 || { echo "ERROR: salmon not found (conda install -c bioconda salmon)"; exit 1; }
[[ -f "$GENOME" && -f "$GTF" ]] || { echo "ERROR: genome FASTA or GTF missing in $REF"; exit 1; }

# 1. transcript sequences extracted from genome + GTF
if [[ ! -s "$TRANSCRIPTS" ]]; then
    command -v gffread >/dev/null 2>&1 || { echo "ERROR: gffread not found (conda install -c bioconda gffread)"; exit 1; }
    echo ">>> Extracting transcripts with gffread"
    gffread -w "$TRANSCRIPTS" -g "$GENOME" "$GTF"
fi

# 2. decoys = names of all genome sequences
if [[ ! -s "$DECOYS" ]]; then
    echo ">>> Writing decoys.txt"
    grep '^>' "$GENOME" | cut -d ' ' -f 1 | sed 's/^>//' > "$DECOYS"
fi

# 3. gentrome = transcripts FIRST, then genome
if [[ ! -s "$GENTROME" ]]; then
    echo ">>> Building gentrome.fa"
    cat "$TRANSCRIPTS" "$GENOME" > "$GENTROME"
fi

# 4. salmon index
if [[ -f "$INDEX/versionInfo.json" && -f "$INDEX/pos.bin" ]] || [[ -f "$INDEX/info.json" && -f "$INDEX/refseq.bin" ]]; then
    echo ">>> Salmon index already exists at $INDEX, skipping"
else
    echo ">>> Building Salmon index (this can take a while)"
    rm -rf "$INDEX"
    salmon index -t "$GENTROME" -d "$DECOYS" -p "$THREADS" -i "$INDEX" -k 31
fi

# 5. tx2gene table for tximport (transcript_id, gene_id, gene_name)
if [[ ! -s "$TX2GENE" ]]; then
    echo ">>> Writing tx2gene.tsv"
    awk -F'\t' '$3=="transcript" {
        tid=""; gid=""; gn="";
        if (match($9, /transcript_id "[^"]+"/)) tid=substr($9, RSTART+15, RLENGTH-16);
        if (match($9, /gene_id "[^"]+"/))       gid=substr($9, RSTART+9,  RLENGTH-10);
        if (match($9, /gene_name "[^"]+"/))     gn =substr($9, RSTART+11, RLENGTH-12);
        print tid "\t" gid "\t" (gn=="" ? gid : gn)
    }' "$GTF" > "$TX2GENE"
fi

echo
echo "Done. Index: $INDEX"
echo "Transcripts: $(grep -c '^>' "$TRANSCRIPTS")   Decoys: $(wc -l < "$DECOYS")   tx2gene rows: $(wc -l < "$TX2GENE")"
