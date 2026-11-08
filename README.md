[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24049830.svg)](https://doi.org/10.5281/zenodo.24049830)
# RNA-seq Differential Expression: WT vs Brca1+/- (HFD mouse liver)

Bulk RNA-seq analysis comparing liver from high-fat-diet (HFD) fed wild-type mice (control, n = 6) against Brca1+/- mice (experimental, n = 6). Reads are quantified with Salmon against the mouse Ensembl transcriptome and tested for differential expression with DESeq2.

**Organism:** *Mus musculus*

**Dataset:** GEO samples GSM9573936-GSM9573947 (series accession: GSE324339), SRA runs SRR37536766-SRR37536777

**Layout:** paired-end

**Mice:** female WT and Brca1+/- mice fed a 60% high-fat diet (D12492i) for 22 weeks starting at 4 weeks of age; liver RNA-seq, n = 6 per group (GEO series description; see Citation).

**Reference:** Ensembl release 102, GRCm38 (mm10)

## **Pipeline**

Study design -> Fetch SRR -> Download/split FASTQ -> FastQC + fastp -> Salmon quant -> tximport -> DESeq2

Reference transcriptome -> Salmon index ----^

| Step | Tool | Script | Output |
| ----- | ----- | ----- | ----- |
| 1. Study design | - | metadata/samples.csv | sample/condition table |
| 2-3. Fetch & split FASTQ | sra-tools | scripts/01_fetch_srr.sh | fastq/*.fastq.gz |
| 5. Salmon index | Salmon | scripts/02_salmon_index.sh | ref/salmon_index_mm10/, ref/tx2gene.tsv |
| 4+6. QC, trim, quantify | FastQC, fastp, Salmon | scripts/03_trim_quant.sh | qc/, results/salmon_quant/ |
| 7-8. Import & DE | tximport, DESeq2 | scripts/04_deseq2.R | results/deseq2/ |
| 9. Pathway analysis (exploratory) | fgsea, msigdbr | scripts/05_enrichment.R | results/enrichment/ |

Steps 4 and 6 run together, one sample at a time, because the machine has limited disk space: each sample's trimmed reads are deleted as soon as its quantification finishes.

## **Samples**

| Group | SRR | GSM | Name |
| ----- | ----- | ----- | ----- |
| control | SRR37536777 | GSM9573936 | WT_HFD_1 |
| control | SRR37536776 | GSM9573937 | WT_HFD_2 |
| control | SRR37536775 | GSM9573938 | WT_HFD_3 |
| control | SRR37536774 | GSM9573939 | WT_HFD_4 |
| control | SRR37536773 | GSM9573940 | WT_HFD_5 |
| control | SRR37536772 | GSM9573941 | WT_HFD_6 |
| experimental | SRR37536771 | GSM9573942 | Brca1het_HFD_1 |
| experimental | SRR37536770 | GSM9573943 | Brca1het_HFD_2 |
| experimental | SRR37536769 | GSM9573944 | Brca1het_HFD_3 |
| experimental | SRR37536768 | GSM9573945 | Brca1het_HFD_4 |
| experimental | SRR37536767 | GSM9573946 | Brca1het_HFD_5 |
| experimental | SRR37536766 | GSM9573947 | Brca1het_HFD_6 |

(Machine-readable version: metadata/samples.csv. The SRR-to-GSM pairing was checked on GEO for rows 6 (GSM9573941 = SRR37536772) and 12 (GSM9573947 = SRR37536766); the other rows follow the same order.)

## **Quantification QC summary**

All 12 samples were quantified with the same pipeline and options (-l A, same index, Salmon 2.8.0; confirmed from each sample's aux_info/cmd_info.json).

**Library type:** all 12 samples were detected as IU (unstranded, paired-end, inward-facing), with no library-type differences between samples.

**Mapping rate:** 84.7% - 89.3% across samples (range of 4.6 percentage points, no outliers). Lowest: SRR37536773 (84.7%), SRR37536771 (85.5%). Rep 6 (SRR37536766) maps at 88.1%, in line with the other samples.

**MultiQC:** qc/multiqc/multiqc_report.html (general statistics table reviewed; per-base quality, adapter and GC plots not yet reviewed). Raw depth: 19.9-29.1 million read pairs per sample (lowest: SRR37536766, rep 6, 19.9 M; highest: SRR37536773, 29.1 M). fastp: 98.7-99.4% of reads passed filtering, 1.5-1.9% adapter-trimmed. FastQC: GC 48-49% in all samples, 66-70% duplicate reads (not deduplicated; common in RNA-seq of liver, where a few genes dominate). Salmon mapping 84.7-89.3%. The QC metrics are similar across samples.

## **Requirements**

Conda environment bioinfo. Record exact versions with:

conda list | grep -E "sra-tools|fastqc|fastp|multiqc|salmon|gffread"

conda env export --from-history > environment.yml

| Tool | Version |
| ----- | ----- |
| Salmon | 2.8.0 |
| sra-tools | 3.0.0 (manual install in ~/bin/sratoolkit.3.0.0-ubuntu64, not in the conda env) |
| FastQC | 0.12.1 |
| fastp | 1.3.6 |
| MultiQC | 1.35 |
| gffread | 0.12.9 |
| R | 4.5.2 |
| tximport | 1.38.2 |
| DESeq2 | 1.50.2 |
| apeglm | 1.32.0 |
| fgsea | 1.36.2 |
| msigdbr | 26.1.1 |

Salmon 2.x is a Rust rewrite of the earlier C++ Salmon. Its index format changed in 2.7.0, so indexes built with Salmon 2.6 or earlier (including all 1.x) cannot be loaded and must be rebuilt.

## **Usage**

Run everything from the project root (~/bioinfo), in order:

bash scripts/01_fetch_srr.sh metadata/accessions.txt fastq 8

bash scripts/02_salmon_index.sh ref 8

bash scripts/03_trim_quant.sh fastq ref/salmon_index_mm10 8

Rscript scripts/04_deseq2.R

# optional sensitivity run without one sample (writes to a separate folder):

EXCLUDE=Brca1het_HFD_5 Rscript scripts/04_deseq2.R

# exploratory pathway analysis (GSEA) and sensitivity version:

Rscript scripts/05_enrichment.R

RES_FILE=results/deseq2_excl_Brca1het_HFD_5/deseq2_all_genes.csv OUT_DIR=results/enrichment_excl_Brca1het_HFD_5 Rscript scripts/05_enrichment.R

The scripts skip completed steps and can be rerun. DELETE_RAW=1 bash scripts/03_trim_quant.sh ... also deletes each raw FASTQ pair after its quantification succeeds (only needed if disk space runs out). 01_fetch_srr.sh prints a harmless "expected 11 accessions" warning (the check is hard-coded).

## **Directory structure**

bioinfo/

- README.md

- environment.yml (from-history) and environment.lock.yml (exact versions)

- metadata/ samples.csv, accessions.txt

- scripts/ numbered, run in order

- fastq/ raw reads (not tracked in git)

- ref/ genome, GTF, transcripts, decoys, Salmon index, tx2gene.tsv

- qc/ FastQC, fastp and MultiQC reports, trim_quant.log

- results/ salmon_quant/, deseq2/, enrichment/, old/ (superseded manual run)

- logs/

## **Methods / parameters**

**Reference:** Ensembl release 102 GRCm38 primary-assembly genome and GTF. Transcripts extracted with gffread.

**Index (Salmon 2.8.0):** decoy-aware index built from a gentrome (transcripts followed by genome) with the genome sequence names as decoys (66 decoys), k = 31. ref/tx2gene.tsv (transcript ID, gene ID, gene name) is built from the same GTF.

**Trimming (fastp):** adapters auto-detected for paired-end, minimum base quality Q20, minimum read length 36 bp.

**Quantification (Salmon 2.8.0):** automatic library-type detection (-l A; resolved to IU for all samples), 8 threads. --validateMappings --gcBias --seqBias are passed by scripts/03_trim_quant.sh. Each sample's aux_info/meta_info.json confirms seq_bias_correct: true and gc_bias_correct: true (positional bias correction off). --validateMappings is not recorded in meta_info.json or cmd_info.json.

**Pathway analysis (exploratory, scripts/05_enrichment.R):** gene set enrichment analysis with fgsea (fgseaMultilevel, seed 42, gene sets of 15-500 genes) on all tested genes ranked by the DESeq2 Wald statistic (one Ensembl gene per symbol, strongest signal kept). Gene sets: mouse MSigDB Hallmark (50 sets) and GO Biological Process (about 3,650 sets) via msigdbr. Positive NES = higher in Brca1+/-. GSEA permutes gene labels and treats genes as independent, so with this much between-mouse variation its p-values are overstated; NES, direction and stability across the sensitivity run are more informative.

**Differential expression (DESeq2, scripts/04_deseq2.R):** gene-level counts imported with tximport (transcripts summed to genes using ref/tx2gene.tsv, version suffixes ignored). Genes with fewer than 10 reads in at least 6 samples were removed. Design ~ condition (experimental vs control; positive log2FC = higher in Brca1+/-). Differentially expressed genes: Benjamini-Hochberg adjusted p < 0.05 and |log2FC| >= 1. Log2 fold changes in plots are shrunken (apeglm if installed, otherwise ashr/normal). PCA and sample-distance heatmap use VST-transformed counts. Brca1 expression is checked as a sanity control. A sensitivity analysis can drop samples via the EXCLUDE environment variable (see Usage); outputs go to a separate results/deseq2_excl_<sample>/ folder.

## **Results summary (DESeq2, all 12 samples)**

13,502 of 54,416 genes passed the low-count filter and were tested.

**Brca1** is reduced in Brca1+/- liver: log2FC -1.55 (apeglm-shrunken -1.43), padj 2.3e-4; all six mutants have lower normalised counts than all six controls.

Differential expression analysis identified 72 statistically significant genes (Genes with padj < 0.05; 43 upregulated, 29 downregulated). Filtering for at least a two-fold change in expression |log2FC| >= 1 yielded 14 genes: 13 upregulated (Dntt, Serpina12, Apoa4, Gm42047, Tubb2a, Arap3, Hhipl2, Fgf21, Gale, Gm47761, Fam131c, Gm47528, Muc1) and one downregulated gene (Brca1). Most of these transcripts fell near the lower boundary of the fold-change cutoff (|log2FC| about 1.0-1.25).

Sample-level quality control via PCA (PC1: 26% variance; PC2: 15% variance) and sample-to-sample distance clustering demonstrated incomplete segregation by genotype, largely obscured by high biological variability between individual mice. Qualitative examination of the top loadings showed that PC1 was enriched for hepatic xenobiotic metabolism (Cyp3a, Sult2a, Fmo3) and lipid homeostasis genes (Scd1, Apoa4, Plin4), while PC2 reflected stress- and fasting-related transcriptional programs (Egr1, Mt1, Mt2, Igfbp1, G6pc). Notably, replicate Brca1het_HFD_5 exhibited the greatest transcriptional distance from the remaining cohort.

While the broader study utilized both sexes, hepatic RNA-seq was restricted to female mice (n = 6 per group; Brca1+/- HFD vs. WT HFD), as documented in the GEO series metadata and published article (Figure 5A). Transcriptomic sex verification confirmed that all 12 samples were female, evidenced by high Xist abundance (9,600–20,300 raw Salmon reads per sample) and zero reads aligned to Y-chromosome genes (Ddx3y, Kdm5d, Eif2s3y, and Uty). Accordingly, no sex covariate was required for linear modeling.

Re-analysis of GSE324339 was performed to evaluate computational reproducibility against the published article (Palikhe et al., Metabolism 2026) rather than independent biological validation. Of the candidate upregulated genes reported in the original study, four (Pnpla5, Adcy1, Sncaip, and Ltc4s) were excluded during low-count filtering. Among the nine retained transcripts, only Apoa4, Serpina12, and Fgf21 (met our significance criteria padj < 0.05), and all three remained significant following the exclusion of outlier sample Brca1het_HFD_5 (Apoa4 padj 1e-8; Serpina12 5e-3; Fgf21 0.03). In contrast, Scd1, Cdkn1a, and Cidec showed non-significant upward trends (padj 0.12-0.19, unchanged by Brca1het_HFD_5 removal), while Egr1, Ppp1r3g, and Pnpla3 exhibited high standard errors and near-zero shrunken log2FC estimates, rendering their raw fold-change values unreliable.

As a robustness check alongside our primary 12-sample analysis, we repeated differential expression modeling after excluding the divergent sample Brca1het_HFD_5 (n = 5 vs. n = 6). In this 11-sample subset, Brca1 downregulation remained stable (log2FC -1.59, padj 6.7e-4), while the total number of significant transcripts (padj < 0.05) narrowed to 38 (21 upregulated, 17 downregulated) and those exceeding |log2FC| >= 1 decreased to 9. Notably, 7 of the 14 high-magnitude genes from the full cohort were retained (Dntt, Brca1, Serpina12, Apoa4, Arap3, Hhipl2, and Fgf21), whereas the other 7 were driven by Brca1het_HFD_5.

Exploratory GSEA across all 12 samples identified 10 significant MSigDB Hallmark gene sets (padj < 0.05). Brca1+/- livers showed positive enrichment for cholesterol homeostasis (NES 2.31, strongest), EMT, mTORC1 signaling, p53 signaling, hypoxia, adipogenesis, and TNFa signaling, alongside negative enrichment for bile acid metabolism, E2F targets, and MYC targets V1. Among 3,640 GO Biological Process terms, 116 reached significance (padj < 0.05; 64 upregulated, 52 downregulated). Accounting for gene-set redundancy within the GO hierarchy, these terms consolidated into two primary axes: upregulation of lipid (triglyceride, glycerolipid, neutral lipid), glutathione, and carbohydrate-derivative catabolic processes, and downregulation of mRNA processing, mRNA catabolism, and RNA stability.

7 of the 10 Hallmark sets remain significant with the same direction: cholesterol homeostasis (NES 2.14), adipogenesis, p53 pathway, EMT, mTORC1 signalling (up) and E2F targets, bile acid metabolism (down). Hypoxia, TNFa signalling and MYC targets V1 are no longer significant and are treated as dependent on that sample. UV response DN (down) becomes significant and oxidative phosphorylation (up, padj 0.06) emerges, but neither was significant with all 12 samples and neither is emphasised. In GO, 83 sets are significant and the mRNA processing / catabolism / RNA stability theme stays at the top (NES about -1.7 to -2.1), so it does not depend on one sample. Leading-edge genes (first 15 listed; full lists in the GSEA CSVs): cholesterol homeostasis (Ethe1, Cd9, Pcyt2, Fbxo6, Fads2, Atf5, Anxa5, Actg1, Gusb, S100a11, Plaur, Gldc, Lgals3, Fdps, Pmvk) contains only a few lipid-synthesis genes (Fdps, Pmvk, Fads2) and several inflammation/stromal-associated genes (Cd9, Lgals3, Plaur), so the label should not be read as increased cholesterol synthesis; p53 pathway (S100a10, Ctsd, Cdkn1a, Ddit4, Ccnd3, Procr, Pidd1, Inhbb, Rad9a, Epha2, Tm7sf3, Mknk2, Baiap2, Bax, Prmt2) includes canonical p53 targets and stress/growth-arrest genes; mRNA processing (Thrap3, Hnrnpf, Ncl, Prpf6, Snw1, Srrm1, Pabpc1, Npm1, Smndc1, Wbp11, Cstf2t, Ilf3, Srsf5, Acin1, Grsf1) is mostly splicing/RNA-binding genes plus nucleolar genes (Ncl, Npm1). Its leading edge has 166 genes, of which only 13 overlap the MYC V1 leading edge and 6 the E2F leading edge (about 2 each expected by chance), so it is largely a separate signal from the E2F/MYC result. Size and per-gene evidence: cholesterol homeostasis has 29 leading-edge genes (median log2FC 0.38, 1 individually padj < 0.05) and the p53 pathway 62 (median log2FC 0.30, 2 individually padj < 0.05), i.e. many small coordinated shifts (about 1.2-1.3-fold) rather than a few strong genes; leading-edge genes are selected for being up, so these medians are biased upward. These are interpretations of gene lists, not tested.

No sample was excluded on QC grounds (mapping rates 84.7-89.3%, uniform IU library type).

## **Notes / decisions**

**SRR37536766 (Brca1+/- rep 6)** was first processed manually before the pipeline scripts existed (Salmon with -l IU and no bias-correction flags). Its raw reads were then deleted. To keep all 12 samples processed identically, it was re-downloaded and re-run through the pipeline (its cmd_info.json matches the other samples). The earlier result is kept in results/old/SRR37536766_manual for comparison only and is not used in the analysis.

Disk space is limited (WSL on Windows). Intermediate files (.sra, gentrome.fa, trimmed reads) are deleted once no longer needed. After large deletions, compact the WSL virtual disk (wsl --shutdown, then diskpart > compact vdisk).

Genome build is GRCm38 (Ensembl 102) rather than the newer GRCm39; kept consistent throughout. Published study used HISAT2 2.2.1, featureCounts 2.0.6, DESeq2 1.42.0 and GRCm39

All FASTQ files were verified with gzip -t after download (no corruption).

The library is unstranded (IU), so strandedness does not affect quantification.

## **Citation & credits**

Original dataset: Palikhe S, Qiao L, Kutz C, Liang X, Dhorajiya P, Koestler DC, Bhardwaj P, McCoin CS, Peterson O, Tse HM, Thyfault JP, Brown KA. *Brca1* heterozygosity leads to hepatic steatosis in male and female mice despite sexually dimorphic effects on systemic metabolism. bioRxiv preprint, 2026. DOI 10.64898/2026.02.25.708005 (not peer reviewed at the time of writing). GEO series: GSE324339.

Tools: Salmon, fastp, FastQC, MultiQC, tximport, DESeq2.

## **Contact**

Gavin Mateo Hajas, gavin.hajas@gmail.com. Last updated: 2026-11-8.

