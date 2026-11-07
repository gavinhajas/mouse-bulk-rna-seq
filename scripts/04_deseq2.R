# Inputs : metadata/samples.csv, results/salmon_quant/<SRR>/quant.sf, ref/tx2gene.tsv
# Outputs: results/deseq2/  (tables + plots + sessionInfo)

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# settings
samples_file <- "metadata/samples.csv"
quant_dir    <- "results/salmon_quant"
tx2gene_file <- "ref/tx2gene.tsv"
out_dir      <- "results/deseq2"
# optional sensitivity run: EXCLUDE=Brca1het_HFD_5 Rscript scripts/04_deseq2.R
# so that the sample is not included in the DESeq2 analysis (e.g. if it failed QC or is an outlier).
exclude <- trimws(strsplit(Sys.getenv("EXCLUDE", ""), ",")[[1]])
if (length(exclude)) out_dir <- paste0(out_dir, "_excl_", paste(exclude, collapse = "_"))
padj_cut     <- 0.05   # adjusted p-value threshold we use benjamini-hochberg FDR so padj_cut = 0.05 is standard
lfc_cut      <- 1      # |log2 fold change| threshold (= 2-fold)
min_count    <- 10     # pre-filter: >= min_count reads ...

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Sample table 
# Column names are matched loosely so this works with slightly different headers.
samples <- read.csv(samples_file, stringsAsFactors = FALSE)
names(samples) <- tolower(names(samples))
pick <- function(cands) {
  hit <- intersect(cands, names(samples))
  if (!length(hit)) stop("samples.csv needs one of these columns: ", paste(cands, collapse = ", "),
                         "\nFound: ", paste(names(samples), collapse = ", "))
  hit[1]
}
srr_col   <- pick(c("srr", "run", "accession"))
name_col  <- pick(c("name", "sample", "sample_name"))
group_col <- pick(c("group", "condition"))

coldata <- data.frame(
  srr       = samples[[srr_col]],
  condition = factor(samples[[group_col]], levels = c("control", "experimental")),
  row.names = samples[[name_col]],
  stringsAsFactors = FALSE
)
if (any(is.na(coldata$condition)))
  stop("Group column must contain only 'control' and 'experimental'.")

if (length(exclude)) {
  coldata <- coldata[!rownames(coldata) %in% exclude, , drop = FALSE]
  message("Excluding: ", paste(exclude, collapse = ", "), " -> ", nrow(coldata), " samples")
}
# quant.sf location: results/salmon_quant/<SRR>/quant.sf (falls back to sample name)
files <- sapply(seq_len(nrow(coldata)), function(i) {
  cands <- file.path(quant_dir, c(coldata$srr[i], rownames(coldata)[i]), "quant.sf")
  cands[file.exists(cands)][1]
})
names(files) <- rownames(coldata)
if (anyNA(files)) stop("Missing quant.sf for: ", paste(names(files)[is.na(files)], collapse = ", "))
message("Found quant.sf for all ", length(files), " samples.")

# 2. tx2gene
first_line <- readLines(tx2gene_file, n = 1)
has_header <- !grepl("^ENSMUST", first_line)
t2g <- read.delim(tx2gene_file, header = has_header, stringsAsFactors = FALSE)
names(t2g)[1:3] <- c("tx", "gene", "symbol")[seq_len(min(3, ncol(t2g)))]
if (!"symbol" %in% names(t2g)) t2g$symbol <- t2g$gene
# Strip Ensembl version suffixes on both sides so IDs always match
t2g$tx   <- sub("\\.\\d+$", "", t2g$tx)
t2g$gene <- sub("\\.\\d+$", "", t2g$gene)
gene_info <- t2g[!duplicated(t2g$gene), c("gene", "symbol")]

## 3. tximport -> DESeq2
txi <- tximport(files, type = "salmon", tx2gene = t2g[, c("tx", "gene")],
                ignoreTxVersion = TRUE, dropInfReps = TRUE)

dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)

# Pre-filter low-count genes
min_samples <- min(table(coldata$condition))
keep <- rowSums(counts(dds) >= min_count) >= min_samples
message(sum(keep), " of ", nrow(dds), " genes kept after filtering.")
dds <- dds[keep, ]

dds <- DESeq(dds)
message("Coefficients: ", paste(resultsNames(dds), collapse = ", "))

# log2FC > 0 means higher in Brca1+/- than in WT
res <- results(dds, contrast = c("condition", "experimental", "control"), alpha = padj_cut)
summary(res)

# Shrunken LFCs: better for ranking/plotting (noisy low-count genes pulled toward 0).
# Significance (padj) always comes from the unshrunken test above.
coef_name <- "condition_experimental_vs_control"
shr_type <- if (requireNamespace("apeglm", quietly = TRUE)) "apeglm" else
            if (requireNamespace("ashr", quietly = TRUE)) "ashr" else "normal"
message("LFC shrinkage method: ", shr_type)
shr <- if (shr_type == "normal") {
  lfcShrink(dds, contrast = c("condition", "experimental", "control"), type = "normal")
} else {
  lfcShrink(dds, coef = coef_name, type = shr_type)
}

# 4. Result tables
tab <- as.data.frame(res)
tab$gene_id   <- rownames(tab)
tab$symbol    <- gene_info$symbol[match(tab$gene_id, gene_info$gene)]
tab$log2FC_shrunk <- shr$log2FoldChange
tab <- tab[order(tab$padj), c("gene_id", "symbol", "baseMean", "log2FoldChange",
                              "log2FC_shrunk", "lfcSE", "stat", "pvalue", "padj")]
tab$status <- "NS"
tab$status[!is.na(tab$padj) & tab$padj < padj_cut & tab$log2FoldChange >=  lfc_cut] <- "Up in Brca1+/-"
tab$status[!is.na(tab$padj) & tab$padj < padj_cut & tab$log2FoldChange <= -lfc_cut] <- "Down in Brca1+/-"

write.csv(tab, file.path(out_dir, "deseq2_all_genes.csv"), row.names = FALSE)
write.csv(subset(tab, status != "NS"), file.path(out_dir, "deseq2_significant.csv"), row.names = FALSE)
write.csv(as.data.frame(counts(dds, normalized = TRUE)), file.path(out_dir, "normalized_counts.csv"))
message("Significant (padj < ", padj_cut, ", |log2FC| >= ", lfc_cut, "): ",
        sum(tab$status != "NS"), "  (up ", sum(tab$status == "Up in Brca1+/-"),
        ", down ", sum(tab$status == "Down in Brca1+/-"), ")")

# 5. Plots
vsd <- vst(dds, blind = TRUE)   # variance-stabilised values for visualisation only
write.csv(assay(vsd), file.path(out_dir, "vst_counts.csv"))
cols <- c(control = "#1b9e77", experimental = "#d95f02")

# PCA: pca is a data frame with PC1, PC2, condition, and name (sample name)
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
p_pca <- ggplot(pca, aes(PC1, PC2, colour = condition, label = name)) +
  geom_point(size = 3) + ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
  scale_colour_manual(values = cols) +
  xlab(paste0("PC1: ", pv[1], "% variance")) + ylab(paste0("PC2: ", pv[2], "% variance")) +
  theme_bw() + ggtitle("PCA (VST counts)")
ggsave(file.path(out_dir, "pca.png"), p_pca, width = 6.5, height = 5, dpi = 200)

# heatmap of sample-to-sample distances (VST)
sd_mat <- as.matrix(dist(t(assay(vsd))))
ann <- data.frame(condition = coldata$condition, row.names = rownames(coldata))
pheatmap(sd_mat, annotation_col = ann, annotation_colors = list(condition = cols),
         clustering_distance_rows = as.dist(sd_mat), clustering_distance_cols = as.dist(sd_mat),
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-sample distances (VST)",
         filename = file.path(out_dir, "sample_distance_heatmap.png"), width = 7, height = 6)

# volcano plot
vol <- tab[!is.na(tab$padj), ]
vol$neglog10 <- -log10(vol$padj)
status_cols <- c("NS" = "grey70", "Up in Brca1+/-" = "#d95f02", "Down in Brca1+/-" = "#1b9e77")
p_vol <- ggplot(vol, aes(log2FoldChange, neglog10, colour = status)) +
  geom_point(size = 0.8, alpha = 0.7) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
  scale_colour_manual(values = status_cols) +
  labs(x = "log2 fold change (Brca1+/- vs WT)", y = "-log10 adjusted p-value", colour = NULL,
       title = "Volcano plot") + theme_bw()
if (requireNamespace("ggrepel", quietly = TRUE)) {
  top <- head(vol[vol$status != "NS", ], 12)       # already sorted by padj
  p_vol <- p_vol + ggrepel::geom_text_repel(data = top, aes(label = symbol),
                                            size = 3, colour = "black", max.overlaps = 20)
}
ggsave(file.path(out_dir, "volcano.png"), p_vol, width = 7, height = 6, dpi = 200)

# MA plot (shrunken LFC)
png(file.path(out_dir, "ma_plot.png"), width = 1400, height = 1100, res = 200)
plotMA(shr, ylim = c(-4, 4), alpha = padj_cut, main = "MA plot (shrunken log2FC)")
dev.off()

# 6. Check Brca1 counts
# A heterozygous loss should show lower Brca1 in the experimental group.
brca1 <- tab[which(tab$symbol == "Brca1"), ]
if (nrow(brca1)) {
  message("Brca1 result:"); print(brca1)
  d <- plotCounts(dds, gene = brca1$gene_id[1], intgroup = "condition", returnData = TRUE)
  d$sample <- rownames(d)
  p_b <- ggplot(d, aes(condition, count, colour = condition, label = sample)) +
    geom_point(size = 3, position = position_jitter(width = 0.1, height = 0, seed = 1)) +
    scale_colour_manual(values = cols) + scale_y_log10() + theme_bw() +
    ggtitle("Brca1 normalised counts") + theme(legend.position = "none")
  ggsave(file.path(out_dir, "brca1_counts.png"), p_b, width = 4.5, height = 4, dpi = 200)
}

# 7. Reproducibility
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Done. Outputs in ", out_dir)
