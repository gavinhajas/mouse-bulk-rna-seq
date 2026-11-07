# Input : results/deseq2/deseq2_all_genes.csv  (from 04_deseq2.R, all 12 samples)
# Output: results/enrichment/
#
# GSEA ranks ALL tested genes by the DESeq2 Wald statistic (sign = direction in Brca1+/-,
# magnitude = strength of evidence) and asks whether gene sets are shifted toward the top or
# bottom of that list. It needs no significance cutoff, which suits a small DE gene list.
# Positive NES = pathway higher in Brca1+/-; negative NES = higher in WT.

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

#settings
res_file <- Sys.getenv("RES_FILE", "results/deseq2/deseq2_all_genes.csv")
out_dir  <- Sys.getenv("OUT_DIR", "results/enrichment")
fdr_cut  <- 0.05       # GSEA is often reported at 0.25 too; both counts are printed
min_size <- 15         # smallest / largest gene set tested (after matching to our genes)
max_size <- 500
# Genes named as upregulated in female Brca1+/- livers in the preprint's RNA-seq figure.
paper_genes <- c("Serpina12", "Fgf21", "Apoa4", "Scd1", "Cdkn1a", "Egr1", "Cidec",
                 "Pnpla3", "Pnpla5", "Ppp1r3g", "Adcy1", "Sncaip", "Ltc4s")


dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

tab <- read.csv(res_file, stringsAsFactors = FALSE)

# 1. Ranked gene list
a <- tab[!is.na(tab$stat) & !is.na(tab$symbol) & tab$symbol != "", ]
a <- a[order(-abs(a$stat)), ]
a <- a[!duplicated(a$symbol), ]            # one Ensembl ID per symbol (strongest signal)
ranks <- sort(setNames(a$stat, a$symbol), decreasing = TRUE)
message("Ranked genes: ", length(ranks))

# 2. Gene sets (mouse, MSigDB via msigdbr)
# this is so that the script can run even if msigdbr is not installed or the MSigDB data package is missing
try_msig <- function(attempts) {
  for (f in attempts) {
    d <- tryCatch(suppressWarnings(f()), error = function(e) NULL)
    if (!is.null(d) && nrow(d) > 0) return(d)
  }
  NULL
}
hallmark <- try_msig(list(
  function() msigdbr(species = "Mus musculus", db_species = "MM", collection = "MH"),
  function() msigdbr(species = "Mus musculus", collection = "H"),
  function() msigdbr(species = "Mus musculus", category = "H")))
gobp <- try_msig(list(
  function() msigdbr(species = "Mus musculus", db_species = "MM", collection = "M5", subcollection = "GO:BP"),
  function() msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP"),
  function() msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP"),
  function() msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")))
if (is.null(hallmark) && is.null(gobp))
  stop("msigdbr returned no gene sets. Check packageVersion('msigdbr') and its documentation ",
       "(newer versions may need the 'msigdbdf' data package).")

# 3. GSEA
run_gsea <- function(sets, label) {
  pw <- lapply(split(sets$gene_symbol, sets$gs_name), unique)
  fg <- as.data.frame(fgseaMultilevel(pathways = pw, stats = ranks,
                                      minSize = min_size, maxSize = max_size))
  fg <- fg[order(fg$padj), ]
  out <- fg
  out$leadingEdge <- vapply(out$leadingEdge, paste, character(1), collapse = ";")
  write.csv(out, file.path(out_dir, paste0("gsea_", label, ".csv")), row.names = FALSE)
  message(label, ": ", nrow(fg), " sets tested; padj < ", fdr_cut, ": ",
          sum(fg$padj < fdr_cut, na.rm = TRUE), " (up ", sum(fg$padj < fdr_cut & fg$NES > 0, na.rm = TRUE),
          ", down ", sum(fg$padj < fdr_cut & fg$NES < 0, na.rm = TRUE),
          "); padj < 0.25: ", sum(fg$padj < 0.25, na.rm = TRUE))
  print(head(fg[, c("pathway", "size", "NES", "pval", "padj")], 10), row.names = FALSE)
  list(fg = fg, pw = pw)
}

plot_nes <- function(fg, label, strip, n = 20) {
  d <- head(fg[order(fg$padj), ], n)
  d$name <- substr(gsub("_", " ", sub(strip, "", d$pathway)), 1, 60)
  d$name <- factor(d$name, levels = d$name[order(d$NES)])
  d$sig  <- factor(d$padj < fdr_cut, levels = c(TRUE, FALSE),
                   labels = c(paste0("padj < ", fdr_cut), "not significant"))
  p <- ggplot(d, aes(NES, name, fill = sig)) +
    geom_col() +
    scale_fill_manual(values = setNames(c("#d95f02", "grey70"), levels(d$sig)), drop = FALSE) +
    labs(x = "Normalised enrichment score (positive = higher in Brca1+/-)", y = NULL, fill = NULL,
         title = paste0("GSEA: ", label, " (top ", nrow(d), " by padj)")) +
    theme_bw()
  ggsave(file.path(out_dir, paste0("gsea_", label, "_top.png")), p,
         width = 9, height = 0.3 * nrow(d) + 2, dpi = 200)
}

if (!is.null(hallmark)) {
  h <- run_gsea(hallmark, "hallmark")
  plot_nes(h$fg, "hallmark", "^HALLMARK_")
  for (nm in head(h$fg$pathway, 4)) {          # classic running-sum plots, top 4 by padj
    p <- plotEnrichment(h$pw[[nm]], ranks) + ggtitle(gsub("_", " ", nm))
    ggsave(file.path(out_dir, paste0("enrichment_", nm, ".png")), p, width = 6, height = 4, dpi = 200)
  }
} else message("Hallmark sets unavailable; skipped.")

if (!is.null(gobp)) {
  g <- run_gsea(gobp, "GO_BP")
  plot_nes(g$fg, "GO_BP", "^GOBP_")
} else message("GO:BP sets unavailable; skipped.")

# 4. Genes named in the preprint
pg <- tab[tab$symbol %in% paper_genes,
          c("symbol", "baseMean", "log2FoldChange", "log2FC_shrunk", "lfcSE", "pvalue", "padj")]
pg <- pg[order(-pg$log2FoldChange), ]
write.csv(pg, file.path(out_dir, "preprint_genes.csv"), row.names = FALSE)
missing <- setdiff(paper_genes, tab$symbol)
if (length(missing)) message("Not in our table (removed by the low-count filter or no symbol): ",
                             paste(missing, collapse = ", "))
message("Preprint genes in our data (", nrow(pg), "):")
print(pg, row.names = FALSE)

if (nrow(pg) > 0) {
  pg$sig <- factor(!is.na(pg$padj) & pg$padj < 0.05, levels = c(TRUE, FALSE),
                   labels = c("padj < 0.05", "padj >= 0.05"))
  pg$symbol <- factor(pg$symbol, levels = pg$symbol[order(pg$log2FoldChange)])
  p_pg <- ggplot(pg, aes(log2FoldChange, symbol, colour = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_pointrange(aes(xmin = log2FoldChange - 1.96 * lfcSE, xmax = log2FoldChange + 1.96 * lfcSE)) +
    scale_colour_manual(values = setNames(c("#d95f02", "grey50"), levels(pg$sig)), drop = FALSE) +
    labs(x = "log2 fold change, Brca1+/- vs WT (unshrunken, 95% CI)", y = NULL, colour = NULL,
         title = "Genes highlighted in the source preprint") + theme_bw()
  ggsave(file.path(out_dir, "preprint_genes.png"), p_pg, width = 6.5, height = 5, dpi = 200)
}

# 5. Reproducibility
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Done. Outputs in ", out_dir)
