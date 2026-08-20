# ==========================================
# 1. INSTALL PACKAGE
# ==========================================

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(
  c("GEOquery", "limma"),
  ask = FALSE,
  update = FALSE
)

install.packages(c(
  "ggplot2",
  "pheatmap",
  "dplyr"
))


# ==========================================
# 2. LOAD PACKAGE
# ==========================================

library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)
library(dplyr)

# ==========================================
# 3. AMBIL DATASET GEO
# ==========================================

gset <- getGEO("GSE10072", GSEMatrix = TRUE)

length(gset)

# ==========================================
# 4. CEK DATASET
# ==========================================

eset <- gset[[1]]

pData(eset)[, 1:5]

table(pData(eset)$characteristics_ch1.1)

colnames(pData(eset))

unique(pData(eset)$source_name_ch1)

table(pData(eset)$source_name_ch1)

# 5. BUAT GROUP
group <- factor(
  ifelse(
    pData(eset)$source_name_ch1 == "Normal Lung Tissue",
    "Normal",
    "Tumor"
  )
)

table(group)

if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma", ask = FALSE, update = FALSE)
}

library(limma)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

design

# 7. ANALISIS DIFFERENTIAL EXPRESSION

expr <- exprs(eset)

fit <- lmFit(expr, design)

contrast.matrix <- makeContrasts(
  Tumor - Normal,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

DEG <- topTable(
  fit2,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

head(DEG)

# 8. KLASIFIKASI DEG

DEG$Status <- "Not Significant"

DEG$Status[DEG$adj.P.Val < 0.05 & DEG$logFC > 1] <- "Upregulated"

DEG$Status[DEG$adj.P.Val < 0.05 & DEG$logFC < -1] <- "Downregulated"

table(DEG$Status)

# 9. VOLCANO PLOT

library(ggplot2)

DEG$negLog10P <- -log10(DEG$adj.P.Val)

ggplot(DEG, aes(x = logFC, y = negLog10P, color = Status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  ) +
  theme_minimal()

# 10. TOP 50 DEG

top50 <- DEG[order(DEG$adj.P.Val), ]

top50 <- top50[1:50, ]

head(top50)

# 11. HEATMAP TOP 50 DEG

library(pheatmap)

# Ambil ID probe dari top 50
top50_ids <- rownames(top50)

# Ambil ekspresi top 50
heatmap_data <- expr[top50_ids, ]

# Buat heatmap
pheatmap(
  heatmap_data,
  scale = "row",
  show_rownames = FALSE,
  show_colnames = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Top 50 Differentially Expressed Genes"
)

# 12. PERSIAPAN GENE ID UNTUK ENRICHMENT

if (!requireNamespace("hgu133a.db", quietly = TRUE)) {
  BiocManager::install("hgu133a.db", ask = FALSE, update = FALSE)
}

library(hgu133a.db)
library(clusterProfiler)
library(org.Hs.eg.db)

# 13. MAPPING PROBE ID KE ENTREZ GENE ID

library(AnnotationDbi)

# Ambil probe yang signifikan
sig_DEG <- DEG[
  DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > 1,
]

# Ambil ID probe
probe_ids <- rownames(sig_DEG)

# Mapping probe Affymetrix ke Entrez ID
gene_map <- AnnotationDbi::select(
  hgu133a.db,
  keys = probe_ids,
  keytype = "PROBEID",
  columns = c("SYMBOL", "ENTREZID")
)

# Hapus ID yang tidak punya Entrez ID
gene_map <- gene_map[
  !is.na(gene_map$ENTREZID),
]

# Hapus duplikasi Entrez ID
gene_ids <- unique(gene_map$ENTREZID)

length(gene_ids)

# 14. GO ENRICHMENT

ego <- enrichGO(
  gene = gene_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Lihat hasil GO
head(as.data.frame(ego))

# Plot GO enrichment

dotplot(
  ego,
  showCategory = 15,
  title = "GO Enrichment Analysis"
)

# 15. KEGG ENRICHMENT

ekegg <- enrichKEGG(
  gene = gene_ids,
  organism = "hsa",
  keyType = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Plot KEGG enrichment
dotplot(
  ekegg,
  showCategory = 15,
  title = "KEGG Pathway Enrichment"
)

# 16. SIMPAN HASIL ANALISIS

dir.create("hasil", showWarnings = FALSE)

# Simpan Volcano Plot
ggsave(
  "hasil/volcano_plot.png",
  width = 8,
  height = 6,
  dpi = 300
)

volcano_plot <- ggplot(
  DEG,
  aes(x = logFC, y = -log10(adj.P.Val), color = status)
) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  ) +
  theme_minimal()

ggsave(
  "hasil/volcano_plot.png",
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

DEG$status <- "Not Significant"

DEG$status[DEG$adj.P.Val < 0.05 & DEG$logFC > 1] <- "Upregulated"

DEG$status[DEG$adj.P.Val < 0.05 & DEG$logFC < -1] <- "Downregulated"

table(DEG$status)

volcano_plot <- ggplot(
  DEG,
  aes(x = logFC, y = -log10(adj.P.Val), color = status)
) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  ) +
  theme_minimal()

ggsave(
  "hasil/volcano_plot.png",
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

DEG$status <- DEGSstatus

DEG$status <- "Not Significant"

DEG$status[DEG$adj.P.Val < 0.05 & DEG$logFC > 1] <- "Upregulated"

DEG$status[DEG$adj.P.Val < 0.05 & DEG$logFC < -1] <- "Downregulated"

table(DEG$status)

DEG$status <- "Not Significant"

DEG$status[DEG$adj.P.Val < 0.05 & DEG$logFC > 1] <- "Upregulated"

DEG$status[DEG$adj.P.Val < 0.05 & DEG$logFC < -1] <- "Downregulated"

table(DEG$status)

library(ggplot2)

volcano_plot <- ggplot(
  DEG,
  aes(x = logFC, y = -log10(adj.P.Val), color = status)
) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  ) +
  theme_minimal()

volcano_plot

dir.create("hasil", showWarnings = FALSE)

ggsave(
  "hasil/volcano_plot.png",
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# 10. TOP 50 DEG

top50 <- DEG[order(DEG$adj.P.Val), ]

top50 <- top50[1:50, ]

head(top50)

library(pheatmap)

top50_ids <- rownames(top50)

heatmap_data <- expr[top50_ids, ]

pheatmap(
  heatmap_data,
  scale = "row",
  show_rownames = FALSE,
  show_colnames = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Top 50 Differentially Expressed Genes"
)

png(
  "hasil/heatmap_top50.png",
  width = 2400,
  height = 1800,
  res = 300
)

pheatmap(
  heatmap_data,
  scale = "row",
  show_rownames = FALSE,
  show_colnames = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Top 50 Differentially Expressed Genes"
)

dev.off()

pheatmap(
  heatmap_data,
  scale = "row",
  show_rownames = FALSE,
  show_colnames = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Top 50 Differentially Expressed Genes",
  filename = "hasil/heatmap_top50.png",
  width = 8,
  height = 6
)

library(pheatmap)
library(grid)

# Buat heatmap sebagai objek
hm <- pheatmap(
  heatmap_data,
  scale = "row",
  show_rownames = FALSE,
  show_colnames = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Top 50 Differentially Expressed Genes",
  silent = TRUE
)

# Simpan heatmap
png(
  "hasil/heatmap_top50.png",
  width = 2400,
  height = 1800,
  res = 300
)

grid.newpage()
grid.draw(hm$gtable)

dev.off()

dir.create("hasil", showWarnings = FALSE)

png(
  "hasil/GO_enrichment.png",
  width = 2400,
  height = 1800,
  res = 300
)

print(dotplot(
  ego,
  showCategory = 15,
  title = "GO Enrichment Analysis"
))

dev.off()

png(
  "hasil/KEGG_enrichment.png",
  width = 2400,
  height = 1800,
  res = 300
)

print(dotplot(
  ekegg,
  showCategory = 15,
  title = "KEGG Pathway Enrichment"
))

dev.off()

GO_result <- as.data.frame(ego)

write.csv(
  GO_result,
  "hasil/GO_enrichment_results.csv",
  row.names = FALSE
)

KEGG_result <- as.data.frame(ekegg)

write.csv(
  KEGG_result,
  "hasil/KEGG_enrichment_results.csv",
  row.names = FALSE
)

GO_result <- as.data.frame(ego)

nrow(GO_result)

sum(GO_result$p.adjust < 0.05)

nrow(KEGG_result)

sum(KEGG_result$p.adjust < 0.05)

KEGG_result[, c("ID", "Description", "GeneRatio", "p.adjust")]

KEGG_result[, c("ID", "Description", "GeneRatio", "p.adjust", "geneID")]

KEGG_result$geneID

library(org.Hs.eg.db)

KEGG_readable <- setReadable(
  KEGG_result,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

KEGG_readable[, c("ID", "Description", "GeneRatio", "p.adjust", "geneID")]

library(clusterProfiler)
library(org.Hs.eg.db)

gene_list <- unique(unlist(strsplit(KEGG_result$geneID, "/")))

gene_map <- bitr(
  gene_list,
  fromType = "ENTREZID",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)

gene_map

KEGG_result$geneSymbol <- sapply(
  KEGG_result$geneID,
  function(x) {
    ids <- unlist(strsplit(x, "/"))
    symbols <- gene_map$SYMBOL[match(ids, gene_map$ENTREZID)]
    paste(symbols[!is.na(symbols)], collapse = "/")
  }
)

KEGG_result[, c("ID", "Description", "GeneRatio", "p.adjust", "geneSymbol")]

write.csv(
  KEGG_result,
  "hasil/KEGG_enrichment_with_geneSymbol.csv",
  row.names = FALSE
)

KEGG_table <- KEGG_result[, c(
  "ID",
  "Description",
  "GeneRatio",
  "p.adjust",
  "geneSymbol"
)]

KEGG_table

KEGG_table <- KEGG_result[, c(
  "ID",
  "Description",
  "GeneRatio",
  "p.adjust",
  "geneSymbol"
)]

KEGG_table

KEGG_top10 <- head(KEGG_table, 10)

KEGG_top10

write.csv(
  KEGG_top10,
  "hasil/KEGG_top10.csv",
  row.names = FALSE
)

library(ggplot2)

KEGG_top10$Description <- factor(
  KEGG_top10$Description,
  levels = rev(KEGG_top10$Description)
)

KEGG_plot <- ggplot(
  KEGG_top10,
  aes(
    x = -log10(p.adjust),
    y = Description
  )
) +
  geom_point(
    aes(size = as.numeric(sub("/.*", "", GeneRatio)))
  ) +
  labs(
    title = "Top 10 KEGG Enriched Pathways",
    x = "-log10 Adjusted P-value",
    y = "Pathway",
    size = "Gene Count"
  ) +
  theme_minimal()

KEGG_plot

print(KEGG_plot)

ggsave(
  "KEGG_top10.png",
  plot = KEGG_plot,
  width = 10,
  height = 7,
  dpi = 300
)

file.exists("KEGG_top10.png")

KEGG_top10[, c("ID", "Description", "GeneRatio", "p.adjust", "geneSymbol")]

View(KEGG_top10)

write.csv(
  KEGG_top10,
  "KEGG_top10_final_with_geneSymbol.csv",
  row.names = FALSE
)

file.exists("KEGG_top10_final_with_geneSymbol.csv")


ls()

write.csv(
  DEG,
  "DEG_results.csv",
  row.names = FALSE
)

up_genes <- DEG[DEG$Status == "Upregulated", ]
down_genes <- DEG[DEG$Status == "Downregulated", ]

write.csv(
  up_genes,
  "upregulated_genes.csv",
  row.names = FALSE
)

write.csv(
  down_genes,
  "downregulated_genes.csv",
  row.names = FALSE
)