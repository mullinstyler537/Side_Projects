library(GEOquery)
library(tidyverse)
library(limma)
library(EnhancedVolcano)
library(pheatmap)
library(patchwork)

gse <- getGEO("GSE45827", GSEMatrix = TRUE)

class(gse)

length(gse)

names(gse)

breast <- gse[[1]]

breast

expr <- exprs(breast)

dim(expr)

head(expr)

metadata <- pData(breast)

dim(metadata)

head(metadata)

colnames(metadata)

metadata$title

metadata$source_name_ch1

genes <- fData(breast)

dim(genes)

head(genes)

colnames(genes)

group <- ifelse(grepl("^Basal", metadata$title), "Basal",
         ifelse(grepl("^Normal", metadata$title), "Normal", NA))

table(group, useNA = "always")

keep <- !is.na(group)

expr2 <- expr[, keep]

metadata2 <- metadata[keep, ]

group2 <- factor(group[keep])

table(group2)

expr_t <- t(expr2)

pca <- prcomp(expr_t, scale = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Group = group2
)

library(ggplot2)

ggplot(pca_df,
       aes(PC1, PC2,
           color = Group)) +

  geom_point(size = 4,
             alpha = 0.9) +

  theme_classic(base_size = 15) +

  labs(
    title = "Principal Component Analysis",
    subtitle = "Basal Breast Cancer vs Normal Tissue",
    x = "PC1",
    y = "PC2"
  ) +

  scale_color_manual(values = c(
    Basal = "#C0392B",
    Normal = "#2980B9"
  )) +

  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

pca_plot <- last_plot()

ggsave(
    "PCA_Basal_vs_Normal.png",
    pca_plot,
    width = 7,
    height = 6,
    dpi = 600
)

design <- model.matrix(~0 + group2)

colnames(design) <- levels(group2)

design

fit <- lmFit(expr2, design)

contrast.matrix <- makeContrasts(
    Basal - Normal,
    levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

results <- topTable(
	fit2,
	number = Inf,
	adjust.method = "BH"
)

head(results)

sig <- subset(results,
	adj.P.Val < 0.05 &
	abs(logFC) > 1)

nrow(sig)

write.csv(
	results,
	"results/all_genes.csv"
)

write.csv(
	sig,
	"results/significant_genes.csv"
)

EnhancedVolcano(
    results,
    lab = rownames(results),
    x = "logFC",
    y = "adj.P.Val",
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.5,
    labSize = 3.5,
    title = "Basal vs Normal",
    subtitle = "Differential Expression",
    caption = "limma analysis"
)

png(
    "figures/VolcanoPlot.png",
    width = 3000,
    height = 2500,
    res = 300
)

EnhancedVolcano(
    results,
    lab = rownames(results),
    x = "logFC",
    y = "adj.P.Val",
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.5,
    labSize = 3.5,
    title = "Basal vs Normal",
    subtitle = "Differential Expression"
)

dev.off()

head(sig, 20)

sig$Gene.symbol[1:20]

colnames(results)

top50 <- results[order(results$adj.P.Val), ]

top50 <- top50[1:50, ]

topGenes <- rownames(top50)

heatmap_data <- expr2[topGenes, ]

dim(heatmap_data)

heatmap_scaled <- t(scale(t(heatmap_data)))

annotation <- data.frame(
    Group = group2
)

rownames(annotation) <- colnames(heatmap_scaled)

ann_colors <- list(
    Group = c(
        Basal = "#C0392B",
        Normal = "#2980B9"
    )
)

library(pheatmap)

pheatmap(
    heatmap_scaled,
    annotation_col = annotation,
    annotation_colors = ann_colors,
    show_rownames = FALSE,
    fontsize_col = 10,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(
        c("navy", "white", "firebrick3")
    )(100),
    main = "Top 50 Differentially Expressed Genes"
)

png(
    "figures/Heatmap_Top50.png",
    width = 3200,
    height = 3000,
    res = 300
)

pheatmap(
    heatmap_scaled,
    annotation_col = annotation,
    annotation_colors = ann_colors,
    show_rownames = FALSE,
    fontsize_col = 10,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(
        c("navy", "white", "firebrick3")
    )(100),
    main = "Top 50 Differentially Expressed Genes"
)

dev.off()

colnames(genes)

results$ProbeID <- rownames(results)

genes$ProbeID <- rownames(genes)

results_annotated <- merge(
    results,
    genes,
    by = "ProbeID"
)

head(results_annotated[, c("ProbeID", "logFC", "adj.P.Val")])

write.csv(
    results_annotated,
    "results/annotated_results.csv",
    row.names = FALSE
)

BiocManager::install(c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
))

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

colnames(results_annotated)

genes_of_interest <- unique(results_annotated$Gene.symbol[
  results_annotated$adj.P.Val < 0.05 &
  !is.na(results_annotated$Gene.symbol)
])

gene_df <- bitr(
  genes_of_interest,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

head(gene_df)

ego <- enrichGO(
  gene          = gene_df$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(ego)

dotplot(ego, showCategory = 15) +
  ggtitle("GO Enrichment: Basal vs Normal")

png(
  "figures/GO_dotplot.png",
  width = 3000,
  height = 2500,
  res = 300
)

dotplot(ego, showCategory = 15) +
  ggtitle("GO Enrichment: Basal vs Normal")

dev.off()

ekegg <- enrichKEGG(
  gene         = gene_df$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

dotplot(ekegg, showCategory = 10) +
  ggtitle("KEGG Pathway Enrichment")

if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("results")) dir.create("results")

gene_annot <- genes[, c("ID", "Gene Symbol", "ENTREZ_GENE_ID")]
colnames(gene_annot) <- c("ProbeID", "Symbol", "ENTREZID")

results$ProbeID <- rownames(results)

results_annotated <- merge(
  results,
  gene_annot,
  by = "ProbeID"
)

results_annotated <- results_annotated[
  !is.na(results_annotated$ENTREZID) &
  results_annotated$ENTREZID != "",
]

sig_entrez <- unique(
  results_annotated$ENTREZID[
    results_annotated$adj.P.Val < 0.05 &
    abs(results_annotated$logFC) > 1
  ]
)

sig_entrez <- as.character(sig_entrez)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

ego <- enrichGO(
  gene = sig_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

summary(ego)

dotplot(ego, showCategory = 15) +
  ggtitle("GO Enrichment: Basal vs Normal")

png(
  "figures/GO_dotplot.png",
  width = 3000,
  height = 2500,
  res = 300
)

dotplot(ego, showCategory = 15) +
  ggtitle("GO Enrichment: Basal vs Normal")

dev.off()

ekegg <- enrichKEGG(
  gene = sig_entrez,
  organism = "hsa",
  pvalueCutoff = 0.05
)

dotplot(ekegg, showCategory = 10) +
  ggtitle("KEGG Pathway Enrichment")

png(
  "figures/KEGG_dotplot.png",
  width = 3000,
  height = 2500,
  res = 300
)

dotplot(ekegg, showCategory = 10) +
  ggtitle("KEGG Pathway Enrichment")

dev.off()

library(ggplot2)
library(patchwork)
library(EnhancedVolcano)
library(pheatmap)

pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_classic(base_size = 14) +
  scale_color_manual(values = c(Basal = "#C0392B", Normal = "#2980B9")) +
  labs(title = "PCA: Basal vs Normal")

volcano_plot <- EnhancedVolcano(
  results,
  lab = rownames(results),
  x = "logFC",
  y = "adj.P.Val",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2,
  labSize = 3,
  title = "Differential Expression",
  subtitle = "Basal vs Normal"
)

heatmap_plot <- pheatmap(
  heatmap_scaled,
  annotation_col = annotation,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  silent = TRUE
)

heatmap_grob <- wrap_elements(full = heatmap_plot$gtable)

go_plot <- dotplot(ego, showCategory = 10) +
  ggtitle("GO Enrichment")

go_grob <- wrap_elements(full = ggplotGrob(go_plot))

final_figure <- (pca_plot | volcano_plot) /
                (heatmap_grob | go_grob)

final_figure

ggsave(
  "figures/Figure1_Publication.png",
  final_figure,
  width = 14,
  height = 10,
  dpi = 600
)

install.packages("gridExtra")
library(gridExtra)
library(ggplot2)
library(patchwork)

pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) +
  theme_classic(base_size = 14) +
  scale_color_manual(values = c(Basal = "#C0392B", Normal = "#2980B9")) +
  labs(title = "A: PCA")

top_labels <- head(rownames(results[order(results$adj.P.Val), ]), 10)

volcano_plot <- EnhancedVolcano(
  results,
  lab = rownames(results),
  x = "logFC",
  y = "adj.P.Val",
  pCutoff = 0.05,
  FCcutoff = 1,
  selectLab = top_labels,
  drawConnectors = TRUE,
  title = "B: Volcano Plot"
)

heatmap_plot <- pheatmap(
  heatmap_scaled,
  annotation_col = annotation,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  silent = TRUE
)

heatmap_grob <- grid::grid.grabExpr(
  grid::grid.draw(heatmap_plot$gtable)
)

go_plot <- dotplot(ego, showCategory = 10) +
  ggtitle("D: GO Enrichment")

final_figure <- (pca_plot | volcano_plot) /
                (heatmap_plot$gtable | go_plot)

ggsave(
  "figures/Figure1_Publication.pdf",
  final_figure,
  width = 14,
  height = 10
)

ggsave(
  "figures/Figure1_Publication.png",
  final_figure,
  width = 14,
  height = 10,
  dpi = 600
)

heatmap_plot <- pheatmap(
  heatmap_scaled,
  annotation_col = annotation,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  silent = TRUE
)

heatmap_panel <- wrap_elements(full = heatmap_plot$gtable)
go_panel <- wrap_elements(full = ggplotGrob(go_plot))

final_figure <- (pca_plot | volcano_plot) /
                (heatmap_panel | go_panel)
