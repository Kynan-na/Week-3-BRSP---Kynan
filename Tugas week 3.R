library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(illuminaHumanv4.db)
library(AnnotationDbi)
library(umap)
library(hgu133a.db)

gset <- getGEO("GSE10072", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

ex <- exprs(gset)

qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))

LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

if (LogTransform) {
  # Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

group_info <- pData(gset)[["source_name_ch1"]]
#make.names(): mengubah teks menjadi format valid untuk R
groups <- make.names(group_info)
#factor():
#Mengubah data kategorik menjadi faktor
#Faktor sangat penting untuk analisis statistik di R
gset$group <- factor(groups)
#levels(): melihat kategori unik dalam faktor
nama_grup <- levels(gset$group)
print(nama_grup)

design <- model.matrix(~0 + gset$group)
#colnames(): memberi nama kolom agar mudah dibaca
colnames(design) <- levels(gset$group)

grup_kanker <- nama_grup[1]
grup_normal <- nama_grup[2]
contrast_formula <- paste(grup_kanker, "-", grup_normal)
print(paste("Kontras yang dianalisis:", contrast_formula))

fit <- lmFit(ex, design)
#makeContrasts(): mendefinisikan perbandingan antar grup
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels
                                 = design)
#contrasts.fit(): menerapkan kontras ke model
fit2 <- contrasts.fit(fit, contrast_matrix)
#eBayes():
#Empirical Bayes untuk menstabilkan estimasi varians
fit2 <- eBayes(fit2)
#topTable():
#Mengambil hasil akhir DEG
#adjust = "fdr" -> koreksi multiple testing
#p.value = 0.01 -> gen sangat signifikan
topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.01
)

head(topTableResults,50)


probe_ids <- rownames(topTableResults)
#Mapping probe -> gene symbol & gene name
gene_annotation <- AnnotationDbi::select(
  hgu133a.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)
#Gabungkan dengan hasil limma
topTableResults$PROBEID <- rownames(topTableResults)
topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)
#Cek hasil anotasi


volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)
#Klasifikasi status gen
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val <
                      0.01] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val <
                      0.01] <- "DOWN"

#Visualisasi
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color =
                           status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot Differentially Expressed Genes Kanker Paru")


# heatmap 
topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]

top50 <- head(topTableResults,50)

mat_heatmap <- ex[top50$PROBEID, ]

gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID, 
  top50$SYMBOL
)

rownames(mat_heatmap) <-gene_label

mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) ==0,
]

gene_variance <- apply(mat_heatmap,1, var)
mat_heatmap <-mat_heatmap[gene_variance >0,]

annotation_col <- data.frame(
  Group = gset$group
)

rownames(annotation_col) <- colnames (mat_heatmap)

pheatmap(
  mat_heatmap,
  annotation_col = annotation_col,
  show_colnames = FALSE,
  show_rownames =TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Top 50 DEGs"
)

library(clusterProfiler)
library(org.Hs.eg.db) # Database anotasi untuk manusia (Homo sapiens)
library(enrichplot)
library(ggplot2)

subset(data.frame())


library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

sig_genes <- subset(topTableResults, !is.na(SYMBOL) & abs(logFC) > 1)

gene_list <- unique(sig_genes$SYMBOL)

gene_entrez <- bitr(gene_list,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)
#GO enrichment
go_enrich <- enrichGO(gene          = gene_entrez$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP", 
                      pAdjustMethod = "BH", 
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

go_dotplot <- dotplot(go_enrich, showCategory=10, title="GO Enrichment (GSE10072)")
print(go_dotplot)

#KEGG Enrichment
kegg_enrich <- enrichKEGG(gene         = gene_entrez$ENTREZID,
                          organism     = 'hsa',
                          pvalueCutoff = 0.05)

kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

kegg_barplot <- barplot(kegg_enrich, showCategory=10, title="KEGG Enrichment (GSE10072)")
print(kegg_barplot)