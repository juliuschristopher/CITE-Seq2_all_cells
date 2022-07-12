## ##Analysis of individual clusters####
####Set up####
setwd("/Volumes/GoogleDrive/Shared drives/Okkengroup/Experiments/Julius/Experiments/CITE-Sequencing/CITE-Seq (2)/Overall_analysis/CITE-Seq2_all_cells")

##Packages
library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(Matrix)
library(RColorBrewer)
library(writexl)
library(ggridges)
library(clustree)
library(scRepertoire)
library(future)
library(alakazam)
library(immunarch)
library(airr)
library(biomaRt)
library(SeuratDisk)
library(SeuratData)
library(stringr)
library(viridis)
library(escape)
library(dittoSeq)
library(SingleCellExperiment)
library(ggsci)
library(pals)
library(harmony)
library(gridExtra)
library(scales)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationHub)
library(org.Mm.eg.db)
library(ggpubr)
library(data.table)
library(Polychrome)
library(tidyr)
library(EnhancedVolcano)

##Functions##
tfidf = function(data,target,universe){
  if(!all(target %in% universe))
    stop('Target must be a subset of universe')
  nObs = Matrix::rowSums(data[,target,drop=FALSE]>0)
  nTot = Matrix::rowSums(data[,universe,drop=FALSE]>0)
  tf = nObs/length(target)
  idf = log(length(universe)/nTot)
  score = tf*idf
  #Calculate p-value for significance based on using a hypergeometric distribution to simulate the results of infinite random sampling
  pvals = phyper(nObs-1,nTot,length(universe)-nTot,length(target),lower.tail=FALSE)
  qvals = p.adjust(pvals,method='BH')
  ntf = (exp(-idf)*length(universe)-tf*length(target))/(length(universe)-length(target))
  return(data.frame(geneFrequency=tf,
                    geneFrequencyOutsideCluster=ntf,
                    geneFrequencyGlobal=exp(-idf),
                    geneExpression=Matrix::rowMeans(data[,target,drop=FALSE]),
                    geneExpressionOutsideCluster = Matrix::rowMeans(data[,universe[!(universe%in%target)],drop=FALSE]),
                    geneExpressionGlobal = Matrix::rowMeans(data),
                    idf=idf,
                    tfidf=score,
                    qval=qvals)[order(score,decreasing=TRUE),])
}

##Colours##
col_con <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
col_con1 <- viridis(50)
col_con2 <- createPalette(50,  c("#ff0000", "#00ff00", "#0000ff"))
col_con2 <-as.character(col_con2)

swatch(col_con2)

##Load data##
All_cells <- LoadH5Seurat("CITE-Seq2_all_cells_annot.h5seurat")
Idents(All_cells) <- All_cells$seurat_clusters

All_cells_p1 <- DimPlot(All_cells, label = TRUE, repel = TRUE, reduction = "rna.umap", pt.size = 1.3, label.size = 3, label.box = TRUE, cols = col_con2) +  ggtitle("Seurat Clusters") + theme_bw() + NoLegend()
All_cells_p1 <- All_cells_p1 + theme(plot.title = element_text(color="black", size=20, face = "bold")) + xlab("UMAP1") + ylab("UMAP2")
All_cells_p1
ggsave("All_cells_p1.pdf", width = 30, height = 20, units = "cm")

All_cells_p2 <- DimPlot(All_cells, label = TRUE, repel = TRUE, reduction = "adt.umap", pt.size = 1.3, label.size = 3, label.box = TRUE, cols = col_con2) +  ggtitle("Seurat Clusters") + theme_bw() + NoLegend()
All_cells_p2 <- All_cells_p2 + theme(plot.title = element_text(color="black", size=20, face = "bold")) + xlab("UMAP1") + ylab("UMAP2")
All_cells_p2
ggsave("All_cells_p2.pdf", width = 30, height = 20, units = "cm")

All_cells_p3 <- DimPlot(All_cells, label = TRUE, repel = TRUE, reduction = "wnn.umap", pt.size = 1.3, label.size = 3, label.box = TRUE, cols = col_con2) +  ggtitle("Seurat Clusters") + theme_bw() + NoLegend()
All_cells_p3 <- All_cells_p3 + theme(plot.title = element_text(color="black", size=20, face = "bold")) + xlab("UMAP1") + ylab("UMAP2")
All_cells_p3
ggsave("All_cells_p3.pdf", width = 30, height = 20, units = "cm")


####CD4+and CD8+ Cytotoxic Effector/Exhausted T cells -  (5)####
Idents(All_cells) <- All_cells$wsnn_res.0.7
All_cells_WNN_c5 <- FindMarkers(All_cells, ident.1 = 5, assay = "ADT")
All_cells_WNN_c5 <- All_cells_WNN_c5 %>%
  filter(p_val_adj <= 0.05) %>%
  arrange(desc(avg_log2FC))

All_cells_WNN_c5_rna <- FindMarkers(All_cells, ident.1 = 5, assay = "RNA")
All_cells_WNN_c5_rna <- All_cells_WNN_c5_rna %>%
  filter(p_val_adj <= 0.05) %>%
  arrange(desc(avg_log2FC))

DefaultAssay(All_cells)<-"RNA"
Allcells_TFIDF_c5_genes <- WhichCells(object = All_cells, ident = "5")
Allcells_TFIDF_c5 <- tfidf(GetAssayData(All_cells), Allcells_TFIDF_c5_genes, colnames(All_cells))

#CD4 and CD8 expression in the same cluster but on different cells - other shared signature (cytotoxic? exhuasted?)
DefaultAssay(All_cells)<-"ADT"
FeaturePlot(All_cells, features = c("Cd4", "Cd8a"), blend = TRUE, reduction = "wnn.umap")

#Cytotoxic signature

#Exhausted signature



