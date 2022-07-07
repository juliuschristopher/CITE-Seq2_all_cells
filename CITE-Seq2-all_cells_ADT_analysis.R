## ##CITE-Seq2 all cells - ADT analysis####
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
All_cells <- LoadH5Seurat("CITE-Seq2_all_cells.h5seurat")
All_cells$seurat_clusters <- All_cells$ADT_snn_res.1.2
Idents(All_cells) <- All_cells$seurat_clusters


All_cells_p1
ggsave("All_cells_p1.pdf", width = 30, height = 20, units = "cm")
All_cells_p2
ggsave("All_cells_p2.pdf", width = 30, height = 20, units = "cm")
All_cells_p3
ggsave("All_cells_p3.pdf", width = 30, height = 20, units = "cm")

####Cluster Identification####
##All ADT markers##
DefaultAssay(All_cells)<-"ADT"
All_cells_p2
CD19 <- FeaturePlot(All_cells, features = "Cd19", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("CD19") + theme(plot.title = element_text(color="black", size=15))
B220 <- FeaturePlot(All_cells, features = "B220", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("B220") + theme(plot.title = element_text(color="black", size=15))
IgM <- FeaturePlot(All_cells, features = "Igm", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("IgM") + theme(plot.title = element_text(color="black", size=15))
IgD <- FeaturePlot(All_cells, features = "Igd", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("IgD") + theme(plot.title = element_text(color="black", size=15))
CD21 <- FeaturePlot(All_cells, features = "Cd21", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("CD21") + theme(plot.title = element_text(color="black", size=15))
CD23 <- FeaturePlot(All_cells, features = "Cd23", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("CD23") + theme(plot.title = element_text(color="black", size=15))
CD93 <- FeaturePlot(All_cells, features = "Cd93", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("CD93") + theme(plot.title = element_text(color="black", size=15))
CD83 <- FeaturePlot(All_cells, features = "Cd83", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("CD83") + theme(plot.title = element_text(color="black", size=15))
CD8a <- FeaturePlot(All_cells, features = "Cd8a", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("CD8a") + theme(plot.title = element_text(color="black", size=15))
CD4 <- FeaturePlot(All_cells, features = "Cd4", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("CD4") + theme(plot.title = element_text(color="black", size=15))
CD40 <- FeaturePlot(All_cells, features = "Cd40", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("CD40") + theme(plot.title = element_text(color="black", size=15))
CD69 <- FeaturePlot(All_cells, features = "Cd69", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("CD69") + theme(plot.title = element_text(color="black", size=15))
CD44 <- FeaturePlot(All_cells, features = "Cd44", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("CD44") + theme(plot.title = element_text(color="black", size=15))
CD62L <- FeaturePlot(All_cells, features = "Cd62l", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("CD62L")+ theme(plot.title = element_text(color="black", size=15))
PD_1 <- FeaturePlot(All_cells, features = "Pd-1", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("PD-1")+ theme(plot.title = element_text(color="black", size=15))
CXCR5 <- FeaturePlot(All_cells, features = "Cxcr5", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("CXCR5")+ theme(plot.title = element_text(color="black", size=15))
CD25 <- FeaturePlot(All_cells, features = "Cd-25", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("CD25")+ theme(plot.title = element_text(color="black", size=15))
GITR <- FeaturePlot(All_cells, features = "Gitr", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("GITR")+ theme(plot.title = element_text(color="black", size=15))
CD69 <- FeaturePlot(All_cells, features = "Cd69", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("CD69")+ theme(plot.title = element_text(color="black", size=15))
PD_L1 <- FeaturePlot(All_cells, features = "Pd-L1", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("PD-L1")+ theme(plot.title = element_text(color="black", size=15))
PD_L2 <- FeaturePlot(All_cells, features = "Pd-L2", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("PD-L2")+ theme(plot.title = element_text(color="black", size=15))
CD86 <- FeaturePlot(All_cells, features = "Cd86", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("CD86")+ theme(plot.title = element_text(color="black", size=15))
CD80 <- FeaturePlot(All_cells, features = "Cd80", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("CD80")+ theme(plot.title = element_text(color="black", size=15))
CD38 <- FeaturePlot(All_cells, features = "Cd38", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("CD38")+ theme(plot.title = element_text(color="black", size=15))
CD95 <- FeaturePlot(All_cells, features = "Cd95", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("CD95")+ theme(plot.title = element_text(color="black", size=15))
ICOS <- FeaturePlot(All_cells, features = "Icos", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("ICOS")+ theme(plot.title = element_text(color="black", size=15))
CTLA4 <- FeaturePlot(All_cells, features = "Ctla4", reduction = "adt.umap", cols = col_con1, pt.size = 1.8) + theme_void() + ggtitle("CTLA4")+ theme(plot.title = element_text(color="black", size=15))


##Find all markers##
All_markers <- FindAllMarkers(All_cells, assay = "ADT", logfc.threshold = 0.25, test.use = "wilcox")

#Cluster 5 - Memory B cells
All_cells_ADT_c5 <- FindMarkers(All_cells, ident.1 = 5, assay = "ADT")
All_cells_ADT_c5 <- All_cells_ADT_c5 %>%
  filter(p_val_adj <= 0.05) %>%
  arrange(desc(avg_log2FC))

All_cells_RNA_c5 <- FindMarkers(All_cells, ident.1 = 5, assay = "RNA")
All_cells_RNA_c5 <- All_cells_RNA_c5 %>%
  filter(p_val_adj <= 0.05) %>%
  arrange(desc(avg_log2FC))

#Cluster 13 - myeloid cells? DCs? Granulocytes?
All_cells_ADT_c13 <- FindMarkers(All_cells, ident.1 = 13, assay = "ADT")
All_cells_ADT_c13 <- All_cells_ADT_c13 %>%
  filter(p_val_adj <= 0.05) %>%
  arrange(desc(avg_log2FC))

All_cells_RNA_c13 <- FindMarkers(All_cells, ident.1 = 13, assay = "RNA")
All_cells_RNA_c13 <- All_cells_RNA_c13 %>%
  filter(p_val_adj <= 0.013) %>%
  arrange(desc(avg_log2FC))

#Cluster 4
All_cells_ADT_c4 <- FindMarkers(All_cells, ident.1 = 4, assay = "ADT")
All_cells_ADT_c4 <- All_cells_ADT_c4 %>%
  filter(p_val_adj <= 0.05) %>%
  arrange(desc(avg_log2FC))

All_cells_RNA_4 <- FindMarkers(All_cells, ident.1 = 4, assay = "RNA")
All_cells_RNA_c4 <- All_cells_RNA_4 %>%
  filter(p_val_adj <= 0.013) %>%
  arrange(desc(avg_log2FC))


DefaultAssay(All_cells)<-"RNA"
FeaturePlot(All_cells, features = "Cyp11a1", reduction = "adt.umap")
