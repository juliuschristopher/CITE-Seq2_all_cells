## ##CITE-Seq2 all cells - basic analysis####
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
All_cells$seurat_clusters <- All_cells$wsnn_res.0.7
Idents(All_cells) <- All_cells$seurat_clusters

All_cells_p3 <- DimPlot(All_cells, label = TRUE, reduction = "wnn.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, cols = col_con2) +  ggtitle("Seurat Clusters") + theme_bw() + NoLegend()
All_cells_p3 <- All_cells_p3 + theme(plot.title = element_text(color="black", size=25, face = "bold")) + xlab("UMAP1") + ylab("UMAP2")

####Analysis by metadata variables####
##Cell cycle##
S.genes <- cc.genes.updated.2019$s.genes
S.genes <- lapply(S.genes, str_to_title)
G2M.genes <-  cc.genes.updated.2019$g2m.genes
G2M.genes <- lapply(G2M.genes, str_to_title)
All_cells <- CellCycleScoring(All_cells, s.features=S.genes, g2m.features=G2M.genes, set.ident = FALSE)
Idents(All_cells)
head(All_cells[[]])

All_cells_cell_cylce <- DimPlot(All_cells, label = FALSE, reduction = "wnn.umap", group.by = "Phase", pt.size = 1.3) +  ggtitle("Highlighted by cell-cycle stage") + theme_void()
All_cells_cell_cylce <- All_cells_cell_cylce + theme(plot.title = element_text(color="black", size=15, face="bold"))
All_cells_cell_cylce
ggsave("Cell_cylce_1.pdf", width = 30, height = 20, units = "cm")

All_cells_cell_cylce_ab <- DimPlot(All_cells, label = FALSE, reduction = "adt.umap", group.by = "Phase", pt.size = 1.3) +  ggtitle("Highlighted by cell-cycle stage") + theme_void()
All_cells_cell_cylce_ab <- All_cells_cell_cylce_ab + theme(plot.title = element_text(color="black", size=15, face="bold"))
All_cells_cell_cylce_ab
ggsave("Cell_cylce_2.pdf", width = 30, height = 20, units = "cm")

All_cells_cell_cylce_rna <- DimPlot(All_cells, label = FALSE, reduction = "rna.umap", group.by = "Phase", pt.size = 1.3) +  ggtitle("Highlighted by cell-cycle stage") + theme_void()
All_cells_cell_cylce_rna <- All_cells_cell_cylce_rna + theme(plot.title = element_text(color="black", size=15, face="bold"))
All_cells_cell_cylce_rna
ggsave("Cell_cylce_3.pdf", width = 30, height = 20, units = "cm")

##Mitochondrial genes##
meta.data <- All_cells@meta.data
colnames(meta.data)
meta.data_mt <- meta.data[, c(13:15, 17)]

meta_data_mt_ordered <- meta.data_mt %>%
  group_by(Genotype, Mouse, seurat_clusters) %>%
  summarise_at(vars(percent.Mt), list(mean_mt = mean))

bpr_mt1 <- ggplot(meta_data_mt_ordered, aes(fill=seurat_clusters, y=mean_mt, x=Genotype)) + 
  geom_bar(position="dodge", stat="identity") + theme_bw() + scale_fill_manual(values = col_con2) +
  xlab("Genotype") + ylab("% of Total cells") + guides(fill=guide_legend(title="Clusters")) +
  theme(legend.title = element_text(face = "bold")) + ggtitle("Cellular Proportions") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
bpr_mt1
ggsave("bpr_mt1.pdf", width = 30, height = 20, units = "cm")

bpr_mt2 <- ggboxplot(meta_data_mt_ordered, x = "seurat_clusters", y = "mean_mt",
                 color = "Genotype", palette = c("WT" = "black","BCL6" = "darkgoldenrod1","E1020K" = "blue", "E1020K_BCL6" = "red"), add = "jitter", shape = "Genotype", outlier.shape = NA
) + labs(x = "Seurat Clusters", y = "% of total cells", color = "Genotype", shape = "Genotype") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        axis.title.x=element_blank()) + theme_bw()
bpr_mt2
ggsave("bpr_mt2.pdf", width = 30, height = 20, units = "cm")

##By Sex##
head(All_cells[[]])
Idents(All_cells) <- All_cells$orig.ident
All_cells <- RenameIdents(All_cells, `a` = "Female", `b` ="Male", `c` ="Female", `d` = "Male", `f` = "Male",`g` ="Male", `h` = "Male")
All_cells[["Sex"]] <- Idents(All_cells)
Idents(All_cells) <- All_cells$seurat_clusters

All_cells_sex <- DimPlot(All_cells, label = FALSE, reduction = "wnn.umap", group.by = "Sex", pt.size = 1.3, order = c("Female", "Male")) +  ggtitle("Highlighted by sex") + theme_bw()
All_cells_sex <- All_cells_sex + theme(plot.title = element_text(color="black", size=15, face="bold"))
All_cells_sex
ggsave("All_cells_sex_1.pdf", width = 30, height = 20, units = "cm")

All_cells_sex_ab <- DimPlot(All_cells, label = FALSE, reduction = "adt.umap", group.by = "Sex", pt.size = 1.3, order = c("Female", "Male")) +  ggtitle("Highlighted by sex") + theme_bw()
All_cells_sex_ab <- All_cells_sex_ab + theme(plot.title = element_text(color="black", size=15, face="bold"))
All_cells_sex_ab
ggsave("All_cells_sex_2.pdf", width = 30, height = 20, units = "cm")

All_cells_sex_rna <- DimPlot(All_cells, label = FALSE, reduction = "rna.umap", group.by = "Sex", pt.size = 1.3, order = c("Female", "Male")) +  ggtitle("Highlighted by sex") + theme_bw()
All_cells_sex_rna <- All_cells_sex_rna + theme(plot.title = element_text(color="black", size=15, face="bold"))
All_cells_sex_rna
ggsave("All_cells_sex_3.pdf", width = 30, height = 20, units = "cm")

##By Genotype##
All_cells_Gen <- DimPlot(All_cells, label = FALSE, reduction = "wnn.umap", group.by = "Genotype", pt.size = 1.3, order = c("Female", "Male")) +  ggtitle("Highlighted by sex") + theme_bw()
All_cells_Gen <- All_cells_Gen + theme(plot.title = element_text(color="black", size=15, face="bold"))
All_cells_Gen
ggsave("All_cells_Gen_1.pdf", width = 30, height = 20, units = "cm")

#WT
Sample.WT <- subset(All_cells, subset = Genotype == "WT")
Sample.WT.plot_sex_1 <- DimPlot(Sample.WT, label = FALSE ,reduction = "wnn.umap", group.by = "Sex", pt.size = 1.2, label.size = 6, order = c("Female", "Male")) +
  theme_bw() +
  ggtitle("WT coloured by Sex")
Sample.WT.plot_sex
ggsave("Sample.WT.plot_sex_1.pdf", width = 30, height = 20, units = "cm")

Sample.WT.plot_sex_2 <- DimPlot(Sample.WT, label = FALSE ,reduction = "adt.umap", group.by = "Sex", pt.size = 1.2, label.size = 6, order = c("Female", "Male")) +
  theme_bw() +
  ggtitle("WT coloured by Sex")
Sample.WT.plot_sex_2
ggsave("Sample.WT.plot_sex_2.pdf", width = 30, height = 20, units = "cm")

Sample.WT.plot_sex_3 <- DimPlot(Sample.WT, label = FALSE ,reduction = "rna.umap", group.by = "Sex", pt.size = 1.2, label.size = 6, order = c("Female", "Male")) +
  theme_bw() +
  ggtitle("WT coloured by Sex")
Sample.WT.plot_sex_3
ggsave("Sample.WT.plot_sex_3.pdf", width = 30, height = 20, units = "cm")

Sample.WT.plot_clus_1 <- DimPlot(Sample.WT, label = FALSE ,reduction = "wnn.umap", group.by = "seurat_clusters", pt.size = 1.2, label.size = 6, order = c("Female", "Male")) +
  theme_bw() +
  ggtitle("WT coloured by clusters")
Sample.WT.plot_clus_1
ggsave("Sample.WT.plot_clus_1.pdf", width = 30, height = 20, units = "cm")

Sample.WT.plot_clus_2 <- DimPlot(Sample.WT, label = FALSE ,reduction = "adt.umap", group.by = "seurat_clusters", pt.size = 1.2, label.size = 6, order = c("Female", "Male")) +
  theme_bw() +
  ggtitle("WT coloured by clusters")
Sample.WT.plot_clus_2
ggsave("Sample.WT.plot_clus_2.pdf", width = 30, height = 20, units = "cm")

Sample.WT.plot_clus_3 <- DimPlot(Sample.WT, label = FALSE ,reduction = "rna.umap", group.by = "seurat_clusters", pt.size = 1.2, label.size = 6, order = c("Female", "Male")) +
  theme_bw() +
  ggtitle("WT coloured by clusters")
Sample.WT.plot_clus_3
ggsave("Sample.WT.plot_clus_3.pdf", width = 30, height = 20, units = "cm")


#BCL6
Sample.BCL6 <- subset(All_cells, subset = Genotype == "BCL6")
Sample.BCL6.plot_sex_1 <- DimPlot(Sample.BCL6, label = FALSE ,reduction = "wnn.umap", group.by = "Sex", pt.size = 1.2, label.size = 6, order = c("Female", "Male")) +
  theme_bw() +
  ggtitle("BCL6 coloured by Sex")
Sample.BCL6.plot_sex_1
ggsave("Sample.BCL6.plot_sex_1.pdf", width = 30, height = 20, units = "cm")

Sample.BCL6.plot_sex_2 <- DimPlot(Sample.BCL6, label = FALSE ,reduction = "adt.umap", group.by = "Sex", pt.size = 1.2, label.size = 6, order = c("Female", "Male")) +
  theme_bw() +
  ggtitle("BCL6 coloured by Sex")
Sample.BCL6.plot_sex_2
ggsave("Sample.BCL6.plot_sex_2.pdf", width = 30, height = 20, units = "cm")

Sample.BCL6.plot_sex_3 <- DimPlot(Sample.BCL6, label = FALSE ,reduction = "rna.umap", group.by = "Sex", pt.size = 1.2, label.size = 6, order = c("Female", "Male")) +
  theme_bw() +
  ggtitle("BCL6 coloured by Sex")
Sample.BCL6.plot_sex_3
ggsave("Sample.BCL6.plot_sex_3.pdf", width = 30, height = 20, units = "cm")

Sample.BCL6.plot_clus_1 <- DimPlot(Sample.BCL6, label = FALSE ,reduction = "wnn.umap", group.by = "seurat_clusters", pt.size = 1.2, label.size = 6, order = c("Female", "Male")) +
  theme_bw() +
  ggtitle("BCL6 coloured by clusters")
Sample.BCL6.plot_clus_1
ggsave("Sample.BCL6.plot_clus_1.pdf", width = 30, height = 20, units = "cm")

Sample.BCL6.plot_clus_2 <- DimPlot(Sample.BCL6, label = FALSE ,reduction = "adt.umap", group.by = "seurat_clusters", pt.size = 1.2, label.size = 6, order = c("Female", "Male")) +
  theme_bw() +
  ggtitle("BCL6 coloured by clusters")
Sample.BCL6.plot_clus_2
ggsave("Sample.BCL6.plot_clus_2.pdf", width = 30, height = 20, units = "cm")

Sample.BCL6.plot_clus_3 <- DimPlot(Sample.BCL6, label = FALSE ,reduction = "rna.umap", group.by = "seurat_clusters", pt.size = 1.2, label.size = 6, order = c("Female", "Male")) +
  theme_bw() +
  ggtitle("BCL6 coloured by clusters")
Sample.BCL6.plot_clus_3
ggsave("Sample.BCL6.plot_clus_3.pdf", width = 30, height = 20, units = "cm")


#E1020K
Sample.E1020K <- subset(All_cells, subset = Genotype == "E1020K")
Sample.E1020K.plot_clus_1 <- DimPlot(Sample.E1020K, label = FALSE ,reduction = "wnn.umap", group.by = "seurat_clusters", pt.size = 1.2, label.size = 6, order = c("Female", "Male")) +
  theme_bw() +
  ggtitle("E1020K coloured by clusters")
Sample.E1020K.plot_clus_1
ggsave("Sample.E1020K.plot_clus_1.pdf", width = 30, height = 20, units = "cm")

Sample.E1020K.plot_clus_2 <- DimPlot(Sample.E1020K, label = FALSE ,reduction = "adt.umap", group.by = "seurat_clusters", pt.size = 1.2, label.size = 6, order = c("Female", "Male")) +
  theme_bw() +
  ggtitle("E1020K coloured by clusters")
Sample.E1020K.plot_clus_2
ggsave("Sample.E1020K.plot_clus_2.pdf", width = 30, height = 20, units = "cm")

Sample.E1020K.plot_clus_3 <- DimPlot(Sample.E1020K, label = FALSE ,reduction = "rna.umap", group.by = "seurat_clusters", pt.size = 1.2, label.size = 6, order = c("Female", "Male")) +
  theme_bw() +
  ggtitle("E1020K coloured by clusters")
Sample.E1020K.plot_clus_3
ggsave("Sample.E1020K.plot_clus_3.pdf", width = 30, height = 20, units = "cm")


#E1020K_BCL6
Sample.E1020K_BCL6 <- subset(All_cells, subset = Genotype == "E1020K_BCL6")
Sample.E1020K_BCL6.plot_mouse <- DimPlot(Sample.E1020K_BCL6, label = FALSE ,reduction = "wnn.umap", group.by = "Mouse", pt.size = 1.2, label.size = 6, order = c("E1020K_BCL6_1", "E1020K_BCL6_2")) +
  theme_bw() +
  ggtitle("E1020K_BCL6 coloured by Mouse")
Sample.E1020K_BCL6.plot_mouse
ggsave("Sample.E1020K_BCL6.plot_mouse.pdf", width = 30, height = 20, units = "cm")

Sample.E1020K_BCL6.plot_clus_1 <- DimPlot(Sample.E1020K_BCL6, label = FALSE ,reduction = "wnn.umap", group.by = "seurat_clusters", pt.size = 1.2, label.size = 6, order = c("E1020K_BCL6_1", "E1020K_BCL6_2")) +
  theme_bw() +
  ggtitle("E1020K_BCL6 coloured by cluster")
Sample.E1020K_BCL6.plot_clus_1
ggsave("Sample.E1020K_BCL6.plot_clus_1.pdf", width = 30, height = 20, units = "cm")

Sample.E1020K_BCL6.plot_clus_2 <- DimPlot(Sample.E1020K_BCL6, label = FALSE ,reduction = "adt.umap", group.by = "seurat_clusters", pt.size = 1.2, label.size = 6, order = c("E1020K_BCL6_1", "E1020K_BCL6_2")) +
  theme_bw() +
  ggtitle("E1020K_BCL6 coloured by cluster")
Sample.E1020K_BCL6.plot_clus_2
ggsave("Sample.E1020K_BCL6.plot_clus_2.pdf", width = 30, height = 20, units = "cm")

Sample.E1020K_BCL6.plot_clus_3 <- DimPlot(Sample.E1020K_BCL6, label = FALSE ,reduction = "rna.umap", group.by = "seurat_clusters", pt.size = 1.2, label.size = 6, order = c("E1020K_BCL6_1", "E1020K_BCL6_2")) +
  theme_bw() +
  ggtitle("E1020K_BCL6 coloured by cluster")
Sample.E1020K_BCL6.plot_clus_3
ggsave("Sample.E1020K_BCL6.plot_clus_3.pdf", width = 30, height = 20, units = "cm")


####Cellular proportions between clusters - wnn####
##Percentage per mouse
cell.numbers <- table(All_cells@meta.data$seurat_clusters, All_cells@meta.data$Mouse)
cell.numbers <- as.data.frame.matrix(cell.numbers)
All_cells_meta <- All_cells@meta.data
genotype.numbers <- All_cells_meta %>% dplyr::count(Mouse)
genotype.numbers.vector <- genotype.numbers %>% pull(n)
str(genotype.numbers.vector)
cell.percent <- sweep(cell.numbers, 2, genotype.numbers.vector, "/")
cell.percent <- cell.percent*100

cell_percent_heatmap <- pheatmap::pheatmap(t(cell.percent), cluster_rows = T, cluster_cols = T,show_rownames = T, show_colnames = T,
                                           cellwidth = 20,cellheight = 20, angle_col = 45, display_numbers = T)

cell.percent <- tibble::rownames_to_column(cell.percent, "Cluster")
write_xlsx(cell.percent, "cell.percent.xlsx")
           
##Percentage per Genotype
cell.numbers1 <- table(All_cells@meta.data$seurat_clusters, All_cells@meta.data$Genotype)
cell.numbers1 <- as.data.frame.matrix(cell.numbers1)
All_cells_meta1 <- All_cells@meta.data
genotype.numbers1 <- All_cells_meta1 %>% dplyr::count(Genotype)
genotype.numbers.vector1 <- genotype.numbers1 %>% pull(n)
str(genotype.numbers.vector1)
cell.percent1 <- sweep(cell.numbers1, 2, genotype.numbers.vector1, "/")
cell.percent1 <- cell.percent1*100

cell_percent_heatmap1 <- pheatmap::pheatmap(t(cell.percent1), cluster_rows = T, cluster_cols = T,show_rownames = T, show_colnames = T,
                                           cellwidth = 20,cellheight = 20, angle_col = 45, display_numbers = T)


##Cluster differences as bar chart##
#Boxplot#
meta.data <- All_cells@meta.data
counts <- meta.data %>% group_by(Genotype, Mouse, seurat_clusters) %>% summarise(count = n())
percentage <- counts %>% group_by(Genotype, Mouse) %>% mutate(percent = count/sum(count)*100)

bxp <- ggboxplot(percentage, x = "seurat_clusters", y = "percent",
                 color = "Genotype", palette = c("WT" = "black","BCL6" = "darkgoldenrod1","E1020K" = "blue", "E1020K_BCL6" = "red"), add = "jitter", shape = "Genotype", outlier.shape = NA
) + labs(x = "Seurat Clusters", y = "% of total cells", color = "Genotype", shape = "Genotype") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        axis.title.x=element_blank()) + theme_bw()
bxp
ggsave("bxp.pdf", width = 30, height = 20, units = "cm")


#Bar Chart#
bpr_1 <- ggplot(percentage, aes(fill=seurat_clusters, y=percent, x=Mouse)) + 
  geom_bar(position="fill", stat="identity") + theme_bw() + scale_fill_manual(values = col_con2) +
  xlab("Genotype") + ylab("% of Total cells") + guides(fill=guide_legend(title="Clusters")) +
  theme(legend.title = element_text(face = "bold")) + ggtitle("Cellular Proportions") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
bpr_1
ggsave("BOX_fill_plot.pdf", width = 30, height = 20, units = "cm")


bpr_2 <- ggplot(percentage, aes(fill=seurat_clusters, y=percent, x=Genotype)) + 
  geom_bar(position="fill", stat="identity") + theme_bw() + scale_fill_manual(values = col_con2) +
  xlab("Genotype") + ylab("% of Total cells") + guides(fill=guide_legend(title="Clusters")) +
  theme(legend.title = element_text(face = "bold")) + ggtitle("Cellular Proportions") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
bpr_2
ggsave("BOX_fill_plot2.pdf", width = 30, height = 20, units = "cm")

bpr_3 <- ggplot(percentage, aes(fill=seurat_clusters, y=percent, x=Mouse)) + 
  geom_bar(position="dodge", stat="identity") + theme_bw() + scale_fill_manual(values = col_con2) +
  xlab("Genotype") + ylab("% of Total cells") + guides(fill=guide_legend(title="Clusters")) +
  theme(legend.title = element_text(face = "bold")) + ggtitle("Cellular Proportions") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
bpr_3
ggsave("BOX_fill_plot3.pdf", width = 30, height = 20, units = "cm")

bpr_4 <- ggplot(percentage, aes(fill=seurat_clusters, y=percent, x=Genotype)) + 
  geom_bar(position="dodge", stat="identity") + theme_bw() + scale_fill_manual(values = col_con2) +
  xlab("Genotype") + ylab("% of Total cells") + guides(fill=guide_legend(title="Clusters")) +
  theme(legend.title = element_text(face = "bold")) + ggtitle("Cellular Proportions") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
bpr_4
ggsave("BOX_fill_plot4.pdf", width = 30, height = 20, units = "cm")

#Pie charts
percentage1 <- percentage %>%
  filter(Genotype == "WT") %>%
  group_by(Genotype, seurat_clusters) %>%
  summarise(mean = mean(percent))

percentage2 <- percentage %>%
  filter(Genotype == "BCL6") %>%
  group_by(Genotype, seurat_clusters) %>%
  summarise(mean = mean(percent))

percentage3 <- percentage %>%
  filter(Genotype == "E1020K") %>%
  group_by(Genotype, seurat_clusters) %>%
  summarise(mean = mean(percent))

percentage4 <- percentage %>%
  filter(Genotype == "E1020K_BCL6") %>%
  group_by(Genotype, seurat_clusters) %>%
  summarise(mean = mean(percent))


pie_WT <- ggplot(percentage1, aes(x=Genotype, y=mean, fill=seurat_clusters)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_classic() + scale_fill_manual(values = col_con2) +
  guides(fill=guide_legend(title="Clusters")) +
  xlab("") + ylab("") +
  theme(legend.title = element_text(size = 25, face = "bold")) + ggtitle("Cellular Proportions - WT") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
pie_WT
ggsave("pie_WT.pdf", width = 30, height = 20, units = "cm")

pie_BCL6 <- ggplot(percentage2, aes(x=Genotype, y=mean, fill=seurat_clusters)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_classic() + scale_fill_manual(values = col_con2) +
  guides(fill=guide_legend(title="Clusters")) +
  xlab("") + ylab("") +
  theme(legend.title = element_text(size = 25, face = "bold")) + ggtitle("Cellular Proportions - E1020K") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
pie_BCL6
ggsave("pie_BCL6.pdf", width = 30, height = 20, units = "cm")

pie_E1020K <- ggplot(percentage3, aes(x=Genotype, y=mean, fill=seurat_clusters)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_classic() + scale_fill_manual(values = col_con2) +
  guides(fill=guide_legend(title="Clusters")) +
  xlab("") + ylab("") +
  theme(legend.title = element_text(size = 25, face = "bold")) + ggtitle("Cellular Proportions - WT") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
pie_E1020K
ggsave("pie_E1020K.pdf", width = 30, height = 20, units = "cm")

pie_E1020K_BCL6 <- ggplot(percentage4, aes(x=Genotype, y=mean, fill=seurat_clusters)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_classic() + scale_fill_manual(values = col_con2) +
  guides(fill=guide_legend(title="Clusters")) +
  xlab("") + ylab("") +
  theme(legend.title = element_text(size = 25, face = "bold")) + ggtitle("Cellular Proportions - WT") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
pie_E1020K_BCL6
ggsave("pie_E1020K_BCL6.pdf", width = 30, height = 20, units = "cm")


####Cellular proportions between clusters - ADT####
##Percentage per mouse
cell.numbers_ab <- table(All_cells@meta.data$ADT_snn_res.1.2, All_cells@meta.data$Mouse)
cell.numbers_ab <- as.data.frame.matrix(cell.numbers_ab)
All_cells_meta_ab <- All_cells@meta.data
genotype.numbers_ab <- All_cells_meta_ab %>% dplyr::count(Mouse)
genotype.numbers.vector_ab <- genotype.numbers_ab %>% pull(n)
str(genotype.numbers.vector_ab)
cell.percent_ab <- sweep(cell.numbers_ab, 2, genotype.numbers.vector_ab, "/")
cell.percent_ab <- cell.percent_ab*100

cell_percent_heatmap_ab <- pheatmap::pheatmap(t(cell.percent_ab), cluster_rows = T, cluster_cols = T,show_rownames = T, show_colnames = T,
                                           cellwidth = 20,cellheight = 20, angle_col = 45, display_numbers = T)

##Percentage per Genotype
cell.numbers1_ab <- table(All_cells@meta.data$ADT_snn_res.1.2, All_cells@meta.data$Genotype)
cell.numbers1_ab <- as.data.frame.matrix(cell.numbers1_ab)
All_cells_meta1_ab <- All_cells@meta.data
genotype.numbers1_ab <- All_cells_meta1_ab %>% dplyr::count(Genotype)
genotype.numbers.vector1_ab <- genotype.numbers1_ab %>% pull(n)
str(genotype.numbers.vector1_ab)
cell.percent1_ab <- sweep(cell.numbers1_ab, 2, genotype.numbers.vector1_ab, "/")
cell.percent1_ab <- cell.percent1_ab*100

cell_percent_heatmap1_ab <- pheatmap::pheatmap(t(cell.percent1_ab), cluster_rows = T, cluster_cols = T,show_rownames = T, show_colnames = T,
                                            cellwidth = 20,cellheight = 20, angle_col = 45, display_numbers = T)


##Cluster differences as bar chart##
#Boxplot#
meta.data <- All_cells@meta.data
counts_ab <- meta.data %>% group_by(Genotype, Mouse, ADT_snn_res.1.2) %>% summarise(count = n())
percentage_ab <- counts_ab %>% group_by(Genotype, Mouse) %>% mutate(percent = count/sum(count)*100)

bxp_ab <- ggboxplot(percentage_ab, x = "ADT_snn_res.1.2", y = "percent",
                 color = "Genotype", palette = c("WT" = "black","BCL6" = "darkgoldenrod1","E1020K" = "blue", "E1020K_BCL6" = "red"), add = "jitter", shape = "Genotype", outlier.shape = NA
) + labs(x = "Ab Seurat Clusters", y = "% of total cells", color = "Genotype", shape = "Genotype") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        axis.title.x=element_blank()) + theme_bw()
bxp_ab
ggsave("bxp_ab.pdf", width = 30, height = 20, units = "cm")


#Bar Chart#
bpr_1_ab <- ggplot(percentage_ab, aes(fill=ADT_snn_res.1.2, y=percent, x=Mouse)) + 
  geom_bar(position="fill", stat="identity") + theme_bw() + scale_fill_manual(values = col_con2) +
  xlab("Genotype") + ylab("% of Total cells") + guides(fill=guide_legend(title="Ab Clusters")) +
  theme(legend.title = element_text(face = "bold")) + ggtitle("Cellular Proportions") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
bpr_1_ab
ggsave("BOX_fill_plot_ab.pdf", width = 30, height = 20, units = "cm")


bpr_2_ab <- ggplot(percentage_ab, aes(fill=ADT_snn_res.1.2, y=percent, x=Genotype)) + 
  geom_bar(position="fill", stat="identity") + theme_bw() + scale_fill_manual(values = col_con2) +
  xlab("Genotype") + ylab("% of Total cells") + guides(fill=guide_legend(title="Ab Clusters")) +
  theme(legend.title = element_text(face = "bold")) + ggtitle("Cellular Proportions") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
bpr_2_ab
ggsave("BOX_fill_plot2_ab.pdf", width = 30, height = 20, units = "cm")

bpr_3_ab <- ggplot(percentage_ab, aes(fill=ADT_snn_res.1.2, y=percent, x=Mouse)) + 
  geom_bar(position="dodge", stat="identity") + theme_bw() + scale_fill_manual(values = col_con2) +
  xlab("Genotype") + ylab("% of Total cells") + guides(fill=guide_legend(title="Ab Clusters")) +
  theme(legend.title = element_text(face = "bold")) + ggtitle("Cellular Proportions") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
bpr_3_ab
ggsave("BOX_fill_plot3_ab.pdf", width = 30, height = 20, units = "cm")

bpr_4_ab <- ggplot(percentage_ab, aes(fill=ADT_snn_res.1.2, y=percent, x=Genotype)) + 
  geom_bar(position="dodge", stat="identity") + theme_bw() + scale_fill_manual(values = col_con2) +
  xlab("Genotype") + ylab("% of Total cells") + guides(fill=guide_legend(title="Ab Clusters")) +
  theme(legend.title = element_text(face = "bold")) + ggtitle("Cellular Proportions") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
bpr_4_ab
ggsave("BOX_fill_plot4_ab.pdf", width = 30, height = 20, units = "cm")

#Pie charts
percentage1_ab <- percentage_ab %>%
  filter(Genotype == "WT") %>%
  group_by(Genotype, ADT_snn_res.1.2) %>%
  summarise(mean = mean(percent))

percentage2_ab <- percentage_ab %>%
  filter(Genotype == "BCL6") %>%
  group_by(Genotype, ADT_snn_res.1.2) %>%
  summarise(mean = mean(percent))

percentage3_ab <- percentage_ab %>%
  filter(Genotype == "E1020K") %>%
  group_by(Genotype, ADT_snn_res.1.2) %>%
  summarise(mean = mean(percent))

percentage4_ab <- percentage_ab %>%
  filter(Genotype == "E1020K_BCL6") %>%
  group_by(Genotype, ADT_snn_res.1.2) %>%
  summarise(mean = mean(percent))


pie_WT_ab <- ggplot(percentage1_ab, aes(x=Genotype, y=mean, fill=ADT_snn_res.1.2)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_classic() + scale_fill_manual(values = col_con2) +
  guides(fill=guide_legend(title="Ab Clusters")) +
  xlab("") + ylab("") +
  theme(legend.title = element_text(size = 25, face = "bold")) + ggtitle("Cellular Proportions - WT") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
pie_WT_ab
ggsave("pie_WT_ab.pdf", width = 30, height = 20, units = "cm")

pie_BCL6_ab <- ggplot(percentage2_ab, aes(x=Genotype, y=mean, fill=ADT_snn_res.1.2)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_classic() + scale_fill_manual(values = col_con2) +
  guides(fill=guide_legend(title="Ab Clusters")) +
  xlab("") + ylab("") +
  theme(legend.title = element_text(size = 25, face = "bold")) + ggtitle("Cellular Proportions - E1020K") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
pie_BCL6_ab
ggsave("pie_BCL6_ab.pdf", width = 30, height = 20, units = "cm")

pie_E1020K_ab <- ggplot(percentage3_ab, aes(x=Genotype, y=mean, fill=ADT_snn_res.1.2)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_classic() + scale_fill_manual(values = col_con2) +
  guides(fill=guide_legend(title="Ab Clusters")) +
  xlab("") + ylab("") +
  theme(legend.title = element_text(size = 25, face = "bold")) + ggtitle("Cellular Proportions - WT") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
pie_E1020K_ab
ggsave("pie_E1020K_ab.pdf", width = 30, height = 20, units = "cm")

pie_E1020K_BCL6_ab <- ggplot(percentage4_ab, aes(x=Genotype, y=mean, fill=ADT_snn_res.1.2)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_classic() + scale_fill_manual(values = col_con2) +
  guides(fill=guide_legend(title="Ab Clusters")) +
  xlab("") + ylab("") +
  theme(legend.title = element_text(size = 25, face = "bold")) + ggtitle("Cellular Proportions - WT") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
pie_E1020K_BCL6_ab
ggsave("pie_E1020K_BCL6_ab.pdf", width = 30, height = 20, units = "cm")

####Cellular proportions between clusters - RNA####
##Percentage per mouse
cell.numbers_rna <- table(All_cells@meta.data$RNA_snn_res.0.5, All_cells@meta.data$Mouse)
cell.numbers_rna <- as.data.frame.matrix(cell.numbers_rna)
All_cells_meta_rna <- All_cells@meta.data
genotype.numbers_rna <- All_cells_meta_rna %>% dplyr::count(Mouse)
genotype.numbers.vector_rna <- genotype.numbers_rna %>% pull(n)
str(genotype.numbers.vector_rna)
cell.percent_rna <- sweep(cell.numbers_rna, 2, genotype.numbers.vector_rna, "/")
cell.percent_rna <- cell.percent_rna*100

cell_percent_heatmap_rna <- pheatmap::pheatmap(t(cell.percent_rna), cluster_rows = T, cluster_cols = T,show_rownames = T, show_colnames = T,
                                              cellwidth = 20,cellheight = 20, angle_col = 45, display_numbers = T)

##Percentage per Genotype
cell.numbers1_rna <- table(All_cells@meta.data$RNA_snn_res.0.5, All_cells@meta.data$Genotype)
cell.numbers1_rna <- as.data.frame.matrix(cell.numbers1_rna)
All_cells_meta1_rna <- All_cells@meta.data
genotype.numbers1_rna <- All_cells_meta1_rna %>% dplyr::count(Genotype)
genotype.numbers.vector1_rna <- genotype.numbers1_rna %>% pull(n)
str(genotype.numbers.vector1_rna)
cell.percent1_rna <- sweep(cell.numbers1_rna, 2, genotype.numbers.vector1_rna, "/")
cell.percent1_rna <- cell.percent1_rna*100

cell_percent_heatmap1_rna <- pheatmap::pheatmap(t(cell.percent1_rna), cluster_rows = T, cluster_cols = T,show_rownames = T, show_colnames = T,
                                               cellwidth = 20,cellheight = 20, angle_col = 45, display_numbers = T)


##Cluster differences as bar chart##
#Boxplot#
meta.data <- All_cells@meta.data
counts_rna <- meta.data %>% group_by(Genotype, Mouse, RNA_snn_res.0.5) %>% summarise(count = n())
percentage_rna <- counts_rna %>% group_by(Genotype, Mouse) %>% mutate(percent = count/sum(count)*100)

bxp_rna <- ggboxplot(percentage_rna, x = "RNA_snn_res.0.5", y = "percent",
                    color = "Genotype", palette = c("WT" = "black","BCL6" = "darkgoldenrod1","E1020K" = "blue", "E1020K_BCL6" = "red"), add = "jitter", shape = "Genotype", outlier.shape = NA
) + labs(x = "RNA Seurat Clusters", y = "% of total cells", color = "Genotype", shape = "Genotype") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        axis.title.x=element_blank()) + theme_bw()
bxp_rna
ggsave("bxp_rna.pdf", width = 30, height = 20, units = "cm")


#Bar Chart#
bpr_1_rna <- ggplot(percentage_rna, aes(fill=RNA_snn_res.0.5, y=percent, x=Mouse)) + 
  geom_bar(position="fill", stat="identity") + theme_bw() + scale_fill_manual(values = col_con2) +
  xlab("Genotype") + ylab("% of Total cells") + guides(fill=guide_legend(title="RNA Clusters")) +
  theme(legend.title = element_text(face = "bold")) + ggtitle("Cellular Proportions") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
bpr_1_rna
ggsave("BOX_fill_plot_rna.pdf", width = 30, height = 20, units = "cm")


bpr_2_rna <- ggplot(percentage_rna, aes(fill=RNA_snn_res.0.5, y=percent, x=Genotype)) + 
  geom_bar(position="fill", stat="identity") + theme_bw() + scale_fill_manual(values = col_con2) +
  xlab("Genotype") + ylab("% of Total cells") + guides(fill=guide_legend(title="RNA Clusters")) +
  theme(legend.title = element_text(face = "bold")) + ggtitle("Cellular Proportions") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
bpr_2_rna
ggsave("BOX_fill_plot2_rna.pdf", width = 30, height = 20, units = "cm")

bpr_3_rna <- ggplot(percentage_rna, aes(fill=RNA_snn_res.0.5, y=percent, x=Mouse)) + 
  geom_bar(position="dodge", stat="identity") + theme_bw() + scale_fill_manual(values = col_con2) +
  xlab("Genotype") + ylab("% of Total cells") + guides(fill=guide_legend(title="RNA Clusters")) +
  theme(legend.title = element_text(face = "bold")) + ggtitle("Cellular Proportions") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
bpr_3_rna
ggsave("BOX_fill_plot3_rna.pdf", width = 30, height = 20, units = "cm")

bpr_4_rna <- ggplot(percentage_rna, aes(fill=RNA_snn_res.0.5, y=percent, x=Genotype)) + 
  geom_bar(position="dodge", stat="identity") + theme_bw() + scale_fill_manual(values = col_con2) +
  xlab("Genotype") + ylab("% of Total cells") + guides(fill=guide_legend(title="RNA Clusters")) +
  theme(legend.title = element_text(face = "bold")) + ggtitle("Cellular Proportions") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
bpr_4_rna
ggsave("BOX_fill_plot4_rna.pdf", width = 30, height = 20, units = "cm")

#Pie charts#
percentage1_rna <- percentage_rna %>%
  filter(Genotype == "WT") %>%
  group_by(Genotype, RNA_snn_res.0.5) %>%
  summarise(mean = mean(percent))

percentage2_rna <- percentage_rna %>%
  filter(Genotype == "BCL6") %>%
  group_by(Genotype, RNA_snn_res.0.5) %>%
  summarise(mean = mean(percent))

percentage3_rna <- percentage_rna %>%
  filter(Genotype == "E1020K") %>%
  group_by(Genotype, RNA_snn_res.0.5) %>%
  summarise(mean = mean(percent))

percentage4_rna <- percentage_rna %>%
  filter(Genotype == "E1020K_BCL6") %>%
  group_by(Genotype, RNA_snn_res.0.5) %>%
  summarise(mean = mean(percent))


pie_WT_rna <- ggplot(percentage1_rna, aes(x=Genotype, y=mean, fill=RNA_snn_res.0.5)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_classic() + scale_fill_manual(values = col_con2) +
  guides(fill=guide_legend(title="RNA Clusters")) +
  xlab("") + ylab("") +
  theme(legend.title = element_text(size = 25, face = "bold")) + ggtitle("Cellular Proportions - WT") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
pie_WT_rna
ggsave("pie_WT_rna.pdf", width = 30, height = 20, units = "cm")

pie_BCL6_rna <- ggplot(percentage2_rna, aes(x=Genotype, y=mean, fill=RNA_snn_res.0.5)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_classic() + scale_fill_manual(values = col_con2) +
  guides(fill=guide_legend(title="RNA Clusters")) +
  xlab("") + ylab("") +
  theme(legend.title = element_text(size = 25, face = "bold")) + ggtitle("Cellular Proportions - E1020K") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
pie_BCL6_rna
ggsave("pie_BCL6_rna.pdf", width = 30, height = 20, units = "cm")

pie_E1020K_rna <- ggplot(percentage3_rna, aes(x=Genotype, y=mean, fill=RNA_snn_res.0.5)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_classic() + scale_fill_manual(values = col_con2) +
  guides(fill=guide_legend(title="RNA Clusters")) +
  xlab("") + ylab("") +
  theme(legend.title = element_text(size = 25, face = "bold")) + ggtitle("Cellular Proportions - WT") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
pie_E1020K_rna
ggsave("pie_E1020K_rna.pdf", width = 30, height = 20, units = "cm")

pie_E1020K_BCL6_rna <- ggplot(percentage4_rna, aes(x=Genotype, y=mean, fill=RNA_snn_res.0.5)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_classic() + scale_fill_manual(values = col_con2) +
  guides(fill=guide_legend(title="RNA Clusters")) +
  xlab("") + ylab("") +
  theme(legend.title = element_text(size = 25, face = "bold")) + ggtitle("Cellular Proportions - WT") +
  theme(plot.title = element_text(size = 35, face = "bold")) +
  theme(legend.text = element_text(size = 15))
pie_E1020K_BCL6_rna
ggsave("pie_E1020K_BCL6_rna.pdf", width = 30, height = 20, units = "cm")

