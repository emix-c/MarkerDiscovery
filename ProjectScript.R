#use data set from GSE197177
#ZC: normal, #YF: primary PDAC tumor, ZY: w/ liver hepatic metastases 
library(Seurat)
library(tidyverse)
library(gridExtra)
library(glmGamPoi)
library(DoubletFinder)
library(annotables)
library(metap)
library(presto)
library(cowplot)
library(patchwork)
library(ggrepel)
source("/Users/Emix/Desktop/RNASeq/functions.R") #load functions script 


#Step 1: Load in the data of my samples. Each sample will be its own Seurat object. 
setwd('/Users/Emix/Desktop/RNASeq/GSE197177_RAW')
files <- list.files(path = '.', recursive = F, full.names = F)
sample_names <- unique(gsub("^([^_]*_[^_]*)_.*$", "\\1", files))


for (x in sample_names){ 
  counts <- ReadMtx(mtx = paste0(x, '_matrix.mtx.gz'), 
                    features = paste0(x, '_features.tsv.gz'), 
                    cells =  paste0(x, '_barcodes.tsv.gz'))
  #seurat object per sample
  newname <- gsub(".*_","",x)
  assign(gsub("-", "", newname), CreateSeuratObject(counts=counts)) 
}


#Step 2: Each sample must be filtered and processed with doublets removed. 

#1YF sample 
filtered_1YF <- seurat_filtered(Case1YF)
processed_1YF <- seurat_prepro(filtered_1YF)
ElbowPlot(processed_1YF) #find PCAs
processed_1YF <- RunUMAP(processed_1YF, dims = 1:20)
DimPlot(processed_1YF)
pK<- seurat_pk(processed_1YF)
final_1YF <- seurat_df(processed_1YF, pK)
View(final_1YF@meta.data)
table(final_1YF@meta.data$DF.classifications_0.25_0.16_620) 
final_1YF <- subset(final_1YF,DF.classifications_0.25_0.16_620 == 'Singlet' )

#1ZY sample 
filtered_1ZY <- seurat_filtered(Case1ZY)
processed_1ZY <- seurat_prepro(filtered_1ZY)
ElbowPlot(processed_1ZY) #find PCAs
processed_1ZY <- RunUMAP(processed_1ZY, dims = 1:20)
DimPlot(processed_1ZY)
pK<- seurat_pk(processed_1ZY)
final_1ZY <- seurat_df(processed_1ZY, pK)
View(final_1ZY@meta.data)
table(final_1ZY@meta.data$DF.classifications_0.25_0.22_646) 
final_1ZY <- subset(final_1ZY,DF.classifications_0.25_0.22_646 == 'Singlet' )

#2YF sample 
filtered_2YF <- seurat_filtered(Case2YF)
processed_2YF <- seurat_prepro(filtered_2YF)
ElbowPlot(processed_2YF) #find PCAs
processed_2YF <- RunUMAP(processed_2YF, dims = 1:20)
DimPlot(processed_2YF)
pK<- seurat_pk(processed_2YF)
final_2YF<- seurat_df(processed_2YF, pK)
View(final_2YF@meta.data)
table(final_2YF@meta.data$DF.classifications_0.25_0.3_863) 
final_2YF <- subset(final_2YF,DF.classifications_0.25_0.3_863 == 'Singlet' )

#2ZC sample 
filtered_2ZC <- seurat_filtered(Case2ZC)
processed_2ZC <- seurat_prepro(filtered_2ZC)
ElbowPlot(processed_2ZC) #find PCAs
processed_2ZC <- RunUMAP(processed_2ZC, dims = 1:20)
DimPlot(processed_2ZC)
pK<- seurat_pk(processed_2ZC)
final_2ZC<- seurat_df(processed_2ZC, pK)
View(final_2ZC@meta.data)
table(final_2ZC@meta.data$DF.classifications_0.25_0.005_486) 
final_2ZC <- subset(final_2ZC,DF.classifications_0.25_0.005_486 == 'Singlet' )


#2ZY sample 
filtered_2ZY <- seurat_filtered(Case2ZY)
processed_2ZY <- seurat_prepro(filtered_2ZY)
ElbowPlot(processed_2ZY) #find PCAs
processed_2ZY <- RunUMAP(processed_2ZY, dims = 1:20)
DimPlot(processed_2ZY)
pK<- seurat_pk(processed_2ZY)
final_2ZY<- seurat_df(processed_2ZY, pK)
View(final_2ZY@meta.data)
table(final_2ZY@meta.data$DF.classifications_0.25_0.09_692) 
final_2ZY <- subset(final_2ZY,DF.classifications_0.25_0.09_692 == 'Singlet' )

#3YF sample 
filtered_3YF <- seurat_filtered(Case3YF)
processed_3YF <- seurat_prepro(filtered_3YF)
ElbowPlot(processed_3YF) #find PCAs
processed_3YF <- RunUMAP(processed_3YF, dims = 1:20)
DimPlot(processed_3YF)
pK<- seurat_pk(processed_3YF)
final_3YF<- seurat_df(processed_3YF, pK)
View(final_3YF@meta.data)
table(final_3YF@meta.data$DF.classifications_0.25_0.3_622) 
final_3YF <- subset(final_3YF,DF.classifications_0.25_0.3_622 == 'Singlet' )

#3ZY sample 
filtered_3ZY <- seurat_filtered(Case3ZY)
processed_3ZY <- seurat_prepro(filtered_3ZY)
ElbowPlot(processed_3ZY) #find PCAs
processed_3ZY <- RunUMAP(processed_3ZY, dims = 1:20)
DimPlot(processed_3ZY)
pK<- seurat_pk(processed_3ZY)
final_3ZY<- seurat_df(processed_3ZY, pK)
View(final_3ZY@meta.data)
table(final_3ZY@meta.data$DF.classifications_0.25_0.3_594) 
final_3ZY <- subset(final_3ZY,DF.classifications_0.25_0.3_594 == 'Singlet' )

#4ZY sample 
filtered_4ZY <- seurat_filtered(Case4ZY)
processed_4ZY <- seurat_prepro(filtered_4ZY)
ElbowPlot(processed_4ZY) #find PCAs
processed_4ZY <- RunUMAP(processed_4ZY, dims = 1:20)
DimPlot(processed_4ZY)
pK<- seurat_pk(processed_4ZY)
final_4ZY<- seurat_df(processed_4ZY, pK)
View(final_4ZY@meta.data)
table(final_4ZY@meta.data$DF.classifications_0.25_0.27_138) 
final_4ZY <- subset(final_4ZY,DF.classifications_0.25_0.27_138 == 'Singlet' )


#Step 3: Now merge the samples together into one Seurat object. 
merged_seurat <- merge(final_1YF, y = c(final_1ZY,final_2YF,final_2ZC,final_2ZY, 
                                              final_3YF, final_3ZY, final_4ZY ), 
                       #keep track of which dir it came from 
                       add.cell.ids = ls()[1:8], 
                       project = 'Pan')

View(merged_seurat@meta.data)


#Step 4: Prepare object for integration. 

class(merged_seurat[['RNA']])
Layers(merged_seurat[['RNA']])
obj_analysis <- seurat_prepro(merged_seurat)
ElbowPlot(obj_analysis)
obj_analysis <-FindNeighbors(obj_analysis)
obj_analysis<- FindClusters(obj_analysis, resolution = 0.8, cluster.name = 'unintegrated_clusters') 
obj_analysis <- RunUMAP(obj_analysis, dims=1:20 )
View(obj_analysis@meta.data)
obj_analysis@meta.data$FullName <- rownames(obj_analysis@meta.data)

obj_analysis@meta.data <- obj_analysis@meta.data %>%  
  separate(FullName, c('Sample', 'Barcode'), sep = '_')%>%  
  separate(Sample, c('Patient', 'Type'), sep = 5)

unintegrated <- DimPlot(obj_analysis, group.by = 'unintegrated_clusters', label = TRUE) + NoLegend()
ggsave("unintegrated_plot.png", unintegrated, width = 10, height = 14, units = 'in')


#check layers in RNA assay 
Layers(obj_analysis[['RNA']])

options(future.globals.maxSize = 8000 * 1024^2)

#Step 5: Perform integration via CCA. 
integrated_obj <- IntegrateLayers(obj_analysis, method = CCAIntegration, orig.reduction = 'pca', 
                                  new.reduction = 'integrated.cca', verbose = FALSE)
integrated_obj <- FindNeighbors(integrated_obj, reduction = 'integrated.cca', dims = 1:30)
integrated_obj <- FindClusters(integrated_obj, resolution = 0.8, cluster.name = 'integratedclusters')
integrated_obj <- RunUMAP(integrated_obj, reduction = 'integrated.cca', dims = 1:30, reduction.name = 'umap.cca')

View(integrated_obj@meta.data)

integrated <- DimPlot(integrated_obj, reduction = 'umap.cca', group.by = 'integratedclusters', label = TRUE)
ggsave("integrated_plot.png", integrated, width = 10, height = 14, units = 'in')

integratedbytype <- DimPlot(integrated_obj, reduction = 'umap.cca', split.by = 'Type', label = TRUE)
ggsave("integratedbytype_plot.png", integratedbytype, width = 14, height = 14, units = 'in')

#Step 6: Look at cells per cluster for each sample type (YF, ZY, ZC). 
#determine the number of cells per cluster per sample type 
cells <- FetchData(integrated_obj, 
                   vars = c('integratedclusters', 'Type')) %>%  
                   count(integratedclusters, Type) %>%
                   spread(Type, n)

write.csv(cells, 'cell_clusterdata.csv')
  
#notably, there is only ZY cells in cluster 20 
#clusters 20 and 21 have very low cell counts 
#ZC has significantly less cells in clusters 0,2 and 4. 



#Step 7: Study the segregation of clusters by other metrics.            
metrics <-  c("nCount_RNA", "nFeature_RNA", "MTpercent")

other_metrics <- FeaturePlot(integrated_obj, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

ggsave("othermetrics.png", other_metrics, width = 14, height = 14, units = 'in')

#Step 8: Visualize how the clusters separate by different PCs. 

PCcolumns <- c(paste0('PC_', 1:20 ), "integratedclusters", "umap_1", "umap_2")
pc_data <- FetchData(integrated_obj, vars=PCcolumns)
umap_label <- FetchData(integrated_obj,
                        vars = c("integratedclusters", "umap_1", "umap_2"))  %>%
  group_by(integratedclusters) %>%
  summarise(x=mean(umap_1), y=mean(umap_2))


#look at top 20 PCs
PC_clusterplot <- plot_grid(pc("PC_1", pc_data, umap_label), pc("PC_2", pc_data, umap_label),
                            pc("PC_3", pc_data, umap_label), pc("PC_4", pc_data, umap_label), 
                            pc("PC_5", pc_data, umap_label), pc("PC_6", pc_data, umap_label), 
                            pc("PC_7", pc_data, umap_label), pc("PC_8", pc_data, umap_label), 
                            pc("PC_9", pc_data, umap_label), pc("PC_10", pc_data, umap_label), 
                            pc("PC_11", pc_data, umap_label), pc("PC_12", pc_data, umap_label), 
                            pc("PC_13", pc_data, umap_label), pc("PC_14", pc_data, umap_label), 
                            pc("PC_15", pc_data, umap_label), pc("PC_16", pc_data, umap_label), 
                            pc("PC_17", pc_data, umap_label), pc("PC_18", pc_data, umap_label), 
                            pc("PC_19", pc_data, umap_label), pc("PC_20", pc_data, umap_label))

#in PC_1, there is high expression in cluster 13, 7, 1 
#in PC_2, there is a lot of high expression in almost all clusters, except for 5 and 16 

ggsave("pcclusterplot.png",PC_clusterplot, width = 18, height = 18, units = 'in')

#Step 9: Identify conserved cell markers that are present in same cluster and are in all conditions. 
rejoin_obj <- JoinLayers(integrated_obj)

#for clusters 20 and 21, there are not enough cells in all 3 groups so exclude 
markers <- lapply(0:19, conserved_markers, obj = rejoin_obj)
markers <- do.call("rbind", markers)


#extract top 10 markers per cluster 
top10 <- markers %>%  
  mutate(avg_fc = (ZC_avg_log2FC + ZY_avg_log2FC + YF_avg_log2FC)/3) %>%  
           group_by(cluster_id) %>%  
           top_n(n = 10, wt = avg_fc)

write.csv(top10, "top10markerspercluster.csv")



#Step 10: Attempt cell type identification given marker data. 
#manually using marker database: PanglaoDB

#cluster 0 : T cells (CD8A, CD8B, CD3G) 
#cluster 1: Ductal cells (MSMB) 
#cluster 2: T cells (CCR7, ILR7) 
#cluster 3: Acinar cells (GSTA1, GP2) 
#cluster 4: NK Cells (KLRD1, TRDC) 
#cluster 5: Fibroblasts (SFRP2)
#cluster 6: Macrophages/Monocytes (CCL18, LPL, CHIT1)
#cluster 7: Ductal cells (KRT16, CEACAM5) 
#cluster 8: Acinar Cells (AMY2B)
#cluster 9: T Cells (TNFRSF4, CTLA4)
#cluster 10: Dendritic Cells (CD1C, FCER1A) 
#cluster 11: Macrophages/Monocytes (FOLR2, STAB1)
#cluster 12: Macrophages/Monocytes (CLEC4E)
#cluster 13: MK167+ Cells 
#cluster 14: Mast Cells (TPSAB1)
#cluster 15: B-cells (CD79A, VPREB3)  
#cluster 16: Endothelial Cells (RGS5) 
#cluster 17: Endothelial Cells (VWF, CDH5)  
#cluster 18: Plasma cells (IGHV1-24, IGLV2-8) 
#cluster 19: Plasma cells (MND1, SHCBP1) 

#Mentioned before, ZC (normal) has significantly less cells in clusters 0, 2 and 4. 
#This makes sense as the cell types in these three clusters generally have higher
#levels of expression when the cell is cancerous. 


rejoin_obj <- RenameIdents(object = rejoin_obj, 
                                  "0" = "T cells",
                                  "1" = "Ductal cells",
                                  "2" = "T cells",
                                  "3" = "Acinar cells",
                                  "4" = "NK cells",
                                  "5" = "Fibroblasts",
                                  "6" = "Macrophages/Monocytes",
                                  "7" = "Ductal cells",
                                  "8" = "Acinar cells",
                                  "9" = "T cells",
                                  "10" = "Dendritic cells",
                                  "11" = "Macrophages/Monocytes",
                                  "12" = "Macrophages/Monocytes",
                                  "13" = "MK167+ cells",
                                  "14" = "Mast cells",
                                  "15" = "B-cells",
                                  "16" = "Endothelial cells",
                                  "17" = "Endothelial cells", 
                                  "18" = "Plasma cells", 
                                  "19" = "Plasma cells")


DimPlot(object = rejoin_obj, 
        reduction = "umap.cca", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)
#removing 20 and 21 b/c not enough cells for analysis 
seurat_final <- subset(rejoin_obj, idents = c('20', '21'), invert = TRUE) 

seurat_final@meta.data$celltype <- Idents(seurat_final)


final_dimplot <- DimPlot(object = seurat_final, 
        reduction = "umap.cca", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)
ggsave('final_dimplot_withcelltypes.png', final_dimplot, width = 10, height = 14, units= 'in')


#Step 11: Compare the samples across conditions and perform differential expression analysis for specific cell type between conditions

#for ZY vs YF: do psuedobulk analysis , ZC has only 1 sample 
bulk <- AggregateExpression(seurat_final, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("Type",
                                                                                                     "Patient", "celltype"))
tail(Cells(bulk))

#T-cells
t.bulk <- subset(bulk, celltype == "T cells")
Idents(t.bulk) <- "Type"
Tcell.de.bulk <- FindMarkers(t.bulk, ident.1 = "ZY", ident.2 = "YF", slot = "counts", test.use = "DESeq2",
                          verbose = F)
Tcell.de.bulk$gene <- rownames(Tcell.de.bulk)
tcell_de <- ggplot(Tcell.de.bulk, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.05, gene,
                                                                          "")), colour = "red", size = 3)
ggsave("Tcell_DE.png", tcell_de, width = 10, height = 10, units = "in") 

#Ductal cells 
duct.bulk <- subset(bulk, celltype == "Ductal cells")
Idents(duct.bulk) <- "Type"
duct.de.bulk <- FindMarkers(duct.bulk, ident.1 = "ZY", ident.2 = "YF", slot = "counts", test.use = "DESeq2",
                             verbose = F)
duct.de.bulk$gene <- rownames(duct.de.bulk)
duct_de <- ggplot(duct.de.bulk, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.05, gene,
                                                                          "")), colour = "red", size = 3)
ggsave("Ductal_DE.png", duct_de, width = 10, height = 10, units = "in") 

#Acinar cells 
acin.bulk <- subset(bulk, celltype == "Acinar cells")
Idents(acin.bulk) <- "Type"
acin.de.bulk <- FindMarkers(acin.bulk, ident.1 = "ZY", ident.2 = "YF", slot = "counts", test.use = "DESeq2",
                            verbose = F)
acin.de.bulk$gene <- rownames(acin.de.bulk)
acin_de <- ggplot(acin.de.bulk, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.05, gene,
                                                                          "")), colour = "red", size = 3)
ggsave("Acinar_DE.png", acin_de, width = 10, height = 10, units = "in")

#Fibroblasts 
fib.bulk <- subset(bulk, celltype == "Fibroblasts")
Idents(fib.bulk) <- "Type"
fib.de.bulk <- FindMarkers(fib.bulk, ident.1 = "ZY", ident.2 = "YF", slot = "counts", test.use = "DESeq2",
                            verbose = F)
fib.de.bulk$gene <- rownames(fib.de.bulk)
fib_de <- ggplot(fib.de.bulk, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.05, gene,
                                                                          "")), colour = "red", size = 3)
ggsave("Fibroblast_DE.png", fib_de, width = 10, height = 10, units = "in")

#Macrophages/Monocytes
mac.bulk <- subset(bulk, celltype == "Macrophages/Monocytes")
Idents(mac.bulk) <- "Type"
mac.de.bulk <- FindMarkers(mac.bulk, ident.1 = "ZY", ident.2 = "YF", slot = "counts", test.use = "DESeq2",
                           verbose = F)
mac.de.bulk$gene <- rownames(mac.de.bulk)
mac_de <- ggplot(mac.de.bulk, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.05, gene,
                                                                          "")), colour = "red", size = 3)
ggsave("MacrophagesMonocytes_DE.png", mac_de, width = 10, height = 10, units = "in")

#Dendritic cells 
den.bulk <- subset(bulk, celltype == "Dendritic cells")
Idents(den.bulk) <- "Type"
den.de.bulk <- FindMarkers(den.bulk, ident.1 = "ZY", ident.2 = "YF", slot = "counts", test.use = "DESeq2",
                           verbose = F)
den.de.bulk$gene <- rownames(den.de.bulk)
den_de <- ggplot(den.de.bulk, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.05, gene,
                                                                          "")), colour = "red", size = 3)
ggsave("Dendritic_DE.png", den_de, width = 10, height = 10, units = "in")

#MK167+ cells 
mk167.bulk <- subset(bulk, celltype == "MK167+ cells")
Idents(mk167.bulk) <- "Type"
mk167.de.bulk <- FindMarkers(mk167.bulk, ident.1 = "ZY", ident.2 = "YF", slot = "counts", test.use = "DESeq2",
                           verbose = F)
mk167.de.bulk$gene <- rownames(mk167.de.bulk)
mkI67_de <- ggplot(mk167.de.bulk, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.05, gene,
                                                                          "")), colour = "red", size = 3)
ggsave("MKI67_DE.png", mkI67_de, width = 10, height = 10, units = "in")

#Mast cells 
mast.bulk <- subset(bulk, celltype == "Mast cells")
Idents(mast.bulk) <- "Type"
mast.de.bulk <- FindMarkers(mast.bulk, ident.1 = "ZY", ident.2 = "YF", slot = "counts", test.use = "DESeq2",
                             verbose = F)
mast.de.bulk$gene <- rownames(mast.de.bulk)
mast_de <- ggplot(mast.de.bulk, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.05, gene,
                                                                          "")), colour = "red", size = 3)
ggsave("Mast_DE.png", mast_de, width = 10, height = 10, units = "in")

#B-cells 
b.bulk <- subset(bulk, celltype == "B-cells")
Idents(b.bulk) <- "Type"
b.de.bulk <- FindMarkers(b.bulk, ident.1 = "ZY", ident.2 = "YF", slot = "counts", test.use = "DESeq2",
                            verbose = F)
b.de.bulk$gene <- rownames(b.de.bulk)
b_de <- ggplot(b.de.bulk, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.05, gene,
                                                                          "")), colour = "red", size = 3)
ggsave("Bcells_DE.png", b_de, width = 10, height = 10, units = "in")

#Endothelial cells 
endo.bulk <- subset(bulk, celltype == "Endothelial cells")
Idents(endo.bulk) <- "Type"
endo.de.bulk <- FindMarkers(endo.bulk, ident.1 = "ZY", ident.2 = "YF", slot = "counts", test.use = "DESeq2",
                         verbose = F)
endo.de.bulk$gene <- rownames(endo.de.bulk)
endo_de <- ggplot(endo.de.bulk, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.05, gene,
                                                                          "")), colour = "red", size = 3)
ggsave("Endothelial_DE.png", endo_de, width = 10, height = 10, units = "in")


#Plasma cells 
plas.bulk <- subset(bulk, celltype == "Plasma cells")
Idents(plas.bulk) <- "Type"
plas.de.bulk <- FindMarkers(plas.bulk, ident.1 = "ZY", ident.2 = "YF", slot = "counts", test.use = "DESeq2",
                            verbose = F)
plas.de.bulk$gene <- rownames(plas.de.bulk)
plas_de <- ggplot(plas.de.bulk, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.05, gene,
                                                                          "")), colour = "red", size = 3)

ggsave("Plasma_DE.png", plas_de, width = 10, height = 10, units = "in")
