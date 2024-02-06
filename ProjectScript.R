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

#Step 10: Attempt cell type identification given marker data. 
#manually using marker database: PanglaoDB

#cluster 0 : T cells (CD8A, CD8B, CD3G)
#cluster 1:  Cholangiocytes 
#cluster 2: T cells (CCR7, ILR7)
#cluster 3: Acinar cells (GSTA1, GP2)
#cluster 4: Gamma Delta T Cells (KLRD1)
#cluster 5: Fibroblasts (SFRP2)
#cluster 6: Macrophages (CCL18, LPL)
#cluster 7: Epithelial cells (KRT16, CEACAM5)
#cluster 8: Acinar Cells (AMY2B)
#cluster 9: Gamma Delta T Cells (LAIR2)
#cluster 10: Dendritic Cells (CLEC10A, FCER1A)
#cluster 11: Dendritic Cells (FOLR2, STAB1)
#cluster 12: Dendritic Cells (CLEC4E)
#cluster 13: Germ cells (PBK)
#cluster 14: Mast Cells (TPSAB1)
#cluster 15: B-cells (CD79A)
#cluster 16: Pancreatic stellac cells (RGS5)
#cluster 17: Endothelial Cells (VWF, CDH5)
#cluster 18: B-cells (IGHV1-24, IGLV2-8)
#cluster 19: Germ cells (CMTM2)

#Mentioned before, ZC (normal) has significantly less cells in clusters 0, 2 and 4. 
#This makes sense as the cell types in these three clusters generally have higher
#levels of expression when the cell is cancerous. 


rejoin_obj <- RenameIdents(object = rejoin_obj, 
                                  "0" = "T cells",
                                  "1" = "Cholangiocytes",
                                  "2" = "T cells",
                                  "3" = "Acinar cells",
                                  "4" = "Gamma Delta T cells",
                                  "5" = "Fibroblasts",
                                  "6" = "Macrophages",
                                  "7" = "Epithelial cells",
                                  "8" = "Acinar cells",
                                  "9" = "Gamma Delta T cells",
                                  "10" = "Dendritic cells",
                                  "11" = "Dendritic cells",
                                  "12" = "Dendritic cells",
                                  "13" = "Germ cells",
                                  "14" = "Mast cells",
                                  "15" = "B-cells",
                                  "16" = "Pancreatic stellac cells",
                                  "17" = "Endothelial cells", 
                                  "18" = "B-cells", 
                                  "19" = "Germ cells")

DimPlot(object = rejoin_obj, 
        reduction = "umap.cca", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)
#removing 20 and 21 b/c not enough cells for analysis 
seurat_final <- subset(rejoin_obj, idents = c('20', '21'), invert = TRUE) 



#DoHeatmap(object = rejoin_obj, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

final_dimplot <- DimPlot(object = seurat_final, 
        reduction = "umap.cca", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)
ggsave('final_dimplot_withcelltypes.png', final_dimplot, width = 10, height = 14, units= 'in')



