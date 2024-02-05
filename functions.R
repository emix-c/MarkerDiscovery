seurat_filtered <- function(obj) { 
  #Store mitochondrial percentage in column 
  obj[['MTpercent']] <- PercentageFeatureSet(obj, pattern = "^MT-")
  #Filter out cells w/ less than 200 genes and >20% MT 
  subset(obj, subset = MTpercent <= 20 & nCount_RNA >= 200 )
}

seurat_prepro<- function(obj){ 
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  RunPCA(obj) }

seurat_pk <- function(obj){ 
  sweep_list <- paramSweep(obj, PCs = 1:20, sct=FALSE) 
  sweep_stats <- summarizeSweep(sweep_list, GT = FALSE)
  #now find optimal pK 
  pKdata <- find.pK(sweep_stats)
  pKdata <- pKdata %>%  
    filter(BCmetric == max(BCmetric)) %>%  
    select(pK)
  as.numeric(as.character(pKdata[[1]]))
}


seurat_df <- function(obj, pK){ 
  annotations <- obj@meta.data$seurat_clusters
  homotypic_est <- modelHomotypic(annotations)
  #in this case 7.5% doublets expected
  exp_doublets <- round(0.076*nrow(obj@meta.data))
  exp_doublets_adj <- round(exp_doublets*(1-homotypic_est))
  #adjusted for homotypic doublets ^ 
  doubletFinder(obj, PCs = 1:20, pN = 0.25, pK = pK, nExp = exp_doublets_adj, 
                reuse.pANN = FALSE) 
}

pc <- function(pc, pc_data, umap_label){
  ggplot(pc_data, aes(umap_1, umap_2)) +
    geom_point(aes(color=.data[[pc]]), alpha = 0.7) +
    scale_color_gradient(guide = 'none' , 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, aes(label=integratedclusters, x, y)) + 
    ggtitle(pc)
}

conserved_markers <- function(cluster, obj){ 
  FindConservedMarkers(obj, ident.1 = cluster, grouping.var = 'Type', only.pos=TRUE) %>%  
    rownames_to_column(var = 'gene') %>%  
    left_join(y = unique(grch38[, c('symbol', 'description')]), by = c('gene' = 'symbol')) %>%  
    cbind(cluster_id = cluster,.)
  
} 

