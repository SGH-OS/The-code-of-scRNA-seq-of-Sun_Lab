
#CCA\harmony----

scRNA[["RNA"]] <- split(scRNA[["RNA"]], f = scRNA$GSE)
#CCA
seurat_merge_v5 <- IntegrateLayers(  object = scRNA, method = CCAIntegration,  orig.reduction = "pca", new.reduction = "cca",  verbose = FALSE)
#RPCA
seurat_merge_v5 <- IntegrateLayers(  object = seurat_merge_v5, method = RPCAIntegration,  orig.reduction = "pca", new.reduction = "rpca",  verbose = FALSE)
#harmony
seurat_merge_v5 <- IntegrateLayers(object = scRNA, method = HarmonyIntegration,  orig.reduction = "pca", new.reduction = "harmony",  verbose = FALSE)
#remotes::install_github("satijalab/seurat-wrappers")
#BiocManager::install('batchelor')
#devtools::install_github("satijalab/seurat-data", force = TRUE)
library(SeuratData)
library(batchelor)
library(SeuratWrappers)
#FastMNN
seurat_merge_v5 <- IntegrateLayers(  object = seurat_merge_v5, method = FastMNNIntegration,  new.reduction = "mnn",  verbose = FALSE)


#####CCA可视化######
seurat_merge_v5 <- FindNeighbors(seurat_merge_v5, reduction = "cca", dims = 1:15)
seurat_merge_v5 <- FindClusters(seurat_merge_v5, resolution = 0.5)
seurat_merge_v5 <- RunUMAP(seurat_merge_v5, reduction = "cca",                
                           dims = 1:15)

p5 <- DimPlot(seurat_merge_v5,  
              reduction = "umap.cca",  
              group.by = c("cca_clusters"))+  
  DimPlot(seurat_merge_v5,    
          reduction = "umap.cca",    
          group.by = c("orig.ident"))

#####RPCA可视化######
seurat_merge_v5 <- FindNeighbors(seurat_merge_v5, reduction = "rpca", dims = 1:20)
seurat_merge_v5 <- FindClusters(seurat_merge_v5, resolution = 0.1, cluster.name = "rpca_clusters")
seurat_merge_v5 <- RunUMAP(seurat_merge_v5, reduction = "rpca",                            
                           dims = 1:20,                            
                           reduction.name = "umap.rpca")
p6 <- DimPlot(  seurat_merge_v5,  reduction = "umap.rpca",  group.by = c("rpca_clusters"))+  DimPlot(    seurat_merge_v5,    reduction = "umap.rpca",    group.by = c("orig.ident"))

#####harmony可视化######
seurat_merge_v5 <- FindNeighbors(seurat_merge_v5, reduction = "harmony", dims = 1:20)
seurat_merge_v5 <- FindClusters(seurat_merge_v5, resolution = 0.5)
seurat_merge_v5 <- RunUMAP(seurat_merge_v5, reduction = "harmony",
                           dims = 1:20)
p7 <- DimPlot(  seurat_merge_v5,  reduction = "umap.harmony",  group.by = c("harmony_clusters"))+  DimPlot(    seurat_merge_v5,    reduction = "umap.harmony",    group.by = c("orig.ident"))

#####FastMNN可视化########
seurat_merge_v5 <- FindNeighbors(seurat_merge_v5, reduction = "mnn", dims = 1:20)
seurat_merge_v5 <- FindClusters(seurat_merge_v5, resolution = 0.1, cluster.name = "mnn_clusters")
seurat_merge_v5 <- RunUMAP(seurat_merge_v5, reduction = "mnn",                            
                           dims = 1:20,                            
                           reduction.name = "umap.mnn")
p8 <- DimPlot(  seurat_merge_v5,  reduction = "umap.mnn",  group.by = c("mnn_clusters"))+  DimPlot(    seurat_merge_v5,    reduction = "umap.mnn",    group.by = c("orig.ident"))













