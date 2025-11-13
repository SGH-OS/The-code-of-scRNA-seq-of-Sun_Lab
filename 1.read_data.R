###########################
#                         #
#     Writen by hora      #
#                         #
###########################

library(Seurat)
library(tidyverse)
library(DoubletFinder)
#Read data
getwd()
files <- basename(dir("/Users/hora/data/"))
file_name <- dput(files)
sc_list <- list()
for (i in 1:length(file_name)) {
  scRNA_data <- Read10X(data.dir = str_c('/Users/hora/data/', file_name[i]))
  sc <- CreateSeuratObject(counts = scRNA_data,
                           min.cells = 3, min.features = 200)
  sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")#"^MT-" "^mt-"
  #HB.genes <- c("Hba-a1","Hba-a2","Hbb-bt",'Hbb-bs',"Hbq1a",'Hbm', 'Hba1', 'Hba2', 'Hbb', 'Hbd', 'Hbe1', 'Hbg1', 'Hbg2', 'Hbq1','Hbz')
  HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB_m <- match(HB.genes, rownames(sc@assays$RNA)) 
  HB.genes <- rownames(sc@assays$RNA)[HB_m] 
  HB.genes <- HB.genes[!is.na(HB.genes)] 
  sc[["percent.HB"]]<-PercentageFeatureSet(sc, pattern =HB.genes) 
  sc[["percent.ribo"]] <- PercentageFeatureSet(sc, pattern = "^RP[SL]")#"^RP[SL]" "^Rp[sl]"
  sc <- subset(sc, 
               subset = nFeature_RNA > 200 & 
                 nFeature_RNA < 5000 & 
                 #nCount_RNA > 800 & 
                 nCount_RNA < 40000 & 
                 percent.mt < 10)
  sc <- NormalizeData(sc,scale.factor = 10000)
  sc <- FindVariableFeatures(sc,selection.method = 'vst',nfeatures = 2000)
  sc <- ScaleData(sc, features = rownames(sc))
  sc <- RunPCA(sc,features = VariableFeatures(object = sc))
  sc <- FindNeighbors(sc,dims = 1:15)
  sc <- FindClusters(sc,resolution = 2)
  sc <- RunUMAP(sc,dims = 1:15)
  #去除双细胞
  sweep.res.list <- paramSweep(sc, PCs = 1:15, sct = F,num.cores = 6) # 扫描不同pK值的效果
  sweep.stats_object <- summarizeSweep(sweep.res.list, GT = FALSE)    # 汇总扫描结果
  bcmvn <- find.pK(sweep.stats_object)                                # 找到最佳pK（BCmetric最高点）
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])         # 提取最佳pK数值
  gc() 
  annotations <- sc@meta.data$seurat_clusters                         # 获取当前分群信息
  homotypic.prop <- modelHomotypic(annotations)                       # 计算同源双细胞比例(相同类型细胞形成的双细胞)
  nExp_poi <- round(0.06 * nrow(sc@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  sc <- doubletFinder(sc, 
                      PCs = 1:15,         # 使用前15个主成分
                      pN = 0.25,          # 生成25%伪双细胞用于训练
                      pK = pK_bcmvn,           # 注入最佳灵敏度参数
                      nExp = nExp_poi.adj,    # 初步预测的双细胞数
                      #reuse.pANN = FALSE, # 首次运行需生成pANN矩阵
                      sct = F)
  
  sc <- subset(sc, cells = colnames(sc)[sc@meta.data[length(colnames(sc@meta.data))] == "Singlet"])
  sc_list[[i]] <- sc
}
names(sc_list) <- file_name
scRNA <- merge(x = sc_list[[1]],
               y = sc_list[c(2:length(sc_list))],
               add.cell.ids = names(sc_list))
scRNA@meta.data$Sample <- substr(rownames(scRNA@meta.data),1,4)
scRNA@meta.data$Group <- substr(rownames(scRNA@meta.data),1,5)

VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = 'Group',pt.size = 0)

plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# reduction----
scRNA <- JoinLayers(scRNA)
scRNA <- NormalizeData(scRNA,scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA,selection.method = 'vst',nfeatures = 2000)
top10 <- head(VariableFeatures(scRNA),10)
plot1 <- VariableFeaturePlot(scRNA)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

scRNA <- ScaleData(scRNA)# vars.to.regress = "percent.mt"
scRNA <- RunPCA(scRNA,features = VariableFeatures(object = scRNA))
PCAPlot(scRNA,group.by="Group")
scRNA <- CellCycleScoring(
  scRNA,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes)
PCAPlot(scRNA,group.by="Phase")

ElbowPlot(scRNA,ndims = 50)

scRNA <- FindNeighbors(scRNA,dims = 1:20)
scRNA <- FindClusters(scRNA,resolution = 0.5)
scRNA <- RunUMAP(scRNA,dims = 1:20)
DimPlot(scRNA,reduction = 'umap',label = T)


