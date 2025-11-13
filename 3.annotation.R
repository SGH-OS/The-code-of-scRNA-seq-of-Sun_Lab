#annotation----
cell_marker <- read.table('scmarker.txt',header = T,sep = '\t')
DotPlot(scRNA,features = cell_marker$Marker)+coord_flip()

#correlation----
library(corrplot)
library(RColorBrewer)
scRNA$subtype <- paste0("Cluster_",scRNA$seurat_clusters)

Idents(scRNA) <- "Sample"
Markers <- FindAllMarkers(scRNA,only.pos= T, min.pct = 0.2, logfc.threshold = 0.25,verbose = T,max.cells.per.ident = 2000)
Markers <- subset(Markers[grep("^RP[L|S]",Markers$gene, ignore.case = FALSE,invert=TRUE),],subset=p_val_adj < 0.05)
Markers_av <- AverageExpression(scRNA,group.by = "Sample",features = unique(Markers$gene),assays = "RNA") 
Markers_av <- Markers_av$RNA
gene_cell_exp <- t(scale(t(Markers_av),scale = T,center = T))

#计算相关性
cell_cor <- cor(as.matrix(gene_cell_exp), method = 'spearman')

#显著性检验
cell_cor_p <- cor.mtest(cell_cor, conf.level = 0.95)$p
corrplot(cell_cor, 
         method = 'circle',#'square'
         order = "hclust", 
         hclust.method='ward.D',
         type = "full", 
         #addrect = 5,
         col = rev(brewer.pal(n=8, name="RdYlBu")),
         tl.cex=0.8,
         tl.col = "black",
         cl.cex = 0.5,
         cl.pos = "r",
         cl.length = 5,
         cl.ratio=0.1,
         mar=c(3, 0, 2, 0)) 

#细胞注释
#----神经（人）----
a <- c('SNAP25','MAP2','RBFOX3','SYP',
       #Astrocyte
       'AQP4','NTSR2','HTRA1',
       #OPC
       'PLP1','MOBP','MAG','MOG',
       #ODC
       'GPR17','SOX10','PDGFRA',
       #Microglia
       'CTSS','CX3CR1','LY86','AIF1',
       #Endothelial
       'CLDN5','FLT1','TEK','CD34','PECAM1','PROM1',
       #Pericyte
       'PDGFRB','VTN','MYL9',
       #Ependyma
       'FOXJ1','SOX2','AK7','RSPH1',
       #Stromal
       'DCN','APOD','GSN','COL1A1','COL3A1',
       #Erythocyte
       
       #Leukocyte
       'MS4A4B','LTB','CTSW','CD3E',
       #Neutrophil
       'TREM1','S100A8','S100A9')
#神经（小鼠）-----
a <- c('Snap25','Map2','Rbfox3','Syp',
       #Astrocyte
       'Aqp4','Ntsr2','Htra1',
       #OPC
       'Plp1','Mobp','Mag','Mog',
       #ODC
       'Gpr17','Sox10','Pdgfra',
       #Microglia
       'Ctss','Cx3cr1','Ly86','Aif1',
       #Endothelial
       'Cldn5','Flt1','Tek','Cd34','Pecam1','Prom1',
       #Pericyte
       'Pdgfrb','Vtn','Myl9',
       #Ependyma
       'Foxj1','Sox2','Ak7','Rsph1',
       #Stromal
       'Dcn','Apod','Gsn','Col1a1','Col3a1',
       #Erythocyte
       
       #Leukocyte
       'Ms4a4b','Ltb','Ctsw','Cd3e',
       #Neutrophil
       'Trem1','S100a8','S100a9')
#皮肤-----
marker <- c('SELE','TM4SF1','PECAM1',#Endothelial_cells
            'COL1A1','COL1A2','COL3A1',#Fibroblasts
            'TAGLN','ACTA2','TPM2',#Smooth_muscle_cells
            'KRT1','KRT5','KRT10','KRT14',#Keratinocytes
            'LYZ','HLA-DRA',#Immune_cells
            'CCL21','LYVE1',#Lymphatic_endothelial_cells
            'NRXN1',#Neural_cells
            'SCGB1B2P','SCGB1D2',#Sweat_gland_cells
            'TYRP1','PMEL' )#melanocytes

marker <- c('Sele','Tm4sf1','Pecam1',#Endothelial_cells
            'Col1a1','Col1a2','Col3a1',#Fibroblasts
            'Tagln','Acta2','Tpm2',#Smooth_muscle_cells
            'Krt1','Krt5','Krt10','Krt14',#Keratinocytes
            'Lyz','Hla-dra',#Immune_cells
            'Ccl21','Lyve1',#Lymphatic_endothelial_cells
            'Nrxn1',#Neural_cells
            'Scgb1b2p','Scgb1d2',#Sweat_gland_cells
            'Tyrp1','Pmel' )#melanocytes
####-----
DotPlot(scRNA,features = marker)+coord_flip()
cell_type <- c("0" = "",
               "1" = "",
               "2" = "",
               "3" = "",
               "4" = "",
               "5" = "",
               "6" = "",
               "7" = "",
               "8" = "",
               "9" = "",
               "10" = "",
               "11" = "",
               "12" = "",
               "13" = "",
               "14" = "",
               "15" = "")

scRNA[['cell_type']]<-unname(cell_type[scRNA@meta.data$seurat_clusters])
DimPlot(scRNA,reduction = 'umap',group.by = 'cell_type',label = T)
