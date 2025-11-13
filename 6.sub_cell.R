
library(Seurat)
library(tidyverse)

load('Result/after_anno.rds')
table(scRNA$cell_type)
i <- scRNA[,scRNA@meta.data$cell_type=='Neuron']
i <- scRNA[,scRNA@meta.data$cell_type=='Oligodendrocytes']
i <- scRNA[,scRNA@meta.data$cell_type=='Microglia']
i <- scRNA[,scRNA@meta.data$cell_type=='Astrocytes']
i <- scRNA[,scRNA@meta.data$cell_type=='Endothelial']
i <- scRNA[,scRNA@meta.data$cell_type=='Pericyte']
i <- scRNA[,scRNA@meta.data$cell_type=='OPC']
i <- scRNA[,scRNA@meta.data$cell_type=='Ependyma']
i <- scRNA[,scRNA@meta.data$cell_type=='Stromal']

ElbowPlot(i,ndims = 50)
i[["RNA"]] <- split(i[["RNA"]], f = i$Tech)
i <- NormalizeData(i,scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = 'vst',nfeatures = 2000) %>% 
  ScaleData() %>% #vars.to.regress = "percent.mt"
  RunPCA(features = VariableFeatures(object = i)) 
i <- IntegrateLayers(object = i, method = CCAIntegration,  orig.reduction = "pca", new.reduction = "cca",  verbose = FALSE)
i[["RNA"]] <- JoinLayers(i[["RNA"]])
i <- FindNeighbors(i,dims=1:15,reduction = "cca") %>% 
  FindClusters(resolution = 0.1) %>% 
  RunUMAP(dims=1:15,reduction = "cca")

#i <- FindClusters(i,resolution = 0.1)
p1 <- DimPlot(i,reduction = 'umap',label = T)
ggsave(plot = p1,filename = './Result/3细胞亚型/9Stromal/1umap.pdf',width = 6, height = 5,create.dir = T)
sc_deg <- FindAllMarkers(i,logfc.threshold = 2,only.pos = T)
sc_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top100
write.csv(top100,file = './Result/3细胞亚型/1Neuron/cell_marker.csv',quote = F)

library(ggalluvial)
Cellratio <- prop.table(table(i$seurat_clusters, i@meta.data$Group), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
ggplot(Cellratio,aes(x =Var2, y= Freq, fill = Var1,stratum = Var1,alluvium = Var1)) + 
  geom_bar(stat = "identity",width = 0.7,size = 0.5,colour = 'white')+ 
  theme_classic() +
  #geom_col(width = 0.2, color='black')+
  geom_flow(width=0.7,alpha=0.3, knot.pos=0.1,colour = 'white')+
  labs(x='Group',y = 'Ratio',fill='Cell_cluster')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))#+
scale_fill_manual(values = brewer.pal(8, "Set3"))

ggsave('./Result/3细胞亚型/1Neuron/2cellratio.pdf',width =7,height = 3,create.dir = T)

Idents(i) <- i$Group
sc_deg <- FindAllMarkers(i,logfc.threshold = 0.2,only.pos = T)
sc_deg %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top100

enrich_go <- list()
enrich_kegg <- list()
for (x in as.numeric(unique(top100$cluster))){
  small_gene_group=top100[top100$cluster==names(table(top100$cluster))[x],]
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
  KEGG <- enrichKEGG(gene = unique(df_name$ENTREZID),
                      organism = "mmu",
                      pvalueCutoff =0.5,
                      qvalueCutoff = 1)
  KEGG <- setReadable(KEGG,OrgDb = "org.Mm.eg.db",keyType = "ENTREZID")
  GO <- enrichGO(gene = unique(df_name$ENTREZID),
                  OrgDb = org.Mm.eg.db,
                  pvalueCutoff =0.5,
                  ont="all",
                  readable =T,
                  qvalueCutoff = 1)
  enrich_go[[x]] <- GO
  enrich_kegg[[x]] <- KEGG
  write.csv(KEGG@result,file = paste0('./Result/3细胞亚型/9Stromal/KEGG_',names(table(top100$cluster))[x],'.csv'),quote = F)
  write.csv(GO@result,file = paste0('./Result/3细胞亚型/9Stromal/GO_',names(table(top100$cluster))[x],'.csv'),quote = F)
}


Neuron <- i
save(Neuron,file = 'Neuron.rds')


edox2 <- pairwise_termsim(GO)
treeplot(edox2,)
library(emapplot)
