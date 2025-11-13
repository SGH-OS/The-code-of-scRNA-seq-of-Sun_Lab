library(monocle)
Tcell <- JoinLayers(Tcell)

sc <- Treg
sc <- JoinLayers(sc)
data <- LayerData(sc, assay="RNA", layer='counts')
ge <- data.frame(gene_id=rownames(data),gene_short_name=rownames(data))
rownames(ge) <- ge$gene_short_name
pd=new("AnnotatedDataFrame",data = sc@meta.data)
fd=new("AnnotatedDataFrame",data = ge)

test=newCellDataSet(as(data,'sparseMatrix'),phenoData = pd,featureData = fd)
#大数据集使用稀疏矩阵，节省内存，加快运算
test <- estimateSizeFactors(test) 
test <- estimateDispersions(test)
test=detectGenes(test,min_expr = 0.1) #计算每个基因在多少细胞中表达

#Idents(Myeloid) <- Myeloid@meta.data$cell_type3
markers1 <- FindAllMarkers(sc, only.pos = TRUE,
                          logfc.threshold = 0.5,
                          test.use = "wilcox")
markers1 %>%
  group_by(cluster) %>%
  top_n(n = 200, wt = avg_log2FC) -> top10

test_ordering_genes=unique(top10$gene)
#test_ordering_genes=rownames(top2000)
test=setOrderingFilter(test,ordering_genes = test_ordering_genes) 
#指明哪些基因用于后续的聚类/排序

test=reduceDimension(test,reduction_method = "DDRTree",max_components = 2, norm_method = 'log',residualModelFormulaStr = "~num_genes_expressed") 
#h
#residualModelFormulaStr减少其他因素的影响，比如不同样本、不同批次
#test=orderCells(test,root_state = 5)
test=orderCells(test)
#绘图
#
p1 <- plot_cell_trajectory(test,color_by = 'Pseudotime')
ggsave(plot = p1,filename = './Result/4monocle/1Pseudotime.pdf',width = 5, height = 4,create.dir = T)

p1 <- plot_cell_trajectory(test,color_by = 'State')
ggsave(plot = p1,filename = './Result/4monocle/1State.pdf',width = 5, height = 4,create.dir = T)

p1 <- plot_cell_trajectory(test,color_by = 'GSE')
ggsave(plot = p1,filename = './Result/4monocle/2cell_type.pdf',width = 5, height = 4,create.dir = T)
p1 <- plot_cell_trajectory(test,color_by = 'cell_type')+facet_wrap(~cell_type,nrow=1)
ggsave(plot = p1,filename = './Result/4monocle/2cell_type2.pdf',width = 10, height = 4,create.dir = T)

p1 <- plot_cell_trajectory(test,color_by = "Pseudotime")+
  #scale_color_manual(values = brewer.pal(6,'Set3'))+
  facet_wrap(~Group,nrow=2)
ggsave(plot = p1,filename = './Result/4monocle/2Pseudotime_group.pdf',width = 10, height = 7,create.dir = T)

plot_cell_trajectory(test,color_by = 'cell_type')
plot_cell_trajectory(test,color_by = 'Group')
plot_cell_trajectory(test,color_by = 'State')
plot_cell_trajectory(test,color_by = 'TTN')
plot_cell_trajectory(test,color_by = "RNA_snn_res.0.2")+scale_color_manual(values = brewer.pal(6,'Set3'))


plot_cell_trajectory(test,color_by = "Pseudotime")+
  #scale_color_manual(values = brewer.pal(6,'Set3'))+
  facet_wrap(~Sample,nrow=4)

plot_cell_trajectory(test,color_by = "RNA_snn_res.0.3")+
  facet_wrap(~RNA_snn_res.0.3,nrow=1)+
  scale_color_manual(values = brewer.pal(5, "Set2"))
plot_cell_trajectory(test,color_by = "dataname")+
  facet_wrap(~dataname,nrow=1)+
  scale_color_manual(values = brewer.pal(4, "Set2"))

plot_genes_in_pseudotime(test[c('TTN','FOXP3','CTLA4'),])
plot_genes_branched_pseudotime(test[c('TTN','FOXP3','CTLA4'),])
plot_genes_branched_heatmap(test[c('TTN','FOXP3','CTLA4'),])

ggsave("celltype.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))

BEAM_res=BEAM(test,branch_point = 1,cores = 6,progenitor_method=c('sequential_split','duplicate')[1])
expressed_genes=row.names(subset(fData(test),num_cells_expressed>=10)) #在部分基因里面找
pseudotime_de <- differentialGeneTest(test[expressed_genes,],cores = 6,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
pseudotime_de <- pseudotime_de[order(pseudotime_de$qval), ]

BEAM_res <- BEAM_res[order(BEAM_res$qval),]
a <- plot_genes_branched_heatmap(test[row.names(subset(head(BEAM_res,20))),],
                                 branch_point = 1,
                                 num_clusters = 3, #这些基因被分成几个group
                                 cores = 4,
                                 # branch_labels = c("Cell fate 1", "Cell fate 2"),
                                 # #hmcols = NULL, #默认值
                                 # hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                 # branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T #是否返回一些重要信息
)
a <- plot_genes_branched_heatmap(test[rownames(head(BEAM_res,300)),],
                            branch_point = 1,
                            num_clusters = 3,
                            show_rownames = T,
                            return_heatmap = T)
ggsave(plot = a$ph_res,filename = './Result/4monocle/4BEAM.pdf',width = 6, height = 6,create.dir = T)

gene_group=a$annotation_row

plot_genes_branched_pseudotime(test[c('TTN'),],
                               branch_point = 2,
                               color_by = 'cluster')

gene_group=a$annotation_row
gene_group$gene=rownames(gene_group)

library(clusterProfiler)
library(org.Hs.eg.db)
allcluster_go=data.frame()
allcluster_kegg=data.frame()
for (i in unique(gene_group$Cluster)) {
  small_gene_group=filter(gene_group,gene_group$Cluster==i)
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  go <- enrichGO(gene         = unique(df_name$ENTREZID),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "all",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 1,
                 readable      = TRUE)
  go_res=go@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    allcluster_go=rbind(allcluster_go,go_res)
  }
}

for (i in unique(gene_group$Cluster)) {
  small_gene_group=filter(gene_group,gene_group$Cluster==i)
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  kegg <- enrichKEGG(gene = unique(df_name$ENTREZID),
                     organism = "hsa", #这里选择物种 人：hsa ；大鼠：rno ；小鼠：mmu
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 1)
  kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  go_res=kegg@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    allcluster_kegg=rbind(allcluster_kegg,go_res)
  }
}

write.csv(allcluster_go,file = './Result/4monocle/GO.csv',quote = F)
write.csv(allcluster_kegg,file = './Result/4monocle/KEGG.csv',quote = F)
a$ph_res


rep <- as.data.frame(pData(test))
a <- as.data.frame(data)
rep$TTN <- as.numeric(a['TTN',])

ggplot(rep, aes(TTN, fill = State)) +geom_density() +facet_wrap(~State,nrow =5, scales = "free_y") +
  theme_bw() +RotatedAxis() +
  theme(
    strip.text = element_blank(),
    strip.background = element_rect(color = "white", fill = "white"),
    panel.grid = element_blank()
  ) 

save(Tcell,test,file = 'Treg.rds')

test@phenoData@data$TTN <- as.numeric(a['TTN',])
