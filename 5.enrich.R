
library(clusterProfiler) #进行富集分析的基础R包
library(org.Rn.eg.db) #大鼠注释信息的R包
library(org.Hs.eg.db) #人注释信息的R包
library(org.Mm.eg.db) #小鼠注释信息的R包
library(ggplot2)#可视化
#如果有没有安装使用BiocManager::install()安装
BiocManager::install('clusterProfiler')



DEG.id <- mapIds(x=org.Rn.eg.db,   #这里选好你需要的物种
                 keys=rownames(DEG), #这里输入你的差异基因
                 keytype ='SYMBOL',
                 column='ENTREZID')
DEG.id <- na.omit(DEG.id)  #去除没有匹配到的基因
KEGG <- enrichKEGG(gene = DEG.id,
                   organism = "rno", #这里选择物种 人：hsa ；大鼠：rno ；小鼠：mmu
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 1)
dotplot(KEGG, showCategory = 15,color='pvalue')  #clusterProfiler自带可视化
GO <- enrichGO(gene = DEG.id,
               OrgDb = org.Rn.eg.db,  #这里选好你需要的物种
               pvalueCutoff =0.05,
               ont="all", #这里可以选择BP CC MF ，选择all就是三个项目全部分析
               readable =T,
               qvalueCutoff = 1)
dotplot(GO,showCategory = 5,color='pvalue',split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')


enrich_go <- list()
enrich_kegg <- list()
for (i in as.numeric(unique(top100$cluster))){
  small_gene_group=top100[top100$cluster==(i-1),]
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
  KEGG4 <- enrichKEGG(gene = unique(df_name$ENTREZID),
                      organism = "mmu",
                      pvalueCutoff =0.05,
                      qvalueCutoff = 0.2)
  
  GO4 <- enrichGO(gene = unique(df_name$ENTREZID),
                  OrgDb = org.Mm.eg.db,
                  pvalueCutoff =0.05,
                  ont="all",
                  readable =T,
                  qvalueCutoff = 0.2)
  enrich_go[[i]] <- GO
  enrich_kegg[[i]] <- KEGG
  # go_res=go@result
  # if (dim(go_res)[1] != 0) {
  #   go_res$cluster=i
  #   allcluster_go=rbind(allcluster_go,go_res)
  # }
}

