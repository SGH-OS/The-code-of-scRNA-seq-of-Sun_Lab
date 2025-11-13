
#visualization
devtools::install_github("zhanghao-njmu/SCP")
#UMAP----
cell <- table(scRNA$cell_type)
annotation_df <- data.frame(
  label = paste0(names(cell), " (n=", cell,')'),
  x = 10,  # 统一x坐标
  y = seq(                          # 动态生成y坐标
    from = 13, 
    by = -2,                     # 控制行间距
    length.out = length(cell)
  )
)

DimPlot(scRNA,reduction = 'umap',group.by = 'cell_type',label = T)+
  geom_text(
    data = annotation_df,
    aes(x = x, y = y, label = label),
    color = "black",          # 统一标签颜色
    hjust = 0,                # 左对齐
    inherit.aes = FALSE       # 避免继承color映射
  )

DimPlot(scRNA,reduction = 'umap',group.by = 'cell_type',label = T)+
  scale_color_discrete(labels=annotation_df$label)

DimPlot(scRNA,reduction = 'umap',group.by = 'cell_type',label = T)+
  theme_bw()+
  #theme(aspect.ratio = 1) +  #图片比例
  theme(panel.grid=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.title = element_text(size = 0))+
  xlab('UMAP_1')+
  ylab('UMAP_2')

library(tidydr)
DimPlot(scRNA,reduction = 'umap',group.by = 'cell_type',label = T)+# 添加reducation="tsne" 改为tsne
  NoAxes() + # 隐藏默认坐标轴标签和刻度线
  labs(x ="UMAP1", y ="UMAP2") + # labs 函数是ggplot2 
  theme_dr() + # 应用带小箭头的坐标轴主题（来自tidydr包）
  theme(panel.grid = element_blank() ) # 去除所有网格线

#dotplot----
scRNA <- JoinLayers(scRNA)
sc_deg <- FindAllMarkers(scRNA,logfc.threshold = 2,only.pos = T)

sc_deg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

library(RColorBrewer)
DotPlot(scRNA,features =c('YAP1'),group.by = 'cell_type')+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  coord_flip()+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )

#vlnplot

VlnPlot(scRNA,features = c('YAP1','SELE','TM4SF1','PECAM1'),pt.size = 0,stack = T)&
  scale_fill_manual(values = c(brewer.pal(11,"Set3"),brewer.pal(3,"Dark2")))&
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    strip.text.x = element_text(angle = 45,size = 10,hjust = 0,vjust = 0,face = 'plain'),
    legend.position = "none"
  )

VlnPlot(scRNA,features = c('YAP1','SELE','TM4SF1','PECAM1'),pt.size = 0)&theme(
  axis.title.x = element_blank(),
  #axis.title.y = element_blank(),
  #axis.ticks.y = element_blank()
)& scale_y_continuous(limits = c(0,5))

geom_boxplot(width=0.2,col="black",fill="white")

a <- M@meta.data
library(ggpubr)
ggviolin(a,x = 'b',y = 'Efferocytosis_score1',fill = 'b',
         add = 'boxplot',
         add.params = list(fill='white'))+
  stat_compare_means(comparisons = com)+theme(legend.position = 'none')

com=list(c('DFU','NC'))

#多个
featuresVlnPlot(human_data, features = c("ANXA1","S100A8"), group.by = "group")&  
  theme_bw()&  
  theme(axis.title.x = element_blank(),        
        axis.text.x = element_text(color = 'black',face = "bold", size = 12),        
        axis.text.y = element_text(color = 'black', face = "bold"),        
        axis.title.y = element_text(color = 'black', face = "bold", size = 15),        
        panel.grid.major = element_blank(),        
        panel.grid.minor = element_blank(),        
        panel.border = element_rect(color="black",size = 1.2, linetype="solid"),        
        panel.spacing = unit(0.12, "cm"),        
        plot.title = element_text(hjust = 0.5, face = "bold.italic"),        
        legend.position = 'none')&  
  stat_compare_means(method="t.test",hide.ns = F,                      
                     comparisons = my_comparisons,                     
                     label="p.signif",                    
                     bracket.size=0.8,                     
                     tip.length=0,                     
                     size=6)&  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))&   
  scale_fill_manual(values = c("#FF5744","#208A42"))


#featureplot----
FeaturePlot(scRNA,
            features=c('YAP1','SELE','TM4SF1','PECAM1'),
            raster=F,#细胞过多时候需要加这个参数
            max.cutoff = 3
)&
  scale_color_viridis_c()&
  theme_bw()&
  theme(panel.grid=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())&
  xlab('UMAP_1')&
  ylab('UMAP_2')


FeaturePlot(scRNA,features = c("YAP1"),
            ,pt.size = 0.5,ncol = 3,cols = c("#ccccca", "#e61a2e"),order = T)&
  theme_bw()&
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 20)
  )

#heatmap----
DoHeatmap(scRNA,
          features = c('YAP1','SELE','TM4SF1','PECAM1'),
          group.by = "cell_type",
          assay = 'RNA',
          label = F)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))


library(pheatmap)
colanno=mye.seu@meta.data[,c("CB","celltype")]
colanno=colanno%>%arrange(celltype)
rownames(colanno)=colanno$CB
colanno$CB=NULL
colanno$celltype=factor(colanno$celltype,levels = unique(colanno$celltype))
rowanno=markerdf1
rowanno=rowanno%>%arrange(celltype)

mat4=mye.seu[["RNA"]]@scale.data[rowanno$gene,rownames(colanno)]
mat4[mat4>=2.5]=2.5
mat4[mat4 < (-1.5)]= -1.5 #小于负数时，加括号！

pheatmap(mat4,cluster_rows = F,cluster_cols = F,
         show_colnames = F,
         annotation_col = colanno,
         gaps_row=as.numeric(cumsum(table(rowanno$celltype))[-6]),
         gaps_col=as.numeric(cumsum(table(colanno$celltype))[-6]),
         filename="heatmap.2.pdf",width=11,height = 7
)

library(ComplexHeatmap)
Heatmap(df_scaled,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 8),
        col = col_fun,
        rect_gp = gpar(col="white"),
        column_split = 4,
        column_title = NULL,
        cluster_rows = TRUE,
        show_row_dend = FALSE,
        top_annotation =group,
        right_annotation = lab2,
        show_row_names = FALSE,
        name = "Exp")

#cellratio----
library(ggalluvial)
Cellratio <- prop.table(table(scRNA$cell_type, scRNA@meta.data$Group), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
ggplot(Cellratio,aes(x =Var2, y= Freq, fill = Var1,stratum = Var1,alluvium = Var1)) + 
  geom_bar(stat = "identity",width = 0.7,size = 0.5,colour = 'white')+ 
  theme_classic() +
  #geom_col(width = 0.2, color='black')+
  geom_flow(width=0.7,alpha=0.3, knot.pos=0.1,colour = 'white')+
  labs(x='Group',y = 'Ratio',fill='Cell_type')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))#+
scale_fill_manual(values = brewer.pal(8, "Set3"))

library(RColorBrewer)
i <- scRNA[,scRNA$Group=='Con']
cell_counts <- table(i$cell_type)
df <- as.data.frame(cell_counts)

# 计算比例
pielabel <- paste0(df$Var1, " (", round(df$Freq / sum(df$Freq) * 100, 2), "%)")

# 绘制2D饼图
pie(df$Freq,
    labels = pielabel,
    col = brewer.pal(length(df$Freq), "Set3"),
    main = "Cell Type Proportion of Con")
ggsave('./Result/2降维注释/10pie_RG3.pdf',width =7,height = 5,create.dir = T)


