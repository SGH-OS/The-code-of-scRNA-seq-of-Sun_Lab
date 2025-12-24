
devtools::install_github("Japrin/STARTRAC")
library(Startrac)
library(pheatmap)

dat.file <- system.file("extdata/example.cloneDat.Zhang2018.txt", package = "Startrac")
in.dat <- Tcell@meta.data
in.dat

# 计算OR比值比时可改为method = "fisher"
R_oe <- calTissueDist(in.dat,
                      byPatient = F,
                      colname.cluster = "cell_type4",
                      colname.patient = "Sample",
                      colname.tissue = "Group",
                      method = "chisq", 
                      min.rowSum = 0) 
R_oe

pheatmap(R_oe,display_numbers = T,scale ='row',color = colorRampPalette(c("#1E90FF", 'white',"#FB8072"))(100))

anno_matrix <-ifelse(
  R_oe >1, "+++",
  ifelse(R_oe >0.8& R_oe <=1, "++",
         ifelse(R_oe >0.2& R_oe <=0.8, "+",
                ifelse(R_oe >0& R_oe <=0.2, "+/−", "−")))
)
max_abs_deviation <-max(abs(R_oe -1), na.rm =TRUE)
breaks <-seq(1- max_abs_deviation, 1+ max_abs_deviation, 
             length.out =101)

pheatmap(R_oe,scale ='row',
         color = colorRampPalette(c("#1E90FF", 'white',"#FB8072"))(100),
         #breaks = breaks,legend_breaks =c(-4,-2,0,1,2,4,6),#c(fivenum(Roe),1)
         display_numbers = anno_matrix  # 传递标注矩阵，显示符号而非数值   
)
pdf(file='./Result/5roe/3heatmap.pdf', width=5, height=5)
pheatmap(R_oe,display_numbers = T,scale ='row',color = colorRampPalette(c("#1E90FF", 'white',"#FB8072"))(100))
dev.off()


