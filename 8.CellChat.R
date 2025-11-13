#cellchat----
library(CellChat)
library(patchwork)
AT.str <- GetAssayData(object = seurat_object,
                       assay = "RNA",
                       slot = "data")
AT.meta.data <- seurat_object@meta.data
cell.ctl = rownames(AT.meta.data)[AT.meta.data$Group == "Con"]
cell.ser = rownames(AT.meta.data)[AT.meta.data$Group == "Ser"]
cell.RG = rownames(AT.meta.data)[AT.meta.data$Group == "RG3"]
#star
data.input = AT.str[, cell.ctl]
meta = AT.meta.data[cell.ctl, ]
data.input = AT.str[, cell.ser]
meta = AT.meta.data[cell.ser, ]
data.input = AT.str[, cell.RG]
meta = AT.meta.data[cell.RG, ]
#meta$cell_type3 <- paste0('Mic',meta$cell_type2)
cellchat <- createCellChat(object = data.input, meta =meta, group.by = "cell_type")

levels(cellchat@idents)
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
str(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB,
                           search = c("Secreted Signaling",'ECM-Receptor','Cell-Cell Contact'))#"Secreted Signaling",'ECM-Receptor','Cell-Cell Contact'
future::plan("multicore", workers = 4)
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
#dplyr::glimpse(CellChatDB$interaction)AT.cellchat@DB <- CellChatDB.use
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

#1
par(mfrow = c(1,2), xpd=TRUE)
#mat <- cellchat@net$count
mat <- cellchat@net$weight
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[, 3] <- mat[,3]
netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[1])

par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#2
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_scatter(cellchat)

#3
netVisual_bubble(cellchat, remove.isolate = T)

netVisual_bubble(cellchat, signaling = c('CCL','PTN'),sources.use = c(5,6,8), targets.use = c(5,6,8), remove.isolate = F)

#4
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",signaling = c("PTN" ,"NCAM","NRXN","PSAP","MK","NEGR","CADM","LAMININ","NRG","JAM","CNTN","CLDN","TENASCIN","MAG","FN1","CCL","COLLAGEN","PDGF","NGL","CDH","VEGF","APP","PTPRM","VTN","ESAM","GAS"    ,"CDH5","IGF"))
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",signaling = c("PTN" ,"NCAM","NRXN","PSAP","MK","NEGR","CADM","LAMININ","NRG","JAM","CNTN","CLDN","TENASCIN","MAG","FN1","CCL","COLLAGEN","PDGF","NGL","CDH","VEGF","APP","PTPRM","VTN","ESAM","GAS"    ,"CDH5","IGF"))
ht1 + ht2

library(RColorBrewer)
mat <- cellchat@netP$pattern$outgoing
n_unique <- length(unique(as.numeric(mat)))
my_color <- colorRampPalette(brewer.pal(9, "BuGn"))(n_unique + 1)
my_color <- colorRampPalette(brewer.pal(9, "BuGn"))(20)
netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "outgoing",
  color.heatmap = my_color
)

ht <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", color.use = colors)

netAnalysis_signalingRole_network(cellchat,signaling = c('CCL','PTN'))

netAnalysis_signalingRole_network(cellchat,
                                  signaling = pathways.show,
                                  width = 8, height = 2.5,
                                  font.size = 10)

cellchat <- cellchat_ser
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = 'CXCL', layout = "chord")
netVisual_bubble(cellchat, sources.use = c(5,6,8), targets.use = c(5,6,8), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(6),remove.isolate = FALSE)
netVisual_bubble(cellchat,  targets.use = c(6), remove.isolate = FALSE)
netVisual_bubble(cellchat,  remove.isolate = T)
plotGeneExpression(cellchat, signaling = "CCL")

cellchat1 <- cellchat

data(CellChatDB.human)
all_interactions_mouse <- CellChatDB.human$interaction
all_pathways_mouse <- unique(all_interactions_mouse$pathway_name)
print(all_pathways_mouse)

pheatmap::pheatmap(b,cluster_rows = F,cluster_cols = F,color = colorRampPalette(brewer.pal(9, "BuGn"))(100),border_color = NA)
