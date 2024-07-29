## ligand-receptor analysis with cellchat v1.4.0
## The code generates figures for fig5A

library(CellChat)
library(ggarchery)

#function to run cellchat 
#parameter int: the column number of annotation used to run cellchat in Seurat object meta data 
cc.fun <- function(int,group,definedGroup=NULL,pair=NULL,trim=NULL,key,genes){
  meta.sub <- meta[meta[,int] %in% group,]
  if(!is.null(definedGroup)){meta.sub <- meta.sub[meta.sub$anno %in% definedGroup,]}else{
    meta.sub <- meta.sub[meta.sub$anno %in% names(table(meta.sub$anno)[table(meta.sub$anno)>9]),]
  }
  data.sub <- data[genes,rownames(meta.sub)]
  cellchat <- createCellChat(object = data.sub, meta = meta.sub, group.by = "anno")
  if(!is.null(pair)){
    CellChatDB.use <- subsetDB(CellChatDB.human,search = pair, key = key)
  }else{CellChatDB.use <- CellChatDB.human}
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  if(!is.null(trim)){
    cellchat <- computeCommunProb(cellchat,type = "truncatedMean",trim = trim,nboot = 50)
  }else{
    cellchat <- computeCommunProb(cellchat,do.fast = T,nboot = 50)#type=triMean; default
  }
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  return(cellchat)
}

##run cellchat (TI as example)
load(file = "./TI/ss_anno.Rdata")
data <- GetAssayData(ss,slot = "data")
meta <- ss@meta.data

ref <- createCellChat(object = data, meta = meta, group.by = "anno")
CellChatDB.use <- CellChatDB.human
ref@DB <- CellChatDB.use
ref <- subsetData(ref)
dim(ref@data.signaling)# 878 82749
genes.all <- rownames(ref@data.signaling)#878 the numebr of genes in within the ligand-receptor database

#run ALL pathway
#int: the annotion column
#genes: use genes within the ligand-receptor database in cellchat to save time.
#trim=0: no trim
cellchat <- cc.fun(int=20,group = unique(meta[,20]),
                   definedGroup=NULL,
                   pair=NULL,trim=0,key="pathway_name",genes=genes.all)
save(cellchat,file = "./TI/Allsample_Allpathway_trim0.Rdata")

##Fig5A incoming and outcoming strength of each celltype in TI

#function to access strength
strengh.plot <- function(){
  centr <- cellchat@netP$centr
  outgoing <- matrix(0, nrow = nlevels(cellchat@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(cellchat@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(cellchat@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[,i] <- centr[[i]][["outdeg"]]
    incoming[,i] <- centr[[i]][["indeg"]]
  }
  signaling <- path.use[path.use %in% cellchat@netP$pathways]
  outgoing <- outgoing[ , signaling, drop = FALSE]
  incoming <- incoming[ , signaling, drop = FALSE]
  outgoing.cells <- rowSums(outgoing)
  incoming.cells <- rowSums(incoming)
  
  num.link <- aggregateNet(cellchat, signaling = signaling, return.object = FALSE, remove.isolate = FALSE)$count
  num.link <- rowSums(num.link) + colSums(num.link)-diag(num.link)
  df <- data.frame(x = outgoing.cells, y = incoming.cells, labels = names(incoming.cells),
                   Count = num.link)
  df$labels <- factor(df$labels, levels = names(incoming.cells))
  df$Group <- ct.group[names(incoming.cells)]
  return(df)
}

load(file = "./TI/Allsample_Allpathway_trim0.Rdata")
cellchat@DB$interaction["PLG_F2R","interaction_name_2"] <- "PLF_F2R"
sort(unique(cellchat@idents))

#cell types group to Immune, Stromal and Epithelial
ct.imm <- c("B Proliferating","B","Plasma","Neutrophils",
            #            "CTL/NK","T Proliferating","T","Mast","Endo.","Fibro.");length(ct.imm)
            "NK","T Proliferating","T","Resident Macro.","Recruited Macro.","Mast","Endo.","Fibro.");length(ct.imm)
ct.epi <- c("Stem","Cycling TA","Early Enterocyte","Intermediate Enterocyte","Mature Enterocyte","LND",
            "Goblet Proliferating","Early Goblet","Mature Goblet","BEST4/OTOP2","EEC","Tuft","Paneth");length(ct.epi)
ct <- c(ct.imm,ct.epi)
ct.group <- c(rep("Immune",length(ct.imm)-2),rep("Stromal",2),rep("Epithelial",length(ct.epi)))
names(ct.group) <- ct

pathway<-cellchat@netP$pathways
pathwaylist<-grep("^CXCL|CCL|^IL",pathway);pathway[pathwaylist]
path.use <- unique(c(pathway[pathwaylist],"CXCL","CCL","SAA","TNF","IFN-II",
                     "IL10","IL2","TRAIL","LIGHT","BAFF","APRIL","VEGI"))

data <- strengh.plot()
ggplot(data = data, aes(x, y)) +
  geom_point(aes(size = Count, colour = labels, fill = labels, shape = Group))+theme_bw()+
  labs(title = "TI", x = "Outgoing (Source) strength", y = "Incoding (Target) strength")+ theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    axis.text.x = element_text(size = 6,color = "black"),
    axis.text.y = element_text(size = 6,color = "black"),
    axis.title = element_text(size = 8,color = "black"),
    legend.text = element_text(colour = "black",size = 6),
    legend.title = element_text(size = 7),
    legend.key.height = unit(0.15,"in"),legend.key.width = unit(0.05,"in"),
    legend.margin = margin(0,0,0,0),legend.box.margin = margin(0,0,0,-5),plot.margin = margin(2,5,2,5))+
  scale_fill_manual(values = cc[levels(cellchat@idents),2])+ guides(fill= "none")+
  scale_colour_manual(values = cc[levels(cellchat@idents),2])+ guides(colour= "none")+
  guides(size = guide_legend(order = 0), shape = guide_legend(order = 1))+
  scale_shape_manual(values = c(21,24,22))+
  scale_size_continuous(range = c(1,3))
