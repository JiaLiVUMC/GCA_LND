## Trajectory analysis of two LND subtypes in TI
## The code generate figures for Fig4A and Fig4D


##Fig4A
#RNA velocity analysis for single sample #GCA062 TI (#6145-YX-2)
#input sub: Seurat object with annotion (anno)
#input ldat: loom file including spliced/unspliced matrix generated from the bam file.

library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

load(file = "6145-YX-2.Rdata")#Seurat object
ldat <- ReadVelocity(file = "6145-YX-2_S1_L005.loom")#loom file

cn <- colnames(ldat$spliced)
cn <- gsub(":","_",cn)
cn2 <- as.data.frame(matrix(unlist(strsplit(cn,"_")),ncol = 5,byrow = T))
cn2$id <- paste(cn2$V1,cn2$V5,sep = "_")
for (i in 1:length(ldat)) {
  colnames(ldat[[i]]) <- cn2$id
  ldat[[i]] <- ldat[[i]][,cn2$id %in% colnames(sub)]
  ldat[[i]] <- ldat[[i]][,match(colnames(sub),colnames(ldat[[i]]))]
}

bm0 <- as.Seurat(x = ldat)
bm <- sub
bm[["spliced"]] <- bm0@assays$spliced
bm[["unspliced"]] <- bm0@assays$unspliced
Idents(bm) <- bm$anno
bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02,verbose = F)

#generate figure
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
names(x = ident.colors) <- levels(x = bm)
cell.colors <- ident.colors[Idents(object = bm)]
names(x = cell.colors) <- colnames(x = bm)
show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), vel = Tool(object = bm, slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), cex = 0.8, arrow.scale = 3, show.grid.flow = T, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, do.par = FALSE, cell.border.alpha = 0.1)


##Fig4D 
#CytoTRACE stemness prediction for single sample #GCA062 TI (#6145-YX-2)
#input sub: Seurat object 

load(file = "6145-YX-2.Rdata")#Seurat object
Idents(sub) <- sub$anno
sub <- subset(sub,idents=c("Stem","Cycling TA","Early Enterocyte","Early LND","Intermediate Enterocyte",
                           "Late LND","Mature Enterocyte"))
cyto.exp <- as.matrix(GetAssayData(sub,slot = "count"))
cyto.cy <- CytoTRACE::CytoTRACE(cyto.exp,ncores = 4)
sub$cytotrace <- cyto.cy$CytoTRACE

tmp <- data.frame(Groups=Idents(sub),Score=sub$cytotrace)
ct <- c("Stem","Cycling TA","Early Enterocyte","Early LND","Intermediate Enterocyte","Late LND","Mature Enterocyte")
tmp$Groups <- factor(tmp$Groups,levels = ct)
levels(tmp$Groups) <- paste0(levels(tmp$Groups),"\n(N=",table(tmp$Groups)[levels(tmp$Groups)],")")
tmp <- melt(tmp,value.name="Score",variable.name="cluster");head(tmp,2)

#p value
stat.test <- tmp %>%
  group_by(cluster) %>%
  pairwise_wilcox_test(Score~ Groups) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>% add_xy_position(x = "Groups",step.increase = 0.25);stat.test

p <- ggplot(tmp,aes(x=Groups, y=Score)) + geom_violin(aes(fill=Groups,col=Groups))+
  labs(x=" ",y="Cytotrace Stemness",title = "GCA062 TI")+
  theme_bw()+theme(
    legend.position = "none",
    plot.title = element_blank(),panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30,hjust = 1,colour = "black",size = 8,margin = margin(0,0,0,0)),
    axis.text.y = element_text(color = "black",size = 6),
    axis.title.y = element_text(color = "black",size = 8),plot.margin = margin(5,0,0,0))+
  scale_fill_manual(values = cc[ct,2])+
  scale_color_manual(values = cc[ct,2])+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.1) +
  stat_summary(fun=mean, geom="point", color="black");p

stat.test$p.adj <- formatC(stat.test$p.adj,format="e",digits=1)
p <- p + stat_pvalue_manual(stat.test, label = "p.adj",hide.ns = T,tip.length = 0,size = 1.5,color = "red",vjust=0.2)