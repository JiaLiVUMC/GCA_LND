## Proportion change of cell types from control to inactive and active CD patients
## The code generates figures for fig2C,2D; fig3E,3F

library(Seurat)
library(ggplot2)
library(patchwork)
library(reshape2)
library(ggpubr)
library(rstatix)
library(readxl)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

#This function generate figures with selected cell type.
prop.plot <- function(input,title,color,hide.ns,method,legendLocation,stat.test,ct,strip.size=8,ncol=7,plot.margin=margin(0,0,0,2),title.y,vjust=-0.3,format=F,mult=c(0, 0.15)){
  input <- cbind(input[,ct],input$`sampleInfo[rownames(prop_ct), 2]`)
  tmp.plot <- melt(input,value.name="Prop",variable.name="cluster");
  colnames(tmp.plot)[1] <- "Groups";dim(tmp.plot)
  
  m <- aggregate(tmp.plot$Prop,by=list(tmp.plot$cluster),max)
  rownames(m) <- m$Group.1
  n <- length(unique(tmp.plot$Groups))
  
  stat.test <- stat.test[stat.test$cluster %in% ct,]
  stat.test$cluster <- factor(stat.test$cluster,levels = ct)
  
  for (i in levels(tmp.plot$cluster)) {
    if(n==2){
      stat.test$y.position[stat.test$cluster %in% i] <- c(m[i,2]+(m[i,2]/10))
    }else{
      stat.test$y.position[stat.test$cluster %in% i] <- c(m[i,2]+(m[i,2]/10)*3,m[i,2]+(m[i,2]/10),m[i,2]+(m[i,2]/10)*6)
    }
  }
  
  stat.test$y.position <- as.numeric(stat.test$y.position)
  
  if(isTRUE(format)){print("Formatted");stat.test$p.adj <- formatC(stat.test$p.adj,format="e",digits=1)}
  p1 <- ggplot(tmp.plot, aes(x=Groups, y=Prop,col=Groups)) + 
    geom_jitter(shape=16,size=1)+facet_wrap(~cluster,scales="free",ncol = ncol)+
    theme_bw()+labs(x=" ",y=title.y,title = title)+theme(
      panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
      panel.spacing.y = unit(0, "lines"),
      #plot.title = element_text(size=8,hjust = 0.5,face = "bold",margin = margin(0,0,0,0)),
      plot.title = element_blank(),
      strip.text = element_text(size=strip.size,margin = margin(0,0,2,0)),
      strip.background = element_blank(),
      legend.position = legendLocation,legend.text = element_text(colour = "black",size = 6),
      legend.title = element_blank(),legend.key.size = unit(0.5, "lines"),
      legend.margin = margin(0,0,0,0),legend.box.margin = margin(0,0,0,-5),legend.spacing.x = unit(0.1,"cm"),
      axis.text.x = element_blank(),axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black",size = 5,margin = margin(r = 0)),
      axis.ticks.y = element_line(color = "black",size = 0.25),
      axis.ticks.length.y = unit(0.05,"cm"),
      axis.title.y = element_text(color = "black",size = 8),
      axis.title.x = element_blank(),plot.margin = plot.margin)+ 
    scale_color_manual(values=color)+
    scale_y_continuous(expand = expansion(mult = mult))+
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.1) +
    stat_summary(fun=mean, geom="point", color="black")
  
  p1 <- p1 + stat_pvalue_manual(stat.test, label = "p.adj",hide.ns = hide.ns,
                                tip.length = 0,size = 1.5,color = "red",vjust = vjust);p1
  return(p1)
}

#This function generate proportion of each cell type for each patient.
prop.fun <- function(bottom,ct.epi){
  if(bottom=="all") {
    meta <- ss@meta.data
    sampleInfo <- unique(data.frame(sample=meta$PatientID,histo3=meta$histo3))
    rownames(sampleInfo) <-sampleInfo$sample
    prop_ct <- table(meta$PatientID,meta$anno)
    prop_ct<- as.data.frame.array(prop_ct/rowSums(prop_ct))
    prop_ct <- prop_ct[,ct.epi]
    input <- cbind(prop_ct,sampleInfo[rownames(prop_ct),2])
  }
  return(input)
}

##TI endoscopy; Fig2C and Fig3E
# input is seurat object with three group histology information (histo3); patient ID (PatientID); annotation (anno)
load(file = "/TI/ss_anno.Rdata") ##seurat object with all TI endoscopy samples.
ss$histo3 <- factor(ss$histo3,levels = c("Normal","inactive_CD","active_CD"))
sample <- unique(ss@meta.data[,c("PatientID","histo3")]);table(sample$histo3)
levels(ss$histo3) <- paste0(levels(ss$histo3)," (N=",table(sample$histo3)[levels(ss$histo3)],")")

input <- prop.fun(bottom = "all",ct.epi = unique(ss$anno))
tmp.plot <- melt(input,value.name="Prop",variable.name="cluster");
colnames(tmp.plot)[1] <- "Groups";dim(tmp.plot)

#wilcoxon test for each histology pair
stat.test <- tmp.plot %>%
  group_by(cluster) %>%
  pairwise_wilcox_test(Prop~ Groups) %>% add_xy_position(x = "Groups")

#FDR adjustment for each histology pair
stat.test$group <- paste0(stat.test$group1,stat.test$group2)
stat.test <- stat.test %>%
  group_by(group) %>%
  adjust_pvalue(method = "fdr") %>% add_significance() 

#Fig2C; TI immune and stromal cell type
ct.imm<- c("B Proliferating","B","Plasma","Resident Macro.","Recruited Macro.","Neutrophils","CTL/NK","T Proliferating","T","Mast","Endo.","Fibro.");length(ct.imm)
prop.plot(input,title = " ",color = c("#B09C85FF", "#00A087FF", "#800080"),hide.ns = T,method = "fdr",legendLocation = "right",stat.test = stat.test,ct=ct.imm,strip.size = 6,ncol = 6,title.y = "Proportion (TI)",format = TRUE)

#Fig3E; TI Epithelial cell type 
ct.epi <- c("Stem","Cycling TA","Early Enterocyte","Intermediate Enterocyte","Mature Enterocyte","LND","Goblet Proliferating","Early Goblet","Mature Goblet","BEST4/OTOP2","EEC","Tuft","Paneth");length(ct.epi)
prop.plot(input,title = " ",color = c("#B09C85FF", "#00A087FF", "#800080"),hide.ns = T,method = "fdr",legendLocation = c(.94,.2),stat.test = stat.test,ct=ct.epi,plot.margin = margin(0,2,0,2),strip.size = 6,title.y = "Proportion (TI)",format = TRUE)


##AC endoscopy; Fig2D and Fig3F
# input is seurat object with three group histology information (histo3); patient ID (PatientID); annotation (anno)
load(file = "/AC/ss_anno.Rdata") ##seurat object with all AC endoscopy samples.
ss$histo3 <- factor(ss$histo3,levels = c("Normal","inactive_CD","active_CD"))
sample <- unique(ss@meta.data[,c("PatientID","histo3")]);table(sample$histo3)
levels(ss$histo3) <- paste0(levels(ss$histo3)," (N=",table(sample$histo3)[levels(ss$histo3)],")")
#levels(ss$histo3) <- c("Control (N=14)","Inactive CD (N=32)","Active CD (N=18)")

input <- prop.fun(bottom = "all",ct.epi = unique(ss$anno))
tmp.plot <- melt(input,value.name="Prop",variable.name="cluster");
colnames(tmp.plot)[1] <- "Groups";dim(tmp.plot)

#wilcoxon test for each histology pair
stat.test <- tmp.plot %>%
  group_by(cluster) %>%
  pairwise_wilcox_test(Prop~ Groups) %>% add_xy_position(x = "Groups")

#FDR adjustment for each histology pair
stat.test$group <- paste0(stat.test$group1,stat.test$group2)
stat.test <- stat.test %>%
  group_by(group) %>%
  adjust_pvalue(method = "fdr") %>% add_significance() 

#Fig2C; AC immune and stromal cell type
ct.imm<- c("B Proliferating","B","Plasma","Resident Macro.","Recruited Macro.","Neutrophils","CTL/NK","T Proliferating","T","Mast","Endo.","Fibro.");length(ct.imm)
prop.plot(input,title = " ",color = c("#B09C85FF", "#00A087FF", "#800080"),hide.ns = T,method = "fdr",legendLocation = "right",stat.test = stat.test,ct=ct.imm,strip.size = 6,ncol = 6,title.y = "Proportion (AC)",format = TRUE)

#Fig3F; AC Epithelial cell type 
ct.epi <- c("Stem","Cycling TA","Early Colonocyte","Intermediate Colonocyte","Mature Colonocyte","LND","Goblet Proliferating","Early Goblet","Mature Goblet","BEST4/OTOP2","EEC","Tuft");length(ct.epi)
prop.plot(input,title = " ",color = c("#B09C85FF", "#00A087FF", "#800080"),hide.ns = T,method = "fdr",legendLocation = c(.94,.2),stat.test = stat.test,ct=ct.epi,plot.margin = margin(0,2,0,2),strip.size = 6,title.y = "Proportion (AC)",format = TRUE)



