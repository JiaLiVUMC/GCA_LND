## Bulk RNA-seq to validate the LND expansion from control to active CD
## Proportion change of cell types from control to inactive and active CD patients in bulk RNAseq
## The code including one bulk RNA-seq example (Fig3G): deconvolution and proportion change.


library(Seurat)
library(ggplot2)
library(patchwork)
library(reshape2)
library(ggpubr)
library(rstatix)
library(dplyr)

## use Cibersort to characterize cell composition of bulk RNAseq data (GSE179285 TI as example)
## signature matrix construction from endoscopy TI samples.
load(file = "TI/ss_anno.Rdata") ##seurat object with annotion (anno)
Idents(ss) <- ss$anno
marker <- FindAllMarkers(ss, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.7)
marker <- marker[marker$p_val_adj < 0.05,];dim(marker)

t100 <- marker %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
gene <- unique(t100$gene);length(gene)#1469
data <- AverageExpression(ss,group.by = "anno")
data <- data$RNA
sig.data <- data[gene,]
save(sig.data, file = "./sig.data_ti.Rdata")

## Cibersort estimation
#GSE179285 expression and information are downloaded from GEO [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179285]
load(file = "./GSE179285/info.Rdata")
load(file = "./GSE179285/exp.Rdata")
info <- info[,c("geo_accession","source_name_ch1","diagnosis:ch1","inflammation:ch1","tissue:ch1")]
colnames(info) <- c("GSM","Details","Diagnosis","Disease","Location")
info <- info[info$Diagnosis %in% c("Crohn's disease","Healthy control"),]
info <- info[! info$Details %in% "Controls uninflamed sigmoid colon",]
info$Histo <- info$Disease
info$Histo[info$Details %in% c("Controls uninflamed terminal ileum","Controls uninflamed ascending/descending colon")] <- "Normal"

cibersort_fun <- function(exp,location){
  if(location=="TI"){load(file="./sig.data_ti.Rdata")}
  if(location=="AC"){load(file = "./sig.data_ac.Rdata")}
  cibersort.data <- CIBERSORT(sig.data,exp,perm = 100,absolute = F,QN=F)##relative
  return(cibersort.data)
}

##prediction matrix: sample * proportion of each cell type
sample <- info[info$Location %in% "Ileum",]
exp <- exp.gene[,rownames(sample)]
cibersort.data <- cibersort_fun(exp,"TI")

int <- which(colnames(info) == "Histo")
sample <- info[info$Location %in% "Ileum",]
sample$Histo2 <- apply(sample,1,function(x){paste0(x[int],"\nN=",table(sample$Histo)[x[int]])})
input <- as.data.frame(cibersort.data[,1:(ncol(cibersort.data)-3)])
input$Type <- sample[rownames(input),int+1]
compare_list <- split(combn(unique(input[,ncol(input)]),2),  col(combn(unique(input[,ncol(input)]),2)))
level = c("Normal\n(N=8)","Uninflamed\n(N=49)","Inflamed\n(N=33)")
input$Type <- factor(input$Type,levels = level)

##########################################################
#proportion change of LND cluster
#Fig5G (GSE179285 TI as example)

#function to plot figure
bulk.plot.LND <- function(input,title,color,method="fdr",cluster="LND",label="p.adj",position,vjust=-0.3,format=T){
  input <- input[,c(which(colSums(input[,1:ncol(input)-1])>0),ncol(input))]
  tmp.plot <- melt(input,value.name="Prop",variable.name="cluster");
  colnames(tmp.plot)[1] <- "Groups";dim(tmp.plot)
  levels(tmp.plot$Groups)[1] <- gsub("Normal","Control",levels(tmp.plot$Groups)[1])
  
  a <- tmp.plot[tmp.plot$cluster %in% unique(tmp.plot$cluster)[1],]
  b <- a %>% group_by(cluster) %>%pairwise_wilcox_test(Prop ~ Groups)%>%
    adjust_pvalue(method = method) %>% add_significance() %>%  add_xy_position(step.increase = 0.05)
  stat.test <- dplyr::bind_rows(replicate(length(unique(tmp.plot$cluster)),b,simplify = F));dim(stat.test)
  n <- length(unique(tmp.plot$Groups))
  stat.test$cluster <- rep(levels(tmp.plot$cluster),each= (n*(n-1)/2));stat.test$cluster <- factor(stat.test$cluster)
  stat.test$p <- 1
  stat.test$p.adj <- 1
  stat.test$p.adj.signif <- "ns";head(stat.test,10)
  m <- aggregate(tmp.plot$Prop,by=list(tmp.plot$cluster),max)
  rownames(m) <- m$Group.1
  
  for (i in levels(tmp.plot$cluster)) {
    #print(i)
    j <- tmp.plot[tmp.plot$cluster %in% i,]
    if(n==2){
      stat.test$p[stat.test$cluster %in% i] <-  as.character(pairwise.wilcox.test(j$Prop,j$Groups,p.adjust.method = "none")$p.value)
      #stat.test$y.position[stat.test$cluster %in% i] <- c(m[i,2]+(m[i,2]/10))
    }else{
      stat.test$p[stat.test$cluster %in% i] <-  as.character(pairwise.wilcox.test(j$Prop,j$Groups,p.adjust.method = "none")$p.value)[c(1,2,4)]
      #stat.test$y.position[stat.test$cluster %in% i] <- c(m[i,2]+(m[i,2]/10),m[i,2]+(m[i,2]/10)*3,m[i,2]+(m[i,2]/10)*6)
    }
  }
  
  stat.test$group <- paste0(stat.test$group1,stat.test$group2)
  stat.test <- stat.test %>%
    group_by(group) %>%
    adjust_pvalue(method = "fdr") %>% add_significance() 
  
  stat.test$p.adj[is.na(stat.test$p.adj)] <- 1
  stat.test$p.adj.signif[stat.test$p.adj.signif %in% ""] <- "ns"
  
  stat.test <- stat.test[stat.test$cluster %in% cluster,]
  tmp.plot <- tmp.plot[tmp.plot$cluster %in% cluster,]
  
  if(isTRUE(format)){print("Formatted");stat.test$p.adj <- formatC(stat.test$p.adj,format="e",digits=1)}
  p1 <- ggplot(tmp.plot, aes(x=Groups, y=Prop,col=Groups)) + 
    geom_jitter(shape=16,size=1)+
    theme_bw()+labs(x=" ",y=" ",title = title)+theme(
      panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
      plot.title = element_text(size=7,hjust = 0.5,margin = margin(0,0,0,0)),
      legend.position = "none",
      axis.text.x = element_text(color = "black",size = 4.5,margin = margin(2,0,0,0)),
      axis.text.y = element_text(color = "black",size = 5,margin = margin(r = 0)),
      axis.ticks.length = unit(0.05,"cm"),
      axis.ticks.y = element_line(color = "black",size = 0.25),
      axis.title.y = element_blank(),axis.title.x = element_blank(),plot.margin = margin(0,0,5,5))+ 
    scale_color_manual(values=color)+
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))+
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.1) +
    stat_summary(fun=mean, geom="point", color="black");p1
  
  print(stat.test[which(stat.test$p.adj.signif < 0.05),])
  stat.test$y.position[which(stat.test$p.adj.signif < 0.05)] <- position
  p1 <- p1 + stat_pvalue_manual(stat.test, label = label,hide.ns = T,
                                tip.length = 0,size = 1.5,color = "red",vjust = vjust);p1
  return(p1)
}

bulk.plot.LND(input,title = "GSE179285 TI",color = c("#B09C85FF", "#00A087FF", "#800080"),position = c(0.28,0.33))
