library(Seurat)
library(ggplot2)
library(patchwork)
library(reshape2)
library(ggpubr)
library(rstatix)
library(dplyr)

bulk.plot.LND <- function(input,title,color,method="fdr",cluster="LND",label="p.adj.signif",position){
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
                                tip.length = 0,size = 2.5,color = "red");p1
  return(p1)
}

##Example
#endoscopy TI signature matrix construction
load(file = "TI/ss_anno.Rdata") ##seurat object with all TI endoscopy samples.
ss$anno[ss$anno %in% c("Early LND","Late LND")] <- "LND"
Idents(ss) <- ss$anno
marker <- FindAllMarkers(ss, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.7)
marker <- marker[marker$p_val_adj < 0.05,];dim(marker)

t100 <- marker %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
gene <- unique(t100$gene);length(gene)#1469
data <- AverageExpression(ss,group.by = "anno")
data <- data$RNA
sig.data <- data[gene,]

#GSE179285 TI
load(file = "/GSE179285/ti_plot.Rdata")
levels(input$Type) <- paste0(gsub("\n","\n(",levels(input$Type)),")")
bulk.plot.LND(input,title = "GSE179285 TI",color = c("#B09C85FF", "#00A087FF", "#DC0000FF"),position = c(0.28,0.33))
