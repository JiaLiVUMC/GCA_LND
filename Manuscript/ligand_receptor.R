library(CellChat)
library(ggarchery)

pair.plot <- function(pairLR){
  net <- cellchat@net
  pairLR.name <- intersect(rownames(pairLR), dimnames(net$prob)[[3]])
  pairLR <- pairLR[pairLR.name, ]
  prob <- net$prob
  pval <- net$pval
  prob[pval > 0.05] <- 0
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name], 3, sum) != 0]
  }else {
    pairLR.name.use <- pairLR.name[sum(prob[, , pairLR.name]) != 0]
  }
  pairLR <- pairLR[pairLR.name.use, ]
  nRow <- length(pairLR.name.use)
  prob <- prob[, , pairLR.name.use]
  pval <- pval[, , pairLR.name.use]
  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify = "array")
    pval <- replicate(1, pval, simplify = "array")
  }
  prob.sum <- apply(prob, c(1, 2), sum)
  return(prob.sum)
}

circle.plot <- function(prob.sum,top=0.005,title,curvature){
  #top=0.005
  thresh <- stats::quantile(prob.sum, probs = 1-top)
  prob.sum[prob.sum < thresh] <- 0
  prob.sum[is.na(prob.sum)] <- 0
  g <- graph_from_adjacency_matrix(prob.sum, mode = "directed", weighted = T);g
  edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE);edge.start
  edge.start <- edge.start[!edge.start[,1] == edge.start[,2],];edge.start
  
  t <- seq(0, 2*pi, length.out = (nrow(prob.sum)+1))
  x <- cos(t);x
  y <- sin(t);y
  df <- data.frame(t, x, y)
  df <- head(df,-1);df
  rownames(df) <- ct
  
  #data <- strengh.plot()
  data <- df
  data$group <- ct.group[rownames(data)]
  data$ct <- rownames(data)
  data$ct <- factor(data$ct,levels = data$ct)
  data <- data[rownames(prob.sum),]
  
  arr <- data.frame(x=data[edge.start[,1],2],y=data[edge.start[,1],3],
                    xend=data[edge.start[,2],2],yend=data[edge.start[,2],3],
                    value=scales::rescale(prob.sum,to=c(0,2))[edge.start],
                    color=cc[rownames(data)[edge.start[,1]],2])
  
  p <- ggplot(data) + 
    geom_point(aes(x=x, y=y,colour = ct, fill = ct, shape = group,label=ct),size=3) + coord_fixed()+
    labs(title = title, x = "", y = "")+
    theme_void()+theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5,size = 10,face = "bold",margin = margin(0,0,15,0)),
      axis.text = element_blank(),
    )+
    scale_fill_manual(values = cc[levels(data$ct),2])+
    scale_colour_manual(values = cc[levels(data$ct),2])+
    scale_shape_manual(values = c(21,24,22))+
    #geom_text_repel(data[as.numeric(unique(as.character(edge.start))),],mapping=aes(x=x,y=y,label = ct, colour = ct))+
    geom_curve(data = arr,aes(x=x,y=y,xend=xend,yend=yend),size=arr$value,arrow = arrow(),
               lineend = "round",linejoin = "round",curvature = curvature,angle = 135,color=arr$color,
               position = position_attractsegment(start_shave = 0.1, end_shave = 0.1, type_shave = "distance"))
  
  
  #geom_arrowsegment(data = arr,aes(x=x,y=y,xend=xend,yend=yend,size=value),color=arr$color,
  #                 arrow_positions = 0.7,lineend = "round",linejoin = "round",
  #           position = position_attractsegment(start_shave = 0.1, end_shave = 0.1, type_shave = "distance"))
  return(p)
}

##Example
#Fig4C TI; LND as source; neutrophil recruitment
db <- CellChatDB.human$interaction
db <- db[db$receptor %in% c("CXCR1","CXCR2"),]
db <- db[rownames(db) %in% dimnames(prob)[[3]],]
db <- db[,1:4]
prob.sum <- pair.plot(db)
circle.plot(prob.sum,title = "Neutrophil Recruitment",curvature = 0.3,top = 0.01)

