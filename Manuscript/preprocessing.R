## preprocess single cell RNA-seq data
## generate proper Seurat object from GEO data [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266546]
## Seurat version v4.1.0

#Endoscopy TI as example
#download gene-cell count matrix: GSE266546_Raw_counts_TI_Endoscopy.txt.gz
#download meta data with annotation, histology information: GSE266546_Annotation_TI_Endoscopy.txt

#create Seurat object
counts <- data.table::fread("GSE266546_Raw_counts_TI_Endoscopy.txt.gz")
rownames(counts) <- counts$V1
counts <- counts[,-1]
meta <- read.table(file = "GSE266546_Annotation_TI_Endoscopy.txt",header = T,row.names = 1,sep = "\t");dim(meta)

ss <- CreateSeuratObject(counts = counts,meta.data = meta)
ss$histo3 <- ss$Histology
ss$histo3 <- ifelse(ss$histo3 %in% c("Normal_CD","Quiescent"),"Inactive CD",ss$histo3)
ss$histo3 <- ifelse(ss$histo3 %in% c("Mild","Moderate","Severe"),"Active CD",ss$histo3)
ss$anno <- ss$Celltypes

save(ss,file = "ss_anno.Rdata")