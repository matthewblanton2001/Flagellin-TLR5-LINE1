library(DESeq2)

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE164704", "file=GSE164704_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")

# load gene annotations 
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID

# sample selection
gsms <- "0XX1XX0XX1XX0XX1XX0X1XX0XX1XX"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
tbl <- tbl[ ,sel]

# group membership for samples
gs <- factor(sml)
groups <- make.names(c("Control","Flagellin"))
levels(gs) <- groups
sample_info <- data.frame(Group = gs, row.names = colnames(tbl))

ds<-DESeqDataSetFromMatrix(countData=tbl,colData=sample_info,design=~Group)
ds<-DESeq(ds)
r<-results(ds,contrast=c("Group",groups[2],groups[1]),alpha=0.05,pAdjustMethod="fdr")
tT<-merge(as.data.frame(r),annot,by=0,sort=F)
tT<-tT[-c(1:6342,6344:9416,9418,9420:9459,9461:18221,18223:20043,20045:30694,30696:31345,31347:39376),]
rownames(tT)<-tT$Symbol
gene_ids <- tT$Row.names
mats <- counts(ds, normalized=TRUE)[gene_ids, ]
rownames(mats) <- tT$Symbol
mat.z<-t(apply(mats,1,scale))
col_fun<-circlize::colorRamp2(c(-2,0,2),c("blue","white","red"))
colnames(mat.z)<-sample_info$Group
col_vec=c("red","blue")
col_vec=stats::setNames(col_vec,levels(sample_info$Group))
col1 <- list(Group=c(Control="blue",Flagellin="red")) 
column_order=rownames(mat.z) 
column_split=factor(as.character(sample_info$Group),levels=c("Flagellin","Control")) 
Heatmap(mat.z,col=col_fun,cluster_rows=FALSE,cluster_columns=FALSE,name="Zscore",column_split=column_split,show_column_names=FALSE,row_names_side="left")
