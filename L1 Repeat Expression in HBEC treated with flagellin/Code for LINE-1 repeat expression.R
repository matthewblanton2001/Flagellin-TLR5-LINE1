# Load Output from HOMER and clinical data file
raw_counts <-read.csv("outputfile.csv",stringsAsFactors=FALSE)
clindata <- read.csv("clindata.csv",stringsAsFactors=FALSE)
rownames(clindata) <- clindata$X 
rownames(raw_counts) <- raw_counts$Transcript 
raw_counts <- raw_counts[-c(1:8)] 
colnames(raw_counts) <- c("Flagellin Treated 1","Flagellin Treated 2", "Flagellin Treated 3", "Flagellin Treated 4", "Untreated 1","Untreated 2","Untreated 3","Untreated 4") 
all(rownames(clindata)%in%colnames(raw_counts)) 
all(rownames(clindata)==colnames(raw_counts)) 

# Run DESeq2 and apply Z-score
dds <- DESeqDataSetFromMatrix(countData=raw_counts,colData=clindata,design=~Group) 
dds <- DESeq(dds) 
resultsNames(dds) 
res <- results(dds,contrast=c("Group","Flagellin","Untreated")) 
res <- as.data.frame(res) 
mats <- counts(dds,normalized=TRUE)[rownames(res),] 
mat.z <- t(apply(mats,1,scale))
colnames(mat.z) <- clindata$X

# ComplexHeatmap Visualization
col_vec=c("red","blue")
col_vec=stats::setNames(col_vec,levels(clindata$Group))
col1 <- list(Group=c(Untreated="blue",Flagellin="red")) 
column_order=rownames(mat_z) 
column_split=factor(as.character(clindata$Group),levels=c("Flagellin","Untreated")) 
Heatmap(my_data,cluster_rows=FALSE,cluster_columns=FALSE,name="Zscore",column_split=column_split,show_column_names=FALSE,row_names_side="left")

