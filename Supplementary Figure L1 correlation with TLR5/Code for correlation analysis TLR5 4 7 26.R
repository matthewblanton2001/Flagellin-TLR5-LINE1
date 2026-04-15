cnts<-read.delim("Flagellin_L1.cntTable",comment.char="#")
rownames(cnts)<-cnts$gene.TE
cnts<-cnts[,-1]
clean_ids<-sub("\\..*","",rownames(cnts))
cnts$clean_id<-clean_ids
library(tidyverse)
cnts_collapsed<-cnts%>%as.data.frame()%>%group_by(clean_id)%>%summarise(across(everything(),sum))
cnts_mat<-as.data.frame(cnts_collapsed)
rownames(cnts_mat)<-cnts_mat$clean_id
cnts_mat<-cnts_mat[,-1]
is_gene<-grepl("^ENSG",rownames(cnts_mat))
gene_counts<-cnts_mat[is_gene,]
te_counts<-cnts_mat[!is_gene,]
library(DESeq2)
coldata<-data.frame(Treatment=factor(c(rep("Flagellin",4),rep("Control",4))))
rownames(coldata)<-colnames(cnts_mat)
dds<-DESeqDataSetFromMatrix(countData=cnts_mat,colData=coldata,design=~Treatment)
vst_data<-vst(dds)
vst_mat<-assay(vst_data)


line_rows<-grep("^L1",rownames(vst_mat),value=TRUE)
line_mat<-vst_mat[line_rows,,drop=FALSE]
rownames(line_mat)<-sub(":.*","",rownames(line_mat))
line_keep<-c("L1HS","L1PA3","L1PA4","L1PA5","L1PA6","L1PA10")
line_mat<-line_mat[rownames(line_mat)%in%line_keep,]
library(rtracklayer)
gtf<-import("gencode.v38.annotation.gtf")
gtf_genes<-gtf[gtf$type=="gene"]
gtf_df<-data.frame(ensembl_gene_id=gtf_genes$gene_id,gene_symbol=gtf_genes$gene_name,stringsAsFactors=FALSE)
gtf_df$ensembl_gene_id<-sub("\\..*","",gtf_df$ensembl_gene_id)
ensembl_ids<-sub("\\..*","",rownames(vst_mat))
symbol_map<-setNames(gtf_df$gene_symbol,gtf_df$ensembl_gene_id)
gene_symbols<-gtf_df$gene_symbol[match(ensembl_ids,gtf_df$ensembl_gene_id)]
gene_symbols[is.na(gene_symbols)|gene_symbols==""]<-ensembl_ids[is.na(gene_symbols)|gene_symbols==""]
rownames(vst_mat)<-make.names(gene_symbols,unique=TRUE)
nrow(vst_mat)
length(gene_symbols)

genes <- c("TLR5","TP53","EGFR","ATM","SETDB1","MAPK1","KEAP1","KRAS")
genes_present<-genes[genes%in%rownames(vst_mat)]
if(length(genes_present)==0){stop("None of you genes found. Check mapping")}
gene_mat<-vst_mat[genes_present,,drop=FALSE]
flag_samples<-rownames(coldata)[coldata$Treatment=="Flagellin"]
gene_flag<-gene_mat[,flag_samples,drop=FALSE]
line_flag<-line_mat[,flag_samples,drop=FALSE]
cors_flag<-cor(t(gene_flag),t(line_flag),method="pearson",use="pairwise.complete.obs")
l1_order<-c("L1HS","L1PA3","L1PA4","L1PA5","L1PA6","L1PA10")
cors_flag<-t(cors_flag)
cors_flag<-cors_flag[l1_order,,drop=FALSE]
cors_flag<-cors_flag[,genes_present,drop=FALSE]
cors_flag<-t(cors_flag)
library(gridExtra)
breaks<-seq(-1,1,length.out=101)
colors<-colorRampPalette(c("blue","white","red"))(100)
library(pheatmap)
pheatmap(cors_flag,color=colors,breaks=breaks,cluster_rows=FALSE,cluster_cols=FALSE,fontsize_row=12,fontsize_col=14)





