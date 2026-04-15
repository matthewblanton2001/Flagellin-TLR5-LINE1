
# Libraries
# ===============================
library(SummarizedExperiment)
library(GEOquery)
library(tidyverse)
library(survival)
library(TCGAbiolinks)
library(survminer)
library(edgeR)
library(limma)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(lme4)
library(emmeans)
library(showtext)
library(sysfonts)

font_add("arial", "C:/Windows/Fonts/arial.ttf")
showtext_auto()
# ===============================
# 1. Download and prepare TCGA-BLCA
# ===============================
projects <- c("TCGA-LUAD")
query <- GDCquery(
  project = projects,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query, method = "api", files.per.chunk = 100)
tcga <- GDCprepare(query)
expr <- assay(tcga)
clin <- as.data.frame(colData(tcga))

# ===============================
# 2. Clinical processing
# ===============================
clin <- clin %>%
  mutate(
    OS_time = ifelse(is.na(days_to_death), days_to_last_follow_up, days_to_death),
    OS_event = ifelse(vital_status == "Dead", 1, 0),
    sample_id = rownames(.)
  )

# ===============================
# 3. Convert Ensembl IDs to gene symbols
# ===============================
ensembl_ids <- gsub("\\..*", "", rownames(expr))
symbols <- mapIds(org.Hs.eg.db,
                  keys = ensembl_ids,
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")
keep<-!is.na(symbols)
expr<-expr[keep, ]
symbols<-symbols[keep]
rownames(expr) <- symbols
expr<-expr%>%as.data.frame()%>%mutate(gene=rownames(.))%>%group_by(gene)%>%summarise(across(everything(),sum))%>%as.data.frame()
rownames(expr)<-expr$gene
expr$gene<-NULL
expr[is.na(expr)]<-0

# 4. Normalize TCGA counts
# ===============================
dge <- DGEList(counts = expr)
dge <- calcNormFactors(dge)
expr_norm <- cpm(dge, log = TRUE, prior.count = 1)

# 5. Align samples
common_samples<-intersect(colnames(expr_norm),clin$sample_id)
expr_norm<-expr_norm[,common_samples]
clin<-clin%>%filter(sample_id%in%common_samples)

# 6. Extract TLR5 Expression
"TLR5"%in%rownames(expr_norm)
tlr5_expr<-expr_norm["TLR5",]
clin$TLR5<-as.numeric(tlr5_expr[clin$sample_id])
summary(clin$TLR5)

# 7. Survial Time in Years
clin<-clin%>%mutate(OS_years=OS_time/365,OS_event=OS_event)

# 8. High/Low TLR5 expression

clin$TLR5_group<-case_when(clin$TLR5>=quantile(clin$TLR5,0.75,na.rm=TRUE)~"High expression",clin$TLR5<=quantile(clin$TLR5,0.25,na.rm=TRUE)~"Low expression",TRUE~NA_character_)
clin<-clin%>%filter(!is.na(TLR5_group))
clin$TLR5_group<-factor(clin$TLR5_group,levels=c("Low expression","High expression"))



# 9. Survival Analysis
surv_obj<-Surv(time=clin$OS_years,event=clin$OS_event)
fit<-survfit(surv_obj~TLR5_group,data=clin)
group_counts<-table(clin$TLR5_group)
n_low<-group_counts["Low expression"]
n_high<-group_counts["High expression"]
legend_labels<-c(paste0("Low expression (n=",n_low,")"),paste0("High expression (n=",n_high,")"))
ggsurvplot(fit,data=clin,pval=TRUE,risk.table=FALSE,conf.int=FALSE,censor=FALSE,ggtheme=theme_classic(),xlab="Time (Years)",ylab="OS Probability",legend="right",legend.labs=legend_labels,palette=c("deepskyblue","orange"),break.time.by=2,xlim=c(0,max(clin$OS_years,na.rm=TRUE)),surv.median.line="hv")



