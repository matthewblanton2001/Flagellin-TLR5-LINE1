#Load Libraries
library(tidyverse)
library(dplyr)
library(purrr)
library(broom)
library(ggplot2)
library(ggpubr)

#Load Outputs from QIIME2 and PICRUSt2
flagellar_kos<-c("K02406")
contrib<-read_tsv("pred_metagenome_contrib.tsv.gz")
flagellar<-contrib%>%filter(`function`%in%flagellar_kos)
tax<-read_tsv("taxonomy.tsv")
flagellar_tax<-flagellar%>%left_join(tax,by=c("taxon"="Feature ID"))
flagellar_tax <- flagellar_tax %>%
  mutate(
    Genus = str_extract(Taxon, "g__[^;]+"),
    Genus = str_remove(Genus, "g__"),
    Genus = str_trim(Genus)
  )

#Analyze the relative fliC abundance from selected genera
genus_fliC<-flagellar_tax%>%group_by(sample,Genus)%>%summarise(fliC_abundance=sum(taxon_function_abun),.groups="drop")%>%filter(!is.na(Genus),Genus !="",Genus !="uncultured",Genus !="metagenome",Genus !="Ambiguous_taxa",fliC_abundance>0)
top_genera<-genus_fliC%>%group_by(Genus)%>%summarise(total=sum(fliC_abundance))%>%arrange(desc(total))
meta<-read.csv("metadata_dataset1345.csv")
genus_fliC_meta<-genus_fliC%>%left_join(meta,by=c("sample"="sample.id"))%>%filter(!is.na(Disease))
genera_of_interest<-c("Escherichia-Shigella","Salmonella","Selenomonas","Pseudomonas","Campylobacter","Delftia","Stenotrophomonas","Sphingomonas")
relative_fliC<-genus_fliC_meta%>%group_by(sample)%>%mutate(percent_fliC=100*fliC_abundance/sum(fliC_abundance))

#Plot
plot_df <- relative_fliC %>%filter(Genus %in% genera_of_interest) %>%group_by(sample, Disease) %>%summarise(percent_fliC = sum(percent_fliC),.groups = "drop")
ggplot(plot_df,aes(x=Disease,y=percent_fliC,color=Disease))+geom_boxplot(fill=NA,outlier.shape=NA,linewidth=1)+geom_jitter(width=0.15,size=2,alpha=0.8)+stat_compare_means(method="wilcox.test",label="p.format",label.x=1.5,label.y=max(plot_df$percent_fliC)*1.05,size=6)+scale_color_manual(values=c("Control"="deepskyblue","NSCLC"="orange"))+labs(x=NULL,y="Relative contribution to predicted fliC abundance (%)")+theme_bw()+theme(legend.position="none",axis.text=element_text(size=12),axis.title=element_text(size=14),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
plot_df %>%group_by(Disease) %>%summarise(n = n_distinct(sample))
