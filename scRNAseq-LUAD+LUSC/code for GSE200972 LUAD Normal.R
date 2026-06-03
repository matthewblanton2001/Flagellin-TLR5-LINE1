library(tidyverse)
library(Seurat)
expression_matrix<-ReadMtx(mtx="GSM6047628_P2_N_R_I_matrix.mtx.gz",features="GSM6047628_P2_N_R_I_features.tsv.gz",cells="GSM6047628_P2_N_R_I_barcodes.tsv.gz")
seurat_object<-CreateSeuratObject(expression_matrix)
expression_matrix<-ReadMtx(mtx="GSM6047630_P2_N_R_M_matrix.mtx.gz",features="GSM6047630_P2_N_R_M_features.tsv.gz",cells="GSM6047630_P2_N_R_M_barcodes.tsv.gz")
seurat_object2<-CreateSeuratObject(expression_matrix)
obj_list<-list(seurat_object,seurat_object2)
# Normalize + find features per object
obj_list <- lapply(obj_list, function(x) {
  x[["percent.mt"]]<-PercentageFeatureSet(x,pattern="^MT-")
  x<-subset(x,subset=nFeature_RNA>200&nFeature_RNA<6000&percent.mt<15)
  return(x)
})
obj_list <- lapply(obj_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)
  return(x)
})
# Find anchors
anchors<-FindIntegrationAnchors(object.list=obj_list,dims=1:30)
combined_seurat<-IntegrateData(anchorset=anchors,dims=1:30)
combined_seurat2<-combined_seurat
DefaultAssay(combined_seurat)<-"integrated"
combined_seurat<-ScaleData(combined_seurat)
combined_seurat<-RunPCA(combined_seurat)
ElbowPlot(combined_seurat,ndims=50)
combined_seurat<-FindNeighbors(combined_seurat,dims=1:15)
combined_seurat<-FindClusters(combined_seurat,resolution=0.5)
combined_seurat<-RunUMAP(combined_seurat,dims=1:15)
DefaultAssay(combined_seurat)<-"RNA"
combined_seurat<-JoinLayers(combined_seurat)
combined_markers<-FindAllMarkers(combined_seurat,only.pos=TRUE)
lung_markers <- list(
  # General epithelial
  Epithelial_cells = c("EPCAM","KRT8","KRT18","KRT19"),
  # Malignant / tumor epithelial
  Tumor_epithelial = c("EPCAM","KRT19","KRT17","SOX2","MSLN"),
  # LUSC markers
  LUSC_cells = c("TP63","KRT5","KRT14","SOX2","DSG3"),
  # LUAD markers
  LUAD_cells = c("SFTPA1","SFTPA2","SFTPB","NAPSA","AGER"),
  # Basal epithelial
  Basal_cells = c("KRT5","KRT14","TP63","KRT17"),
  # Club cells
  Club_cells = c("SCGB1A1","SCGB3A1","KRT19","CYP2F1"),
  # Alveolar type I
  AT1_cells = c("AGER","CAV1","EMP2","PDPN"),
  # Alveolar type II
  AT2_cells = c("SFTPA1","SFTPA2","SFTPB","SFTPC"),
  # Ciliated cells
  Ciliated_cells = c("FOXJ1","PIFO","TPPP3","DNAH5"),
  # Fibroblasts
  Fibroblasts = c("COL1A1","COL1A2","DCN","LUM"),
  # Cancer-associated fibroblasts
  CAFs = c("FAP","PDPN","ACTA2","TAGLN"),
  # Myofibroblasts
  Myofibroblasts = c("ACTA2","TAGLN","COL1A1","FN1"),
  # Endothelial
  Endothelial_cells = c("PECAM1","VWF","KDR","CDH5"),
  # Lymphatic endothelial
  Lymphatic_endothelial = c("PROX1","LYVE1","PDPN","FLT4"),
  # Monocytes
  CD14_Monocytes = c("CD14","LYZ","S100A8","S100A9"),
  CD16_Monocytes = c("FCGR3A","MS4A7","LST1"),
  # Macrophages
  Macrophages = c("CD68","APOE","C1QA","C1QB"),
  # Alveolar macrophages
  Alveolar_macrophages = c("FABP4","PPARG","MARCO","INHBA"),
  # Neutrophils
  Neutrophils = c("CXCL8","S100A8","S100A9","FCGR3B"),
  # Dendritic cells
  cDC1 = c("CLEC9A","BATF3","XCR1"),
  cDC2 = c("CD1C","FCER1A","CLEC10A"),
  pDCs = c("LILRA4","IRF7","GZMB"),
  # T cells
  CD4_Tcells = c("CD3D","CD4","IL7R","LTB"),
  CD8_Tcells = c("CD3D","CD8A","NKG7","GZMB"),
  # Exhausted T cells
  Exhausted_Tcells = c("PDCD1","LAG3","TIGIT","HAVCR2"),
  # Tregs
  Tregs = c("FOXP3","IL2RA","CTLA4","TIGIT"),
  # NK
  NK_cells = c("NKG7","PRF1","GNLY","KLRD1"),
  # B cells
  B_cells = c("CD79A","MS4A1","CD19","CD74"),
  # Plasma cells
  Plasma_cells = c("MZB1","JCHAIN","IGHG1","SDC1"),
  # Mast cells
  Mast_cells = c("TPSAB1","CPA3","KIT","MS4A2")
)
markers<-combined_markers
filtered_markers<-markers%>%filter(avg_log2FC>0.25,p_val_adj<0.05)
cluster_genes<-filtered_markers%>%group_by(cluster)%>%summarise(all_marker_genes=list(unique(gene)),.groups="drop")
annotations<-cluster_genes%>%rowwise()%>%mutate(scores=list(sapply(lung_markers,function(ref){sum(all_marker_genes%in%ref)/length(ref)})),Cell_Type=names(which.max(scores)),Score=max(unlist(scores)))%>%ungroup()
DimPlot(combined_seurat,reduction="umap",label=TRUE)
new_ids<-c("0"="CD8 T","1"="Macrophage","2"="CD4 T","3"="Macrophage","4"="Endothelial","5"="Macrophage","6"="Monocyte","7"="Neutrophil","8"="Epithelial","9"="Macrophage","10"="Fibroblast","11"="Epithelial","12"="Mast Cells","13"="Epithelial","14"="Epithelial","15"="Epithelial","16"="Endothelial","17"="Plasma Cells")
major_cols<-c("CD4 T"="red3","Epithelial"="forestgreen","CD8 T"="red","Fibroblast"="orange","Plasma Cells"="yellow2","Macrophage"="deepskyblue","Mast Cells"="purple","Endothelial"="magenta1","Neutrophil"="orangered1","Monocyte"="darkorange2")
combined_seurat<-RenameIdents(combined_seurat,new_ids)
combined_seurat$celltype<-Idents(combined_seurat)
library(ggrepel)
centers<-as.data.frame(Embeddings(combined_seurat,"umap"))%>%dplyr::mutate(cluster=Idents(combined_seurat))%>%dplyr::group_by(cluster)%>%dplyr::summarize(UMAP_1=mean(umap_1),UMAP_2=mean(umap_2))
p_umap<-DimPlot(combined_seurat,reduction="umap",cols=major_cols,pt.size=0.5)+labs(x="UMAP_1",y="UMAP_2")+theme_bw(base_size=14)+theme(panel.grid=element_blank(),plot.title=element_text(hjust=0.5))
library(patchwork)
tlr_genes<-c("TLR1","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7","TLR8","TLR9","TLR10")
present_tlrs<-tlr_genes[tlr_genes%in%rownames(combined_seurat)]
min_val<-0
max_val<-3
tlrfeatureplots<-lapply(present_tlrs,function(gene){FeaturePlot(combined_seurat,features=gene,reduction="umap",pt.size=0.3,order=TRUE)+labs(x="UMAP_1",y="UMAP_2")+scale_color_gradientn(colors=c("lightgrey","yellow","orange","red"),limits=c(min_val,max_val))+ggtitle(gene)+theme_bw(base_size=12)+theme(panel.grid=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_text(hjust=0.5))})
p_umap+geom_label_repel(data=centers,aes(x=UMAP_1,y=UMAP_2,label=cluster,fill=cluster),fontface="bold",size=4,color="black",box.padding=0.5,point.padding=0.5,segment.color="black",show.legend=FALSE)+scale_fill_manual(values=major_cols)+wrap_plots(tlrfeatureplots,ncol=3)
