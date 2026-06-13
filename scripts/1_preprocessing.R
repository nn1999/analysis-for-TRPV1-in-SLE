library(Matrix)
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)

x=list.files()
dir = c('~/Data/LYYlab2024/slejyj/matrix/H9/outs/filtered_feature_bc_matrix/','~/Data/LYYlab2024/slejyj/matrix/h10/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/h11/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/h12/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/h13/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/h14/filter_matrix/',
      '~/Data/LYYlab2024/slejyj_batch2/matrix/h17/filter_matrix/','~/Data/LYYlab2024/slejyj_batch2/matrix/h19/filter_matrix/','~/Data/LYYlab2024/slejyj_batch2/matrix/h20/filter_matrix/',
        '~/Data/LYYlab2024/slejyj/matrix/P8/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/P10/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/P12/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/P13/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/P14/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/P15/filter_matrix/',
        '~/Data/LYYlab2024/slejyj/matrix/P16/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/P17/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/P19/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/P20/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/P21/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/P22/filter_matrix/',
        '~/Data/LYYlab2024/slejyj/matrix/P23/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/P25/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/P26/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/P27/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/P28/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/P29/filter_matrix/',
        '~/Data/LYYlab2024/slejyj/matrix/P32/filter_matrix/','~/Data/LYYlab2024/slejyj/matrix/P33/filter_matrix/',
        '~/Data/LYYlab2024/slejyj_batch2/matrix/p34/filter_matrix/','~/Data/LYYlab2024/slejyj_batch2/matrix/p35/filter_matrix/','~/Data/LYYlab2024/slejyj_batch2/matrix/p36/filter_matrix/','~/Data/LYYlab2024/slejyj_batch2/matrix/p37/filter_matrix/','~/Data/LYYlab2024/slejyj_batch2/matrix/p38/filter_matrix/','~/Data/LYYlab2024/slejyj_batch2/matrix/p39/filter_matrix/')

names(dir) = c('h9',"h10","h11","h12","h13","h14","h17","h19","h20","P8","P10","P12","P13","P14","P15","P16","P17","P19","P20","P21","P22","P23","P25","P26","P27","P28","P29","P32","P33","p34","p35","p36","p37","p38","p39")


#####################1:create a list to read different samples data##########################
scRNAlist <- list()
for(i in (1:length(dir))){
  counts <- Read10X(data.dir = dir[i],gene.column=1)
  scRNAlist[[i]] <- CreateSeuratObject(counts, min.cells = 3, min.features =300)
  print(i)
}





######################2.qc######################
# for (i in 1:length(scRNAlist)) {
#   scRNAlist[[i]][["percent.mt"]]<- PercentageFeatureSet( scRNAlist[[i]], pattern = "^MT-")
#   col.num <- length(levels(scRNAlist[[i]]@active.ident))
#   violin <- VlnPlot(scRNAlist[[i]],
#                     features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
#                     cols =rainbow(col.num),
#                     pt.size = 0.01, #                     ncol = 3) +
#     theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
#   ggsave(paste("~/2023projects/LYY/sle.tcell.v2/results1/plotsqc/vlnplot_before_qc_",i,".jpg",sep = ""), plot =violin, width = 12, height = 6)
#   plot1=FeatureScatter(scRNAlist[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
#   plot2=FeatureScatter(scRNAlist[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#   pearplot <- CombinePlots(plots = list(plot1, plot2), nrow=1, legend="none")
#   ggsave(paste("~/2023projects/LYY/sle.tcell.v2/results1/plotsqc/pearplot_before_qc_",i,".jpg",sep = ""), plot =pearplot, width = 12, height = 6)
# }

scRNAlist.qc <- list()
for (i in 1:length(scRNAlist)) {
  print(i)
  scRNAlist[[i]][["percent.mt"]]<- PercentageFeatureSet( scRNAlist[[i]], pattern = "^MT-")
  
  scRNAlist.qc[[i]]<- subset(scRNAlist[[i]], subset = nFeature_RNA > 200& nFeature_RNA < 10000 & percent.mt < 10 & nCount_RNA < 100000)

  print(table(scRNAlist[[i]]@meta.data$orig.ident))
  print(table(scRNAlist.qc[[i]]@meta.data$orig.ident))
}








##################3.harmony###########################

scRNA1=scRNAlist.qc[[1]]
for (i in c(2:35)) {
  scRNA1 <- merge(scRNA1, scRNAlist.qc[[i]])}



#add meta data
table(scRNA1@meta.data$orig.ident)
sample_info = as.data.frame(colnames(scRNA1))
colnames(sample_info) = c('ID')
rownames(sample_info) = sample_info$ID
sample_info$sample = scRNA1@meta.data$orig.ident

library(stringr)
sample_info$group = 0


sample_info[which(sample_info$sample %in%names(dir)[1:12]),]$group='HC'
sample_info[which(sample_info$sample %in%names(dir)[13:35]),]$group='SLE'
scRNA1= AddMetaData(object = scRNA1, metadata = sample_info)

library(harmony)
DefaultAssay(scRNA1) <- "RNA"
scRNA1[['percent.mito']] <- PercentageFeatureSet(scRNA1, pattern = "^MT-")
scRNA1 <- NormalizeData(scRNA1,, normalization.method = "LogNormalize", scale.factor = 1e4) %>% FindVariableFeatures(selection.method = "vst",nfeatures = 2000) %>% ScaleData(verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito")) %>% RunPCA(verbose=FALSE)
system.time({scRNA1 <- RunHarmony(scRNA1, group.by.vars = "sample")})










#########3.3 PCA and UMAP############


# 
# scRNA1 <- ScaleData(scRNA1, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
# scRNA1 <- RunPCA(scRNA1,verbose = FALSE,npcs = 50)
# scRNA1 <- ProjectDim(object = scRNA1)
#ElbowPlot(object = scRNA1,ndims = 50)




scRNA1 <- FindNeighbors(scRNA1, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.5)
scRNA1 <- RunUMAP(scRNA1, reduction = "harmony", dims = 1:20)
DimPlot(scRNA1, reduction = "umap", group.by = "sample")
DimPlot(scRNA1, reduction = "umap", group.by = "group")

DimPlot(scRNA1, reduction = "umap", label = TRUE)
DimPlot(scRNA1.1, reduction = "umap", split.by="group")

#################################4.0 remove doublets#################################################################


library(DoubletFinder)

sce.all.list <- SplitObject(scRNA1, split.by = "orig.ident")
phe_lt <- lapply(names(dir), function(x){
 
  sce.all.filt=sce.all.list[[x]]
  sce.all.filt = FindVariableFeatures(sce.all.filt)
  sce.all.filt = ScaleData(sce.all.filt, 
                           vars.to.regress = c("nFeature_RNA", "percent_mito"))
  sce.all.filt = RunPCA(sce.all.filt, npcs = 20)
  sce.all.filt = RunTSNE(sce.all.filt, npcs = 20)
  sce.all.filt = RunUMAP(sce.all.filt, dims = 1:10)
  

  sweep.res.list <- paramSweep_v3(sce.all.filt, PCs = 1:10, sct = F)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
  bcmvn <- find.pK(sweep.stats)  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()   ## 
 DoubletRate = 0.076     
  homotypic.prop <- modelHomotypic(sce.all.filt$seurat_clusters) 
   nExp_poi <- round(DoubletRate*ncol(sce.all.filt)) 
   nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  

  
  sce.all.filt <- doubletFinder_v3( sce.all.filt, PCs = 1:10, pN = 0.25, pK = pK_bcmvn, 
                                    nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
  
  
  # name of the DF prediction can change, so extract the correct column name.
  DF.name = colnames(sce.all.filt@meta.data)[grepl("DF.classification", 
                                                   colnames(sce.all.filt@meta.data))]
  p5.dimplot=cowplot::plot_grid(ncol = 2, DimPlot(sce.all.filt, group.by = "orig.ident") + NoAxes(), 
                                DimPlot(sce.all.filt, group.by = DF.name) + NoAxes())
 
  ggsave(filename=paste0("~/2023projects/LYY/sle.tcell.v2/results2/plotsqc/doublet/doublet_dimplot_",x,".png"),
         plot=p5.dimplot)
  
  p5.vlnplot=VlnPlot(sce.all.filt, features = "nFeature_RNA", 
                     group.by = DF.name, pt.size = 0.1)
  
  ggsave(paste0("~/2023projects/LYY/sle.tcell.v2/results2/plotsqc/doublet/doublet_vlnplot_",x,".png"),
         plot=p5.vlnplot)
  print(table(sce.all.filt@meta.data[, DF.name] ))

  phe=sce.all.filt@meta.data
  phe
  
})

kpCells=unlist(lapply(phe_lt, function(x){
  table(x[,ncol(x)])
  rownames(x[ x[,ncol(x)]=='Singlet', ])    
}))


kp = colnames(scRNA1) %in% kpCells
scRNA1@meta.data$doublet =ifelse(colnames(scRNA1) %in% kpCells,  "Singlet","Doublet")
table(kp)
scRNA1.1=scRNA1[,kp]




DimPlot(scRNA1.1, reduction = "umap", label = TRUE)
scRNA1.1= RunTSNE(scRNA1.1, dims = 20)
DimPlot(scRNA1.1, reduction = "tsne")
DimPlot(scRNA1.1, reduction = "tsne", label = TRUE)
DimPlot(scRNA1.1, reduction = "umap", label = TRUE,group.by = "RNA_snn_res.0.5")





################4.annotation


markergenelist=read.table("~/2023projects/LYY/sle.tcell/markergeneplot_1level_pbmc_collection.csv",sep = ",",header = T)

a=DotPlot(scRNA1.1,features =markergenelist$gene,assay = 'RNA')+ RotatedAxis()+scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) 
TEMP=as.data.frame(a$data)
FeaturePlot(scRNA1.1, features = c("CD40LG", "CD4"),raster=FALSE,
            cols = c("lightgrey", "blue", "blue"),
            blend=T,blend.threshold=0)


Idents(scRNA1.1)="seurat_clusters"
test=FindMarkers(scRNA1.1,ident.1 = "16",test.use = "MAST",assay = "RNA")
#HEXIM1 




celltype=read.table("~/2023projects/LYY/sle.tcell.v2/results2/celltype_level1.csv",header = T,sep = ",")
scRNA1.1@meta.data$celltypelevel1 ="NA"
for(i in 1:nrow(celltype)){
  scRNA1.1@meta.data[which(scRNA1.1@meta.data$RNA_snn_res.0.5== celltype$ClusterID[i]),'celltypelevel1'] <- celltype$celltypelevel1[i]
}
library(tidydr)
DimPlot(scRNA1.1, group.by="celltypelevel1", label=F, label.size=3)+theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
 
#############FORMAL PLOTS#####################################
markergenelistgene=c("SLC4A10","SIGLEC6", "CST3","LYZ","CD79A", "CD79B" ,"CD3D","CD3E","CD4","CD40LG","CD8A" ,"CD8B")
#Idents(scRNA1.1)="RNA_snn_res.0.5"
Idents(scRNA1.1)="celltypelevel1"
scRNA1.1=subset(scRNA1.1,idents="NK",invert=T)
DotPlot(scRNA1.1,features =markergenelistgene,assay = "RNA")+ RotatedAxis()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45))+ 
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )
mycol=c("#9B5B33","#F5CFE4","#B383B9","#8FA4AE","#FCED82","#F5D2A8","#BBDD78","#393b79")

library(reshape2)
Idents(scRNA1.1)="sample"
scRNA1.1=subset(scRNA1.1,sample %in% c("P16","h19","P13"),invert=T) 
pB2_df <- table(scRNA1.1@meta.data$group,scRNA1.1@meta.data$celltypelevel1) %>% melt()
colnames(pB2_df) <- c("group","Cluster","Number")


pB2_df$group <- factor(pB2_df$group)
pB2_df$Cluster <- factor(pB2_df$Cluster)
pB2_df=filter(pB2_df,is.na(pB2_df$Cluster)==F)



library(ggplot2)
mycol=c("#FCED82","#F5CFE4","#B383B9","#F5D2A8","#8FA4AE","#9B5B33","#BBDD78","#393b79")
pB2_df$Cluster <- factor(pB2_df$Cluster, levels = c("CD8 T","Monocyte", "CD4 T" ,"DC","B","DN T"))  
pB4 <- ggplot(data = pB2_df, aes(x =group, y = Number, fill =Cluster )) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=mycol ) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  #ggsci::scale_color_d3("category20")+
  
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )
pB4


Idents(scRNA1.1)="celltypelevel1"
DimPlot(scRNA1.1, reduction = "umap", label =F, pt.size = 0.2,raster=F,cols =mycol ) + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.1, "inches"),angle=20,type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 1.5,hjust = 0.03)) 
#4.0x3.5
DimPlot(scRNA1.1, reduction = "umap", label =F, pt.size = 0.2,raster=F,cols =mycol,split.by = "group" ) + 
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.1, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 1.5,hjust = 0.03)) 
#8.0x3.5
save(scRNA1.1,file = "~/Data/LYYlab2024/slejyj_batch2/Rdatalog/scRNA1.1.35.Rdata")
