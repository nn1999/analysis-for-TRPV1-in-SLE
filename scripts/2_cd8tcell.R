load("~/Data/LYYlab2024/slejyj_batch2/Rdatalog/scRNA1.1.35.Rdata")
Idents(scRNA1.1)="celltypelevel1"
scRNA1_cd8= subset(scRNA1.1,idents = "CD8 T")

DefaultAssay(scRNA1_cd8) <- "RNA"
Subset.list <- SplitObject(scRNA1_cd8, split.by = "sample")
for (i in 1:length(Subset.list)) {
  Subset.list[[i]] <- NormalizeData(Subset.list[[i]], verbose = FALSE)
  Subset.list[[i]] <- FindVariableFeatures(Subset.list[[i]], selection.method = "vst", nfeatures = 2000,verbose = FALSE)
}

#f your cells are grouping mostly by cell type and not by sample, batch correction is not necessary. If you find your cells are mostly clustering by sample, then batch correction is recommended.
scobj_cd8=Subset.list[[1]]
for (i in c(2:35)) {
 scobj_cd8 <- merge(scobj_cd8, Subset.list[[i]])} 


library(harmony)
DefaultAssay(scobj_cd8) <- "RNA"
scobj_cd8[['percent.mito']] <- PercentageFeatureSet(scobj_cd8, pattern = "^MT-")
scobj_cd8  <- NormalizeData(scobj_cd8, normalization.method = "LogNormalize", scale.factor = 1e4 ) %>% FindVariableFeatures(selection.method = "vst",nfeatures = 2000) %>% ScaleData(verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito")) %>% RunPCA(verbose=FALSE)
system.time({scobj_cd8  <- RunHarmony(scobj_cd8 , group.by.vars = "sample")})






#########3.3 PCA and UMAP############



ElbowPlot(object =scobj_cd8,ndims = 50)
scobj_cd8 <- FindNeighbors(scobj_cd8, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 1)
scobj_cd8 <- RunUMAP(scobj_cd8, reduction = "harmony", dims = 1:30)

DimPlot(scobj_cd8, reduction = "umap", label = TRUE)

FeaturePlot(scobj_cd8, features = c("CD40LG", "CD8A",'CD8B'),raster=FALSE,
            cols = c("lightgrey", "blue"),
            blend=F,blend.threshold=0)


###############annotation MNN################



markergenelist=read.table("~/2023projects/LYY/sle.tcell/markergene_cd8t_collection.csv",sep = ",",header = T)  

sum(duplicated(markergenelist$gene))
markergenelist=distinct(markergenelist,gene,keep_all = T)
Idents(scobj_cd8)="RNA_snn_res.1"
DotPlot(scobj_cd8,features =markergenelist$gene )+ RotatedAxis()
DimPlot(scobj_cd8, reduction = 'umap',label=T,raster=FALSE,pt.size = 0.3)



Idents(scobj_cd8)="RNA_snn_res.1"
tempmarkers=FindMarkers(object = scobj_cd8, ident.1 = "12",assay = 'RNA',only.pos = FALSE, test.use = 'MAST')
tempmarkers2=FindMarkers(object = scobj_cd8, ident.1 = "16",assay = 'RNA',only.pos = FALSE, test.use = 'MAST')

####################添加celltype信息##################


celltype=read.table("~/2023projects/LYY/sle.tcell.v2/results2/celltype_level3.csv",header = T,sep = ",")
#add a new column to metadata
scobj_cd8@meta.data$celltypelevel3="NA"

for(i in 1:nrow(celltype)){
  scobj_cd8@meta.data[which(scobj_cd8@meta.data$RNA_snn_res.1 == celltype$ClusterID[i]),'celltypelevel3'] <- celltype$celltypelevel2[i]
}




# ####1.markergeneplot##############



Idents(scobj_cd8)="celltypelevel3"
markergenelistgene=c("SELL","LEF1","TCF7","PRF1","GZMB","GZMK","EOMES","CD58","ITGAL","HNRNPLL","ADAM19")

DotPlot(scobj_cd8,features =markergenelistgene,assay = "RNA")+ RotatedAxis()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45))+ 
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )


  
  

library(tidydr)
mycol=c("#FCED82","#BBDD78","#EE934E")

Idents(scobj_cd8)=scobj_cd8$celltypelevel3
DimPlot(scobj_cd8, reduction = 'umap',label=F,raster=FALSE,pt.size = 0.3,cols = mycol,split.by = "group") +
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.1, "inches"),angle=20,type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 1.5,hjust = 0.03)) 


pB2_df <- table(scobj_cd8@meta.data$group,scobj_cd8@meta.data$celltypelevel3) %>% melt()
colnames(pB2_df) <- c("group","Cluster","Number")


pB2_df$group <- factor(pB2_df$group)
pB2_df$Cluster <- factor(pB2_df$Cluster)
pB2_df=filter(pB2_df,is.na(pB2_df$Cluster)==F)

mycol=c("#EE934E","#FCED82","#BBDD78")

library(ggplot2)
pB4 <- ggplot(data = pB2_df, aes(x =group, y = Number, fill =Cluster )) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=mycol ) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+


  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )
pB4




##########################corr plots############################

library(readxl)
gmt_fib<- readLines(paste0("~/2023projects/LYY/sle.tcell/","Cytotoxicity.gmt"),skipNul = T)
gmt_fib <- strsplit(gmt_fib, "\t")

names(gmt_fib) <- vapply(gmt_fib, function(y) y[1], character(1))
gmt_fib <- lapply(gmt_fib, "[", -c(1:2))

scobj_cd8<- AddModuleScore(scobj_cd8,
                           features = gmt_fib,
                           ctrl = 5,
                           name = "Cytotoxicity")     
allgene=c("TRPV1")
gene_cell_exp <- AverageExpression(scobj_cd8,
                                   features = allgene,
                                   group.by = 'sample',
                                   slot = 'data') 
average_cytotoxicity <- scobj_cd8@meta.data %>%
  group_by(sample) %>%               
  summarise(mean_cytotoxicity = mean(Cytotoxicity1, na.rm = TRUE)) 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
gene_cell_exp <-t(gene_cell_exp)
gene_cell_exp <-as.data.frame(gene_cell_exp)

gene_cell_exp <-rownames_to_column(gene_cell_exp,var="sample")
gene_cell_exp <-left_join(gene_cell_exp, average_cytotoxicity,by="sample")
gene_cell_exp$groupvec=c(rep("HC",9),rep("SLE",26))


data=gene_cell_exp

res <- rcorr( data$mean_cytotoxicity, data$V1)
p_value <- round(res$P[1,2],3)
cor_value <- round(res$r[1,2], 2)
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

data$mean_cytotoxicity=normalize(data$mean_cytotoxicity)

colnames(data)[2]="TRPV1"
ggplot(data,aes(mean_cytotoxicity,TRPV1))+
  geom_point(colour = "darkgrey", size = 3,position="identity",shape = 20)+
  #theme_bw()+   theme_classic()+  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.title.y = element_text(size = 12, face = "italic"),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size =12, color = "black"),  
    axis.text.y = element_text(size = 12, color = "black")  
  )+
  labs(title = paste0("R = ", cor_value, ", p = ", p_value), x = "Cytotoxicity score")+
  geom_smooth(method=lm, fill="#191516",colour="#0707CB",level=0.95)



#####MAFK#################

test=read.csv(file = "~/2023projects/LYY/sle.tcell.v2/results2/OUTPUT/t.deg.35samples.csv")
tf=read.csv(file = "~/2023projects/LYY/sle.tcell/humanTF.csv")
test=mutate(test,istf=ifelse(X %in% tf$Name,"yes","no"))
test1=filter(test,istf=="yes")

# test1=filter(test,test$X %in% genes$SYMBOL)
library(ggpubr)
library(ggthemes)
deg.data=test1
deg.data$p_val_adj[which(deg.data$p_val_adj<1.572166e-319)]=0
deg.data$logP<--log10(1e-323+deg.data$p_val_adj)
deg.data$Group ="not-significant"
deg.data$Group[which( (deg.data$p_val_adj <0.05)&(deg.data$avg_log2FC >0))]= "up-regulated"
deg.data$Group[which((deg.data$p_val_adj<0.05)&(deg.data$avg_log2FC <0))]="down-regulated"
table(deg.data$Group)

deg.data$Label="" 
colnames(deg.data)[1]="gene"

deg.top10.genes<-c("MAFK","TCF7","NFKB2","RUNX1","KLF2","CDC5L","MAF")
deg.data$Label[match(deg.top10.genes, deg.data$gene)]<- deg.top10.genes



ggscatter(deg.data,x="avg_log2FC",y="logP",color ="Group",
          palette =c("#2f5688","#BBBBBB","#CC0000"),size =1,label = deg.data$Label,font.label =12,repel =T,
          xlab = "log2FoldChange",ylab ="-log10(Adjust P-value)")+
  theme_base()+
  geom_hline(yintercept=1.30,linetype="dashed",color = "darkgrey",linewidth = 0.5)+
  geom_vline(xintercept=c(-0.5,0.5),linetype="dashed",color = "white",linewidth = 0)+
  theme(panel.border = element_blank())+
  ylim(0, 350) 





library(clusterProfiler)
library(org.Hs.eg.db)
slehcmarker=read.csv(file = "~/2023projects/LYY/sle.tcell.v2/results2/OUTPUT/cd8t.deg.csv")
colnames(slehcmarker)[1]="Symbol"
DEGs_up<- slehcmarker%>%filter(avg_log2FC>log2(1)& p_val_adj<0.05)

DEGs_down<- slehcmarker%>%filter(avg_log2FC<log2(1)& p_val_adj<0.05)

IDs_up<- bitr(DEGs_up$Symbol,
              fromType ='SYMBOL',
              toType =c('ENTREZID'),
              OrgDb ='org.Hs.eg.db')#
IDs_down<- bitr(DEGs_down$Symbol,
                fromType ='SYMBOL',
                toType = c('ENTREZID'),
                OrgDb ='org.Hs.eg.db')
enrichGo_result_up<- enrichGO(IDs_up$ENTREZID,
                              OrgDb ="org.Hs.eg.db",qvalueCutoff = 0.05,pvalueCutoff = 0.05,
                              ont="all",
                              readable = T)

dotplot(enrichGo_result_up,showCategory = 20,decreasing  = T,font.size = 8)

original_timeout=getOption("timeout")
options(timeout=180)
enrichGo_result_up_1 <- enrichKEGG(gene = IDs_up$ENTREZID, #需要分析的基因的EntrezID
                                   organism = "hsa",  #x
                                   pvalueCutoff =0.05, #设置pvalue界值
                                   qvalueCutoff = 0.05) #设置qvalue界值(FDR校正后的p值）
dotplot(enrichGo_result_up_1,showCategory = 20,decreasing  = T,font.size = 8)

enrichGo_result_down<- enrichGO(IDs_down$ENTREZID,
                                OrgDb ="org.Hs.eg.db",qvalueCutoff = 0.05,pvalueCutoff = 0.05,
                                ont="all",
                                readable = T)
dotplot(enrichGo_result_down,showCategory = 20,decreasing  = T,font.size = 8)
enrichGo_result_down_1 <- enrichKEGG(gene = IDs_down$ENTREZID, #需要分析的基因的EntrezID
                                     organism = "hsa",  #人
                                     pvalueCutoff =0.05, #设置pvalue界值
                                     qvalueCutoff = 0.05)
dotplot(enrichGo_result_down_1,showCategory = 20,decreasing  = T,font.size = 8)


library("enrichplot")
enrichGo_result_up_data<- data.frame(enrichGo_result_up)
enrichGo_result_down_data<- data.frame(enrichGo_result_down)
enrichGo_result_up_data_1<- data.frame(enrichGo_result_up_1)
enrichGo_result_down_data_1<- data.frame(enrichGo_result_down_1)

enrichGo_result_up_data <- mutate(enrichGo_result_up_data,up_down = "Up-regulated")
enrichGo_result_down_data<- enrichGo_result_down_data %>%mutate(up_down="Down-regulated", Count =-Count)

enrichGo_result_up_data_1 <- mutate(enrichGo_result_up_data_1,up_down = "Up-regulated")
enrichGo_result_down_data_1<- enrichGo_result_down_data_1 %>%mutate(up_down="Down-regulated", Count =-Count)

enrichGo_result <-rbind(enrichGo_result_up_data,enrichGo_result_down_data)
enrichGo_result_1 <-rbind(enrichGo_result_up_data_1,enrichGo_result_down_data_1)
write.csv(enrichGo_result,"~/2023projects/LYY/sle.tcell.v2/results2/OUTPUT/GO_SLEvsHC.CD8T.csv",row.names = F)
write.csv(enrichGo_result_1,"~/2023projects/LYY/sle.tcell.v2/results2/OUTPUT/KEGG_SLEvsHC.CD8T.csv",row.names = F)



