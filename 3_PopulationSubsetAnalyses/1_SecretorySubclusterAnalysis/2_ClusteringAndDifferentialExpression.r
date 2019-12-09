######################################################################### SUBCLUSTERING
library(Seurat, lib.loc="/Users/jacksonna/R-dev")
library(plotrix)
library(scales)
library(ggplot2)
library(cowplot)
library(openxlsx)
source("CommonFunctions.r")




#################################### Luminal and Secretory cell clustering (after re-integrating)

#Bring in the integrated dataset (T15_lumSec)
load("subcluster_c4_5_reintegrate_3000.rda")

#Pick num PCs
ElbowPlot(T15_lumSec,ndims=30) #How many PCs

#Look at heat maps
pdf("DimHeatmaps.pdf")
DimHeatmap(T15_lumSec, reduction="pca", dims = 1:6, cells = 500, balanced = TRUE)
DimHeatmap(T15_lumSec, reduction="pca", dims = 7:12, cells = 500, balanced = TRUE)
DimHeatmap(T15_lumSec, reduction="pca", dims = 13:18, cells = 500, balanced = TRUE)
DimHeatmap(T15_lumSec, reduction="pca", dims = 19:24, cells = 500, balanced = TRUE)
DimHeatmap(T15_lumSec, reduction="pca", dims = 25:30, cells = 500, balanced = TRUE)
dev.off()




#########Clustering
#Input
nPCs<-30
colorVec<-c("darkorange1","black","yellowgreen","green","darkolivegreen","blue","turquoise","grey","midnightblue","yellow","red","purple")

#Best parameters for the 3000 gene integration
#min_dist - Sensible values are in the range 0.001 to 0.5 (lower values optimize for local structure; default is 0.3).
T15_lumSec <- RunUMAP(T15_lumSec, reduction.use = "pca", dims=1:nPCs, n.neighbors = 50, min.dist = 0.3)
T15_lumSec <- FindNeighbors(object = T15_lumSec, dims = 1:nPCs, k.param = 50)
T15_lumSec <- FindClusters(T15_lumSec, reduction.type="pca", resolution=0.54, algorithm=1)

#Plot preliminary UMAP
#pdf("T15_lumSec_UMAP_subclusters_3000.pdf",height=5,width=7)
dev.new(height=5,width=7)
DimPlot(T15_lumSec, reduction='umap', label=T, pt.size = 0.5, cols=colorVec)



#Look at overlay of the old T15 secretory and luminal clusters
#pdf("T15_lumSec_UMAP_subclusters_oldT15clusterColors.pdf")
dev.new(height=5,width=7)
DimPlot(T15_lumSec, reduction='umap', label=F, pt.size = 0.8, cols=c("tomato","darkgoldenrod","darkorange2","brown","black"),group.by="clusters_18")





#Now I've got final clusters, rename them according to standard naming
#Overlay original clusters
treat.col<-data.frame("subclusters_10"=rep(NA,ncol(GetAssayData(T15_lumSec))),row.names=rownames(T15_lumSec@meta.data),stringsAsFactors=F)
treat.col[which(T15_lumSec@meta.data$integrated_snn_res.0.54 == 8),]<-"c4a"
treat.col[which(T15_lumSec@meta.data$integrated_snn_res.0.54 == 5),]<-"c4b"
treat.col[which(T15_lumSec@meta.data$integrated_snn_res.0.54 == 2),]<-"c4c"
treat.col[which(T15_lumSec@meta.data$integrated_snn_res.0.54 == 3),]<-"c4d"
treat.col[which(T15_lumSec@meta.data$integrated_snn_res.0.54 == 1),]<-"c4e"
treat.col[which(T15_lumSec@meta.data$integrated_snn_res.0.54 == 0),]<-"c4f"
treat.col[which(T15_lumSec@meta.data$integrated_snn_res.0.54 == 4),]<-"c5a"
treat.col[which(T15_lumSec@meta.data$integrated_snn_res.0.54 == 6),]<-"c5b"
treat.col[which(T15_lumSec@meta.data$integrated_snn_res.0.54 == 7),]<-"c5c"
treat.col[which(T15_lumSec@meta.data$integrated_snn_res.0.54 == 9),]<-"c4g_5d"
T15_lumSec<-AddMetaData(T15_lumSec,metadata=treat.col,col.name="subclusters_10")

colorVec<-c("lightcoral","deeppink","deeppink4","rosybrown","gold","pink","black","darkgoldenrod","darkorange2","brown")

#Plot UMAP
#pdf("T15_lumSec_UMAP_subclusters.pdf")
dev.new(height=5,width=7)
DimPlot(T15_lumSec, reduction='umap', label=F, pt.size = 0.8, cols=colorVec, group.by="subclusters_10")

#png("T15_lumSec_UMAP_subclusters.pdf")
dev.new(height=5,width=7)
DimPlot(T15_lumSec, reduction='umap', label=F, pt.size = 0.8, cols=colorVec, group.by="subclusters_10")

save(list="T15_lumSec",file="T15_lumSec_3000.rda")




















############################################# Do differential expression among these clusters

#Do DE
Idents(T15_lumSec) <- T15_lumSec@meta.data$subclusters_10
clusterList<-sort(as.character(unique(Idents(T15_lumSec))))
markers.list <- list()
for(i in 1:length(clusterList)) {
    markers.list[[i]] <- FindMarkers(T15_lumSec, ident.1=clusterList[i], min.pct=0.1, logfc.threshold=0.25, only.pos=T, 
    	assay="SCT", test.use = "LR", latent.vars = "smoke")
}

#Combine tables
comparisonVec<-clusterList

#Get vector of clusters for calculating pct values from RNA
pctIndexList<-list()
for(i in 1:length(clusterList)){
	currClusterIndex<-which(clusterList == clusterList[i])
	pctIndexList[[length(pctIndexList) + 1]]<-list(clusterList[i],clusterList[-currClusterIndex])
}

#Get RNA counts and cluster assignments
counts<-as.matrix(GetAssayData(T15_lumSec,assay="RNA",slot="data"))
clusterAssignments<-T15_lumSec@meta.data$subclusters_10

#Run OrganizeDF_list
T15_1againstAll_DEGs<-OrganizeDF_list(datasetList=markers.list,comparisonVec=comparisonVec,pctIndexList=pctIndexList,
	counts=counts,clusterAssignments=clusterAssignments)

#Write all to single txt file
write.table(T15_1againstAll_DEGs,sep="\t",quote=FALSE,row.names=FALSE,file="DEG_Tables/T15_1againstAll_DEGs.txt")

#Write to excel
writeDEGsToExcel(DEG_table = T15_1againstAll_DEGs, savedFile = "DEG_Tables/T15_1againstAll_DEGs.xlsx")

#Export table for Enrichr
T15_1againstAll_DEGs_forEnrichr<-T15_1againstAll_DEGs[which(T15_1againstAll_DEGs$p_adj_FDR < 0.05),] #Only good DEGs

#Write all to single txt file
write.table(T15_1againstAll_DEGs_forEnrichr,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,
	file="Enrichr/T15_1againstAll_DEGs_forEnrichr.txt")

#Do enrichments
DEG_table<-read.table("Enrichr/T15_1againstAll_DEGs_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"T15_1againstAll_DEGs"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)









