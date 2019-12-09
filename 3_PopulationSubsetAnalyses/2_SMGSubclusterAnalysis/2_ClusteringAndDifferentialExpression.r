######################################################################### SUBCLUSTERING
library(Seurat, lib.loc="/Users/jacksonna/R-dev")
library(plotrix)
library(scales)
library(ggplot2)
library(cowplot)
library(openxlsx)
source("../../CommonFunctions.r")



#################################### SMG cells (after re-integrating)

#Bring in the integrated dataset (T15_smg)
load("../../Subclustering_reintegrated_rda/T15_smg_noC1c_allGenes.rda")

#Pick num PCs
ElbowPlot(T15_smg,ndims=30) #How many PCs

#Look at heat maps
pdf("DimHeatmaps.pdf")
DimHeatmap(T15_smg, reduction="pca", dims = 1:6, cells = 500, balanced = TRUE)
DimHeatmap(T15_smg, reduction="pca", dims = 7:12, cells = 500, balanced = TRUE)
DimHeatmap(T15_smg, reduction="pca", dims = 13:18, cells = 500, balanced = TRUE)
DimHeatmap(T15_smg, reduction="pca", dims = 19:24, cells = 500, balanced = TRUE)
DimHeatmap(T15_smg, reduction="pca", dims = 25:30, cells = 500, balanced = TRUE)
dev.off()


################################Clustering
#Input
nPCs<-15
colorVec<-c("darkorange1","black","grey","pink","darkolivegreen","blue","turquoise","grey","midnightblue","yellow","red","purple","tan","pink")

#min_dist - Sensible values are in the range 0.001 to 0.5 (lower values optimize for local structure; default is 0.3).
T15_smg <- RunUMAP(T15_smg, reduction.use = "pca", dims=1:nPCs, n.neighbors = 50, min.dist = 0.3)
T15_smg <- FindNeighbors(object = T15_smg, dims = 1:nPCs, k.param = 30)
T15_smg <- FindClusters(T15_smg, reduction.type="pca", resolution=0.1, algorithm=1)

#pdf("T15_smg_UMAP_subclusters_1500.pdf",height=5,width=7)
dev.new(height=5,width=7)
DimPlot(T15_smg, reduction='umap', label=T, pt.size = 1.5, cols=colorVec)



#Look at the old T15 secretory and luminal clusters
#pdf("T15_smg_UMAP_subclusters_oldT15clusterColors.pdf")
dev.new(height=5,width=7)
DimPlot(T15_smg, reduction='umap', label=F, pt.size = 1.5, cols=c("slategrey","aquamarine","turquoise4","saddlebrown"),group.by="clusters_18")










###############Now rename subclusters we're using (5 subclusters version)
#Overlay original clusters
treat.col<-data.frame("subclusters_4"=rep(NA,ncol(GetAssayData(T15_smg))),row.names=rownames(T15_smg@meta.data),stringsAsFactors=F)
treat.col[which(T15_smg@meta.data$integrated_snn_res.0.1 == 0),]<-"c3a"
treat.col[which(T15_smg@meta.data$integrated_snn_res.0.1 == 1),]<-"c7"
treat.col[which(T15_smg@meta.data$integrated_snn_res.0.1 == 2),]<-"c3b"
treat.col[which(T15_smg@meta.data$integrated_snn_res.0.1 == 3),]<-"c3c"
T15_smg<-AddMetaData(T15_smg,metadata=treat.col,col.name="subclusters_4")

colorVec<-c("slategrey","aquamarine","turquoise4","saddlebrown")

#pdf("T15_smg_UMAP_5subclusters.pdf")
dev.new(height=5,width=7)
DimPlot(T15_smg, reduction='umap', label=F, pt.size = 1.5, cols=colorVec, group.by="subclusters_4")

save(list="T15_smg",file="T15_smg_noC1c_allGenes.rda")



#####Overlay these new clusters onto the global UMAP
treat.col<-data.frame("smg_subclusters"=T15_int@meta.data$clusters_10,row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col[,1]<-as.character(treat.col[,1])
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_smg@meta.data)[which(T15_smg@meta.data$subclusters_4 == "c3a")])],]<-"c3a"
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_smg@meta.data)[which(T15_smg@meta.data$subclusters_4 == "c7")])],]<-"c7"
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_smg@meta.data)[which(T15_smg@meta.data$subclusters_4 == "c3b")])],]<-"c3b"
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_smg@meta.data)[which(T15_smg@meta.data$subclusters_4 == "c3c")])],]<-"c3c"
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="smg_subclusters")
#Colors
colorVec<-c("#1619BA","mediumpurple3","#3281FF","slategrey","aquamarine","turquoise4","tomato","orange","gold4","tan","saddlebrown","green3")
#pdf("T15_UMAP_smg_5subclusters.pdf",height=5,width=7)
dev.new(height=5,width=7)
DimPlot(T15_int, reduction='umap', label=F, pt.size = 0.5, cols=colorVec,group.by="smg_subclusters")























############################################# Do differential expression among these clusters

#Do DE
Idents(T15_smg) <- T15_smg@meta.data$subclusters_4
clusterList<-sort(as.character(unique(Idents(T15_smg))))
markers.list <- list()
for(i in 1:length(clusterList)) {
    markers.list[[i]] <- FindMarkers(T15_smg, ident.1=clusterList[i], min.pct=0.1, logfc.threshold=0.25, only.pos=T, 
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
counts<-as.matrix(GetAssayData(T15_smg,assay="RNA",slot="data"))
clusterAssignments<-T15_smg@meta.data$subclusters_4

#Run OrganizeDF_list
T15_1againstAll_DEGs<-OrganizeDF_list(datasetList=markers.list,comparisonVec=comparisonVec,pctIndexList=pctIndexList,
	counts=counts,clusterAssignments=clusterAssignments)

#Write all to single txt file
write.table(T15_1againstAll_DEGs,sep="\t",quote=FALSE,row.names=FALSE,file="DEG_Tables/T15_1againstAll_DEGs_4subclusters.txt")

#Write to excel
writeDEGsToExcel(DEG_table = T15_1againstAll_DEGs, savedFile = "DEG_Tables/T15_1againstAll_DEGs_4subclusters.xlsx")

#Export table for Enrichr
T15_1againstAll_DEGs_forEnrichr<-T15_1againstAll_DEGs[which(T15_1againstAll_DEGs$p_adj_FDR < 0.05),] #Only good DEGs

#Write all to single txt file
write.table(T15_1againstAll_DEGs_forEnrichr,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,
	file="Enrichr/T15_1againstAll_DEGs_4subclusters_forEnrichr.txt")

#Do enrichments
DEG_table<-read.table("Enrichr/T15_1againstAll_DEGs_4subclusters_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"T15_1againstAll_DEGs_4subclusters"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)

















############################################# Repete DE excluding c7

#Do DE
Idents(T15_smg) <- T15_smg@meta.data$subclusters_4
clusterList<-sort(as.character(unique(Idents(T15_smg))))[-4]
markers.list <- list()
for(i in 1:length(clusterList)) {
    markers.list[[i]] <- FindMarkers(T15_smg, ident.1=clusterList[i], ident.2=clusterList[-i], min.pct=0.1, 
    	logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "LR", latent.vars = "smoke")
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
counts<-as.matrix(GetAssayData(T15_smg,assay="RNA",slot="data"))
clusterAssignments<-T15_smg@meta.data$subclusters_4

#Run OrganizeDF_list
T15_1againstAll_minusC7_DEGs<-OrganizeDF_list(datasetList=markers.list,comparisonVec=comparisonVec,pctIndexList=pctIndexList,
	counts=counts,clusterAssignments=clusterAssignments)

#Write all to single txt file
write.table(T15_1againstAll_minusC7_DEGs,sep="\t",quote=FALSE,row.names=FALSE,file="DEG_Tables/T15_1againstAll_minusC7_DEGs_4subclusters.txt")

#Write to excel
writeDEGsToExcel(DEG_table = T15_1againstAll_minusC7_DEGs, savedFile = "DEG_Tables/T15_1againstAll_minusC7_DEGs_4subclusters.xlsx")

#Export table for Enrichr
T15_1againstAll_minusC7_DEGs_forEnrichr<-T15_1againstAll_minusC7_DEGs[which(T15_1againstAll_minusC7_DEGs$p_adj_FDR < 0.05),] #Only good DEGs

#Write all to single txt file
write.table(T15_1againstAll_minusC7_DEGs_forEnrichr,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,
	file="Enrichr/T15_1againstAll_minusC7_DEGs_4subclusters_forEnrichr.txt")

#Do enrichments
DEG_table<-read.table("Enrichr/T15_1againstAll_minusC7_DEGs_4subclusters_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"T15_1againstAll_minusC7_DEGs_4subclusters"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)



#Muscle genes
unique(c("TPM2","TPM1","LMOD1","SORBS1","ACTG2","MYLK","ACTA2","MYL6","CALD1","ANXA6","MYH11","TLN1","MYL9","TPM2","TPM1",
"LMOD1","SORBS1","ACTG2","MYLK","ACTA2","GJA1","EDNRA","MYL6","SMTN","CALD1","SSPN","MYH11","DMD","VIM","TLN1","MYL9","CRYAB",
"PPP1R14A","PPP1R12A","EDN3","ITPR1","ACTG2","MYLK","ACTA2","EDNRA","MYL6","CALD1","KCNMB1","MYH11","PPP1R12B","MYL9","SGCE",
"MEF2C","TAGLN","SMTN","MRAS","FHL1","DMD","AEBP1","CNN1","CAV1","TPM1","PPP1R12B","MYL9"))

























##################################### Make heat map of DEGs among the three subclusters (including c7)


#Bring in c5M_c7 genes
T15_1againstAll_DEGs_forEnrichr<-read.table("Enrichr/T15_1againstAll_DEGs_4subclusters_forEnrichr.txt",sep="\t",header=T,stringsAsFactors=F)

c3c<-T15_1againstAll_DEGs_forEnrichr$gene[which(T15_1againstAll_DEGs_forEnrichr$comparison == "c3c")]
c3a<-T15_1againstAll_DEGs_forEnrichr$gene[which(T15_1againstAll_DEGs_forEnrichr$comparison == "c3a")]
c3b<-T15_1againstAll_DEGs_forEnrichr$gene[which(T15_1againstAll_DEGs_forEnrichr$comparison == "c3b")]
c7<-T15_1againstAll_DEGs_forEnrichr$gene[which(T15_1againstAll_DEGs_forEnrichr$comparison == "c7")]

#Randomize order of genes
c3c<-sample(c3c,length(c3c),replace=FALSE)
c3a<-sample(c3a,length(c3a),replace=FALSE)
c3b<-sample(c3b,length(c3b),replace=FALSE)
c7<-sample(c7,length(c7),replace=FALSE)

#Set colors
cols_clusters<-c("slategrey","aquamarine","turquoise4","saddlebrown")


#First, need to isolate those genes we want and cells we want in the order that we want them
#Randomize (to minimize donor differences)
cellOrdering<-c(
	sample(which(T15_smg@meta.data$subclusters_4 == "c3c"),
		length(which(T15_smg@meta.data$subclusters_4 == "c3c")),replace=FALSE),
	sample(which(T15_smg@meta.data$subclusters_4 == "c3a"),
		length(which(T15_smg@meta.data$subclusters_4 == "c3a")),replace=FALSE),
	sample(which(T15_smg@meta.data$subclusters_4 == "c3b"),
		length(which(T15_smg@meta.data$subclusters_4 == "c3b")),replace=FALSE),
	sample(which(T15_smg@meta.data$subclusters_4 == "c7"),
		length(which(T15_smg@meta.data$subclusters_4 == "c7")),replace=FALSE))

#Get gene ordering based on DEGs
geneOrdering<-rev(c(c3c,c3a,c3b,c7))

#Get module ordering (cluster for targeted DEGs)
moduleOrdering<-c(rep("turquoise4",length(c3c)),rep("slategrey",length(c3a)),rep("aquamarine",length(c3b)),rep("saddlebrown",length(c7)))
moduleOrdering<-rev(moduleOrdering)

#Now make a heatmap using heatmap3 of expression across genes and cells
mat<-as.matrix(GetAssayData(T15_smg,assay="SCT"))[geneOrdering,cellOrdering]

#Get side colors for population
#Note that the colors need to be in numerical order, not in the order of cellOrdering
clusters<-color.scale(as.numeric(as.factor(T15_smg@meta.data[cellOrdering,]$subclusters_4)),extremes=cols_clusters,color.spec="rgb")
ColSideColors<-data.frame(clusters=clusters)

#Get row colors for modules
modules<-color.scale(as.numeric(as.factor(moduleOrdering)),extremes=rev(c("turquoise4","slategrey","aquamarine","saddlebrown")), color.spec="rgb")
RowSideColors<-data.frame(modules=modules)

mycols = colorRampPalette(c("blue","white","red"))(1000)
breakscale = c(-8,seq(-1.2,1.2, length.out=length(mycols)-1), 8)

pdf('T15_smg_markersSCT_heatmap.pdf',height=24,width=18)
heatmap3(mat, scale="row",Rowv=NA,Colv=NA, 
	col=mycols,ColSideColors=as.matrix(ColSideColors),RowSideColors=as.matrix(RowSideColors),
	cexRow=0.1,cexCol=1,breaks=breakscale,labCol="",useRaster=T)
dev.off()














