######################################################################### SUBCLUSTERING

#Download the processed Seurat R object, "Processed_invitro_seurat.Rdata" from GEO,
#obtainable from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134174
#Then load the object it contains ("T15_int") into the session.
#T15_int is a Seurat dataset that contains the processed, integrated dataset
load("Processed_invitro_seurat.Rdata")


#################################### Ciliated Cells

################### Subset and then use the original integration, but re-select the top varient genes
#Get subset
Idents(T15_int) <- T15_int@meta.data$clusters_10
T15_cil<-subset(T15_int,idents=c("c8"))
T15_cil<-ScaleData(T15_cil, assay="integrated")
T15_cil<-FindVariableFeatures(T15_cil, selection.method = "vst", nfeatures = 1500,assay="integrated")


#Run pca
T15_cil<- RunPCA(T15_cil, npcs = 30, verbose = TRUE)
ElbowPlot(T15_cil,ndims=30) #How many PCs

#Look at heat maps
pdf("DimHeatmaps.pdf")
DimHeatmap(T15_cil, reduction="pca", dims = 1:6, cells = 500, balanced = TRUE)
DimHeatmap(T15_cil, reduction="pca", dims = 7:12, cells = 500, balanced = TRUE)
DimHeatmap(T15_cil, reduction="pca", dims = 13:18, cells = 500, balanced = TRUE)
DimHeatmap(T15_cil, reduction="pca", dims = 19:24, cells = 500, balanced = TRUE)
DimHeatmap(T15_cil, reduction="pca", dims = 25:30, cells = 500, balanced = TRUE)
dev.off()


################################Clustering
#Input
nPCs<-22
colorVec<-c("green3","darkolivegreen","greenyellow","red","purple","tan","black","green","magenta","saddlebrown")
	
#Run umap across different values of hyperparameters
#n_neighbors - In general this parameter should often be in the range 5 to 50 (lower values optimize for local structure; default is 30)
#min_dist - Sensible values are in the range 0.001 to 0.5 (lower values optimize for local structure; default is 0.3).
T15_cil <- RunUMAP(T15_cil, reduction.use = "pca", dims=1:nPCs, n.neighbors = 10, min.dist = 0.5)
T15_cil <- FindNeighbors(object = T15_cil, dims = 1:nPCs)
T15_cil <- FindClusters(T15_cil, reduction.type="pca", resolution=0.1, algorithm=1)

#pdf("T15_cil_UMAP_subclusters.pdf",height=5,width=7)
dev.new(height=5,width=7)
DimPlot(T15_cil, reduction='umap', label=T, pt.size = 2, cols=colorVec)
#FeaturePlot(T151, c("nCount_RNA"), cols=c("gray", "blue"), pt.size=1, min.cutoff="q5", max.cutoff="q95")
dev.off()


#Now rename subclusters we're using
#Overlay original clusters
treat.col<-data.frame("clusters_3"=rep(NA,ncol(GetAssayData(T15_cil))),row.names=rownames(T15_cil@meta.data),stringsAsFactors=F)
treat.col[which(T15_cil@meta.data$integrated_snn_res.0.1 == 0),]<-"c8c"
treat.col[which(T15_cil@meta.data$integrated_snn_res.0.1 == 1),]<-"c8b"
treat.col[which(T15_cil@meta.data$integrated_snn_res.0.1 == 2),]<-"c8a"
T15_cil<-AddMetaData(T15_cil,metadata=treat.col,col.name="clusters_3")

colorVec<-c("greenyellow","darkolivegreen","green3")

#pdf("T15_cil_UMAP_subclusters.pdf")
dev.new(height=5,width=7)
DimPlot(T15_cil, reduction='umap', label=F, pt.size = 2, cols=colorVec, group.by="clusters_3")












#####Overlay these new clusters onto the global UMAP
treat.col<-data.frame("cil_subclusters"=T15_int@meta.data$clusters_10,row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col[,1]<-as.character(treat.col[,1])
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_cil@meta.data)[which(T15_cil@meta.data$clusters_3 == "c8a")])],]<-"c8a"
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_cil@meta.data)[which(T15_cil@meta.data$clusters_3 == "c8b")])],]<-"c8b"
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_cil@meta.data)[which(T15_cil@meta.data$clusters_3 == "c8c")])],]<-"c8c"
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="cil_subclusters")
#Colors
colorVec<-c("#1619BA","mediumpurple3","#3281FF","darkturquoise","tomato","orange","gold4","tan","saddlebrown","greenyellow","darkolivegreen","green3")
#pdf("T15_UMAP_cil_3subclusters.pdf",height=5,width=7)
dev.new(height=5,width=7)
DimPlot(T15_int, reduction='umap', label=F, pt.size = 0.5, cols=colorVec,group.by="cil_subclusters")


















################################## Do differential expression #########################################

#Do DE
Idents(T15_cil) <- T15_cil@meta.data$clusters_3
clusterList<-sort(as.character(unique(Idents(T15_cil))))
markers.list <- list()
for(i in 1:length(clusterList)) {
    markers.list[[i]] <- FindMarkers(T15_cil, ident.1=clusterList[i], min.pct=0.1, logfc.threshold=0.25, only.pos=T, 
    	assay="SCT", test.use = "LR", latent.vars = "smoke")
}

#Secondly, do DE for each rare cell against all non-rare cells
#For this, need to bring in T15_int
Idents(T15_int) <- T15_int@meta.data$clusters_16a
markers.list[[length(markers.list) + 1]] <- FindMarkers(T15_int, ident.1="c8a", ident.2=
	unique(T15_int@meta.data$clusters_16a)[-grep("c8",unique(T15_int@meta.data$clusters_16a))],
	min.pct=0.1, logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "LR", latent.vars = "smoke")
markers.list[[length(markers.list) + 1]] <- FindMarkers(T15_int, ident.1="c8b", ident.2=
	unique(T15_int@meta.data$clusters_16a)[-grep("c8",unique(T15_int@meta.data$clusters_16a))],
	min.pct=0.1, logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "LR", latent.vars = "smoke")
markers.list[[length(markers.list) + 1]] <- FindMarkers(T15_int, ident.1="c8c", ident.2=
	unique(T15_int@meta.data$clusters_16a)[-grep("c8",unique(T15_int@meta.data$clusters_16a))],
	min.pct=0.1, logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "LR", latent.vars = "smoke")

#Combine tables
comparisonVec<-c("c8a_vs_cil","c8b_vs_cil","c8c_vs_cil","c8a_vs_noncil","c8b_vs_noncil","c8c_vs_noncil")

#Get vector of clusters for calculating pct values from RNA
pctIndexList<-list()
for(i in 1:length(clusterList)){
	currClusterIndex<-which(clusterList == clusterList[i])
	pctIndexList[[length(pctIndexList) + 1]]<-list(clusterList[i],clusterList[-currClusterIndex])
}
pctIndexList[[length(pctIndexList) + 1]]<-list("c8a",unique(T15_int@meta.data$clusters_16a)[-grep("c8",unique(T15_int@meta.data$clusters_16a))])
pctIndexList[[length(pctIndexList) + 1]]<-list("c8b",unique(T15_int@meta.data$clusters_16a)[-grep("c8",unique(T15_int@meta.data$clusters_16a))])
pctIndexList[[length(pctIndexList) + 1]]<-list("c8c",unique(T15_int@meta.data$clusters_16a)[-grep("c8",unique(T15_int@meta.data$clusters_16a))])


#Get RNA counts and cluster assignments
counts<-as.matrix(GetAssayData(T15_int,assay="RNA",slot="data"))
clusterAssignments<-T15_int@meta.data$clusters_16a

#Run OrganizeDF_list
T15_cil_DEGs<-OrganizeDF_list(datasetList=markers.list,comparisonVec=comparisonVec,pctIndexList=pctIndexList,
	counts=counts,clusterAssignments=clusterAssignments)

#Write all to single txt file
write.table(T15_cil_DEGs,sep="\t",quote=FALSE,row.names=FALSE,file="DEG_Tables/T15_cil_DEGs.txt")

#Write to excel
writeDEGsToExcel(DEG_table = T15_cil_DEGs, savedFile = "DEG_Tables/T15_cil_DEGs.xlsx")

#Export table for Enrichr
T15_cil_DEGs_forEnrichr<-T15_cil_DEGs[which(T15_cil_DEGs$p_adj_FDR < 0.05),] #Only good DEGs

#Write all to single txt file
write.table(T15_cil_DEGs_forEnrichr,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,
	file="Enrichr/T15_cil_DEGs_forEnrichr.txt")

#Do enrichments
DEG_table<-read.table("Enrichr/T15_cil_DEGs_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"T15_cil_DEGs"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)









