library(Seurat, lib.loc="/Users/jacksonna/R-dev")
library(plotrix)
library(scales)
library(ggplot2)
library(cowplot)
library(openxlsx)
source("CommonFunctions.r")



##################################################### 1. SEURAT analysis using the full filtered dataset
#####Read in processed data file (where gene expressed in fewer than 0.1% of cell were removed)
##### Available here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134174
T15f<-read.table("Processed_invivo_raw.txt",sep="\t",header=T)
#Read in metadata
T15_metadata<-read.table("../Metadata.txt",sep="\t",header=T,stringsAsFactors=F)
T15_metadata<-data.frame("donor"=T15_metadata$Donor,"Smoke_status"=T15_metadata$Smoke_status,
	"cluster_ident"=T15_metadata$cluster_ident,"subcluster_ident"=T15_metadata$subcluster_ident,
	stringsAsFactors=F,row.names=T15_metadata$Cell)

#Create Seurat object and carry out integration
#bsub -e test.err -o test.out -R "rusage[mem=64000]" "module load R/3.5.1 && env R_MAX_VSIZE=700Gb R CMD BATCH reintegrateFull_scTransform_3000Features_run2_0.3propMT.r"
min.cells<-floor(0.001 * dim(T15f)[2])
T15 <- CreateSeuratObject(counts = T15f, min.cells = min.cells, min.features = 0, 
	meta.data = T15_metadata, project = "T15")




###########################Integration

#To construct a reference, we will identify ‘anchors’ between the individual datasets. First, we split the combined object into a list, with each dataset as an element.
T15.list<-SplitObject(T15,split.by="donor")

#Do Seurat's scTransform plug in in liu of the normal normalization and feature selection (~20 minutes on cluster)
for (i in 1:length(T15.list)) {
    T15.list[[i]] <- SCTransform(T15.list[[i]], verbose = FALSE)
}

#Next, select features for downstream integration, and run PrepSCTIntegration, which ensures that 
#all necessary Pearson residuals have been calculated
T15.features <- SelectIntegrationFeatures(object.list = T15.list, nfeatures = 3000)
T15.list <- PrepSCTIntegration(object.list = T15.list, anchor.features = T15.features, verbose = TRUE)

#Now find integration "anchors" (uses CCA) based on 30 dimensions (~3 hours on cluster)
T15.anchors <- FindIntegrationAnchors(T15.list, anchor.features = T15.features, dims = 1:30, normalization.method = "SCT")

#Now integrate the datasets (~2 hours on cluster)
T15_int <- IntegrateData(anchorset = T15.anchors, dims = 1:30, normalization.method = "SCT")

#Switch to integrated assay. The variable features of this assay are automatically set during IntegrateData
DefaultAssay(T15_int) <- "integrated"

#Run pca
T15_int<- RunPCA(T15_int, npcs = 30, verbose = TRUE)

#Elbow
ElbowPlot(T15_int,ndims=30) #How many PCs

#Look at loadings for each dim
VizDimLoadings(object = T15_int, dims = 1:2)

#Look at PCs
DimPlot(T15_int, reduction="pca")

#Look at heat maps
pdf("DimHeatmaps.pdf")
DimHeatmap(T15_int, reduction="pca", dims = 1:6, cells = 500, balanced = TRUE)
DimHeatmap(T15_int, reduction="pca", dims = 7:12, cells = 500, balanced = TRUE)
DimHeatmap(T15_int, reduction="pca", dims = 13:18, cells = 500, balanced = TRUE)
DimHeatmap(T15_int, reduction="pca", dims = 19:24, cells = 500, balanced = TRUE)
DimHeatmap(T15_int, reduction="pca", dims = 25:30, cells = 500, balanced = TRUE)
dev.off()


#Don't forget to normalize the RNA data in case I need it
T15_int<-NormalizeData(T15_int, verbose = FALSE,assay="RNA")









################################################################# 2. Clustering and UMAP
#Input
nPCs<-30
colorVec<-c("blue","royalblue4","mediumpurple3","paleturquoise3","paleturquoise2","lightsteelblue1","lightskyblue","purple",
	"slategrey","darkseagreen","turquoise4","peachpuff4","brown","brown2","yellow","deeppink","deeppink4","lightcoral","pink","yellowgreen",
	"rosybrown","black","darkgoldenrod","chocolate","gold","greenyellow","yellow4","yellow3","wheat","azure3","sienna4",
	"darkolivegreen","darkolivegreen1","forestgreen","darkslategray","springgreen")
	
#Run umap across different values of hyperparameters
#n_neighbors - In general this parameter should often be in the range 5 to 50 (lower values optimize for local structure; default is 30)
#min_dist - Sensible values are in the range 0.001 to 0.5 (lower values optimize for local structure; default is 0.3).
T15_int <- RunUMAP(T15_int, reduction.use = "pca", dims=1:nPCs, n.neighbors = 10, min.dist = 0.35)


#Now do clustering
T15_int <- FindNeighbors(object = T15_int, dims = 1:nPCs, k.param = 20)

#To break up c1a, c1b, and c1c do res = 0.3
#To break out the old c9 cluster, look at cluster 34 at res=2
T15_int <- FindClusters(T15_int, reduction.type="pca", resolution=0.22, algorithm=1)


#pdf("T15_UMAP_clusters.pdf",height=5,width=7)
dev.new(height=5,width=7)
DimPlot(T15_int, reduction='umap', label=T, pt.size = 0.05, cols=colorVec)

dev.off()


#Plot by donor, nGene, etc
p1<-DimPlot(T15_int, reduction = "umap", group.by = "donor",pt.size=0.2, cols=c("midnightblue","yellow","orange","pink","tan","darkgreen",
	"violetred4","slateblue4","blue","black","orangered4","yellowgreen","greenyellow","wheat","darkslategray")) + NoAxes()
p2<-DimPlot(T15_int, reduction = "umap", group.by = "Smoke_status",pt.size=0.2, cols = c("red","orange","black")) + NoAxes()
p3<-FeaturePlot(T15_int, "nFeature_RNA",cols=c("gray", "blue"),pt.size=0.05,min.cutoff="q5",max.cutoff="q95") + NoAxes()
p4<-FeaturePlot(T15_int, "nCount_RNA",cols=c("gray", "blue"),pt.size=0.05,min.cutoff="q5",max.cutoff="q95") + NoAxes()
pdf("T15_UMAP_donor.pdf",height=10,width=14)
plot_grid(p1)
dev.off()
pdf("T15_UMAP_smoke.pdf",height=10,width=14)
plot_grid(p2)
dev.off()
pdf("T15_UMAP_nGene.pdf",height=5,width=8)
plot_grid(p3)
dev.off()
pdf("T15_UMAP_nUMI.pdf",height=5,width=8)
plot_grid(p4)
dev.off()




#Rename clusters to be previous names
treat.col<-data.frame("clusters_10"=T15_int@meta.data$integrated_snn_res.0.22,row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col[,1]<-as.character(treat.col[,1])
treat.col[rownames(T15_int@meta.data)[which(T15_int@meta.data$integrated_snn_res.0.22 == "2")],]<-"c1"
treat.col[rownames(T15_int@meta.data)[which(T15_int@meta.data$integrated_snn_res.0.22 == "6")],]<-"c1c"
treat.col[rownames(T15_int@meta.data)[which(T15_int@meta.data$integrated_snn_res.0.22 == "1")],]<-"c2"
treat.col[rownames(T15_int@meta.data)[which(T15_int@meta.data$integrated_snn_res.0.22 == "5")],]<-"c3"
treat.col[rownames(T15_int@meta.data)[which(T15_int@meta.data$integrated_snn_res.0.22 == "0")],]<-"c4"
treat.col[rownames(T15_int@meta.data)[which(T15_int@meta.data$integrated_snn_res.0.22 == "3")],]<-"c5"
treat.col[rownames(T15_int@meta.data)[which(T15_int@meta.data$integrated_snn_res.0.22 == "8")],]<-"c6ab"
treat.col[rownames(T15_int@meta.data)[which(T15_int@meta.data$integrated_snn_res.0.22 == "9")],]<-"c6c"
treat.col[rownames(T15_int@meta.data)[which(T15_int@meta.data$integrated_snn_res.0.22 == "7")],]<-"c7"
treat.col[rownames(T15_int@meta.data)[which(T15_int@meta.data$integrated_snn_res.0.22 == "4")],]<-"c8"
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="clusters_10")

#Colors
colorVec<-c("#1619BA","mediumpurple3","#3281FF","darkturquoise","tomato","orange","gold4","tan","saddlebrown","green3")

#pdf("T15_UMAP_clusters_10.pdf",height=5,width=7)
dev.new(height=5,width=6)
DimPlot(T15_int, reduction='umap', label=T, pt.size = 0.08, cols=colorVec,group.by="clusters_10") + NoLegend()

#Repeat for just SMG
colorVec<-c("grey","grey","grey","darkturquoise","grey","grey","grey","grey","saddlebrown","grey")
dev.new(height=5,width=5.5)
DimPlot(T15_int, reduction='umap', label=F, pt.size = 0.08, cols=colorVec,group.by="clusters_10") + NoLegend()





#UMAP with smoker IDs
#Colors
colorVec<-c("grey","red","grey","black")

#pdf("T15_UMAP_smokers_FIGURE.pdf",height=5,width=7)
dev.new(height=5,width=5.3)
DimPlot(T15_int, reduction='umap', label=F, pt.size = 0.01, cols=colorVec,group.by="Smoke_status") + NoLegend()

png("T15_UMAP_smokers_FIGURE.png",height=5,width=5.3,units="in",res=400)
DimPlot(T15_int, reduction='umap', label=F, pt.size = 0.01, cols=colorVec,group.by="Smoke_status") + NoLegend()
dev.off()

#Do for donors as well
#Get new ordering of donors, with first nonsmokers, then smokers, then other
T15_int@meta.data$donor<-factor(T15_int@meta.data$donor,levels=c("T126","T137","T153","T164","T165","T166",
	"T85","T90","T101","T120","T154","T167","T84","T89","T121"))
colorVec<-c("pink","tan","greenyellow","yellow","lightblue","grey","purple4","slateblue","midnightblue",
	"black","orangered4","darkgreen","deeppink","darkorange1","deepskyblue")

#pdf("T15_UMAP_donor_FIGURE.pdf",height=5,width=7)
dev.new(height=5,width=5.3)
DimPlot(T15_int, reduction='umap', label=F, pt.size = 0.01, cols=colorVec,group.by="donor") + NoLegend()

png("T15_UMAP_donor_FIGURE.png",height=5,width=5.3,units="in",res=400)
DimPlot(T15_int, reduction='umap', label=F, pt.size = 0.01, cols=colorVec,group.by="donor") + NoLegend()
dev.off()



