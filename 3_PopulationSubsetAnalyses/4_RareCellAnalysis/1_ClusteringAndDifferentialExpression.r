######################################################################### SUBCLUSTERING

#Download the processed Seurat R object, "Processed_invitro_seurat.Rdata" from GEO,
#obtainable from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134174
#Then load the object it contains ("T15_int") into the session.
#T15_int is a Seurat dataset that contains the processed, integrated dataset
load("Processed_invitro_seurat.Rdata")

#################################### Rare cells (tuft, ionocytes, and pnecs)

################### Subset and then use the original integration, but re-select the top varient genes
#Get subset
Idents(T15_int) <- T15_int@meta.data$clusters_10
T15_rare3<-subset(T15_int,idents=c("c6ab","c6c"))
T15_rare3<-ScaleData(T15_rare3,assay="integrated")
T15_rare3<-FindVariableFeatures(T15_rare3, selection.method = "vst", nfeatures = 1500,assay="integrated")

#Run pca
T15_rare3<- RunPCA(T15_rare3, npcs = 30, verbose = TRUE)
ElbowPlot(T15_rare3,ndims=30) #How many PCs

#Look at heat maps
pdf("DimHeatmaps.pdf")
DimHeatmap(T15_rare3, reduction="pca", dims = 1:6, cells = 500, balanced = TRUE)
DimHeatmap(T15_rare3, reduction="pca", dims = 7:12, cells = 500, balanced = TRUE)
DimHeatmap(T15_rare3, reduction="pca", dims = 13:18, cells = 500, balanced = TRUE)
DimHeatmap(T15_rare3, reduction="pca", dims = 19:24, cells = 500, balanced = TRUE)
DimHeatmap(T15_rare3, reduction="pca", dims = 25:30, cells = 500, balanced = TRUE)
dev.off()


################################Clustering
#Input
nPCs<-30
colorVec<-c("yellow3","wheat","yellowgreen")
	
#Run umap across different values of hyperparameters
#n_neighbors - In general this parameter should often be in the range 5 to 50 (lower values optimize for local structure; default is 30)
#min_dist - Sensible values are in the range 0.001 to 0.5 (lower values optimize for local structure; default is 0.3).
T15_rare3 <- RunUMAP(T15_rare3, reduction.use = "pca", dims=1:nPCs, n.neighbors = 10, min.dist = 0.3)
T15_rare3 <- FindNeighbors(object = T15_rare3, dims = 1:nPCs)
T15_rare3 <- FindClusters(T15_rare3, reduction.type="pca", resolution=0.1, algorithm=1)

#pdf("T15_rare3_UMAP_subclusters.pdf",height=5,width=7)
dev.new(height=5,width=7)
DimPlot(T15_rare3, reduction='umap', label=T, pt.size = 2, cols=colorVec)
#FeaturePlot(T151, c("nCount_RNA"), cols=c("gray", "blue"), pt.size=1, min.cutoff="q5", max.cutoff="q95")


#Now rename subclusters we're using
#Overlay original clusters
treat.col<-data.frame("clusters_3"=rep(NA,ncol(GetAssayData(T15_rare3))),row.names=rownames(T15_rare3@meta.data),stringsAsFactors=F)
treat.col[which(T15_rare3@meta.data$integrated_snn_res.0.1 == 0),]<-"c6a"
treat.col[which(T15_rare3@meta.data$integrated_snn_res.0.1 == 1),]<-"c6b"
treat.col[which(T15_rare3@meta.data$integrated_snn_res.0.1 == 2),]<-"c6c"
T15_rare3<-AddMetaData(T15_rare3,metadata=treat.col,col.name="clusters_3")

colorVec<-c("yellowgreen","yellow3","wheat")

#pdf("T15_rare3_UMAP_subclusters.pdf")
dev.new(height=5,width=7)
DimPlot(T15_rare3, reduction='umap', label=F, pt.size = 4, cols=colorVec, group.by="clusters_3")

FindMarkers(T15_rare3,ident.1 = c("11"), min.pct = 0.1, only.pos=T, logfc.threshold = 0.25, test.use = "wilcox",assay="SCT")


















################################## Do differential expression #########################################

#Do DE (each against the other rare cells)
Idents(T15_rare3) <- T15_rare3@meta.data$clusters_3
clusterList<-sort(as.character(unique(Idents(T15_rare3))))
markers.list <- list()
for(i in 1:length(clusterList)) {
    markers.list[[i]] <- FindMarkers(T15_rare3, ident.1=clusterList[i], min.pct=0.1, logfc.threshold=0.25, only.pos=T, 
    	assay="SCT", test.use = "LR", latent.vars = "smoke")
}

#Secondly, do DE for each rare cell against all non-rare cells
#For this, need to bring in T15_int
Idents(T15_int) <- T15_int@meta.data$clusters_18
markers.list[[length(markers.list) + 1]] <- FindMarkers(T15_int, ident.1="c6a", ident.2=
	unique(T15_int@meta.data$clusters_18)[-grep("c6",unique(T15_int@meta.data$clusters_18))],
	min.pct=0.1, logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "LR", latent.vars = "smoke")
markers.list[[length(markers.list) + 1]] <- FindMarkers(T15_int, ident.1="c6b", ident.2=
	unique(T15_int@meta.data$clusters_18)[-grep("c6",unique(T15_int@meta.data$clusters_18))],
	min.pct=0.1, logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "LR", latent.vars = "smoke")
markers.list[[length(markers.list) + 1]] <- FindMarkers(T15_int, ident.1="c6c", ident.2=
	unique(T15_int@meta.data$clusters_18)[-grep("c6",unique(T15_int@meta.data$clusters_18))],
	min.pct=0.1, logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "LR", latent.vars = "smoke")

#Combine tables
comparisonVec<-c("Ion_vs_rare","Tuft_vs_rare","PNEC_vs_rare","Ion_vs_nonrare","Tuft_vs_nonrare","PNEC_vs_nonrare")

#Get vector of clusters for calculating pct values from RNA
pctIndexList<-list()
for(i in 1:length(clusterList)){
	currClusterIndex<-which(clusterList == clusterList[i])
	pctIndexList[[length(pctIndexList) + 1]]<-list(clusterList[i],clusterList[-currClusterIndex])
}
pctIndexList[[length(pctIndexList) + 1]]<-list("c6a",unique(T15_int@meta.data$clusters_18)[-grep("c6",unique(T15_int@meta.data$clusters_18))])
pctIndexList[[length(pctIndexList) + 1]]<-list("c6b",unique(T15_int@meta.data$clusters_18)[-grep("c6",unique(T15_int@meta.data$clusters_18))])
pctIndexList[[length(pctIndexList) + 1]]<-list("c6c",unique(T15_int@meta.data$clusters_18)[-grep("c6",unique(T15_int@meta.data$clusters_18))])

#Get RNA counts and cluster assignments
counts<-as.matrix(GetAssayData(T15_int,assay="RNA",slot="data"))
clusterAssignments<-T15_int@meta.data$clusters_18

#Run OrganizeDF_list
Rare3_DEGs<-OrganizeDF_list(datasetList=markers.list,comparisonVec=comparisonVec,pctIndexList=pctIndexList,
	counts=counts,clusterAssignments=clusterAssignments)

#Write all to single txt file
write.table(Rare3_DEGs,sep="\t",quote=FALSE,row.names=FALSE,file="DEG_Tables/Rare3_DEGs.txt")

#Write to excel
writeDEGsToExcel(DEG_table = Rare3_DEGs, savedFile = "DEG_Tables/Rare3_DEGs.xlsx")

#Export table for Enrichr
Rare3_DEGs_forEnrichr<-Rare3_DEGs[which(Rare3_DEGs$p_adj_FDR < 0.05),] #Only good DEGs

#Write all to single txt file
write.table(Rare3_DEGs_forEnrichr,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,
	file="Enrichr/Rare3_DEGs_forEnrichr.txt")

#Do enrichments
DEG_table<-read.table("Enrichr/Rare3_DEGs_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"Rare3_DEGs"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)



















#Finally, merge the log fold changes and FDRs from across the six tables above.
temp<-Rare3_DEGs_forEnrichr[,c(1:2,4,9)]
colnames(temp)<-c("gene","LFC","FDR","comparison")

#Make list of the six tables    
DEGList<-list(temp[which(temp$comparison==unique(temp$comparison)[1]),],
	temp[which(temp$comparison==unique(temp$comparison)[2]),],temp[which(temp$comparison==unique(temp$comparison)[3]),],
	temp[which(temp$comparison==unique(temp$comparison)[4]),],temp[which(temp$comparison==unique(temp$comparison)[5]),],
	temp[which(temp$comparison==unique(temp$comparison)[6]),])	

#Loop through merging
rare_merged<-DEGList[[1]]
for(i in 2:length(unique(temp$comparison))){
	rare_merged<-merge(rare_merged,DEGList[[i]],by="gene",all=T,
		suffixes=c(paste(".",unique(as.character(temp$comparison))[i-1],sep=""),
		paste(".",unique(as.character(temp$comparison))[i],sep="")))
}

#Add a column that lists all the comparisons significant for a gene
rare_merged<-cbind(rare_merged,"Which_DEG"=apply(rare_merged[,grep("comparison",colnames(rare_merged))],1,
	function(x)paste(x[which(!is.na(x))],collapse=",")))

#Filter out the comparison columns
rare_merged<-rare_merged[,-grep("comparison",colnames(rare_merged))]

#Write to xcel
write.xlsx(rare_merged,file="Combined_rareCell_DEG_table.xlsx")
	
	
	
	
	
	





























################################## Repeat differential expression above, just using the 6 nonsmokers #########################################

#Bring in T15_int and then define T15 to include only smokers
T15_ns<-subset(T15_int,cells=rownames(T15_int@meta.data)[which(T15_int@meta.data$smoke_noT89 == "never")])

#Do DE (each against the other rare cells)
Idents(T15_ns) <- T15_ns@meta.data$clusters_16a
clusterList<-c("c6a","c6b","c6c")
markers.list <- list()
for(i in 1:length(clusterList)) {
    markers.list[[i]] <- FindMarkers(T15_ns, ident.1=clusterList[i], ident.2=clusterList[-i], min.pct=0.1, logfc.threshold=0.25, only.pos=T, 
    	assay="SCT", test.use = "LR", latent.vars = "donor")
}

#Secondly, do DE for each rare cell against all non-rare cells
#For this, need to bring in T15_ns
Idents(T15_ns) <- T15_ns@meta.data$clusters_16a
markers.list[[length(markers.list) + 1]] <- FindMarkers(T15_ns, ident.1="c6a", ident.2=
	unique(T15_ns@meta.data$clusters_16a)[-grep("c6",unique(T15_ns@meta.data$clusters_16a))],
	min.pct=0.1, logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "LR", latent.vars = "donor")
markers.list[[length(markers.list) + 1]] <- FindMarkers(T15_ns, ident.1="c6b", ident.2=
	unique(T15_ns@meta.data$clusters_16a)[-grep("c6",unique(T15_ns@meta.data$clusters_16a))],
	min.pct=0.1, logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "LR", latent.vars = "donor")
markers.list[[length(markers.list) + 1]] <- FindMarkers(T15_ns, ident.1="c6c", ident.2=
	unique(T15_ns@meta.data$clusters_16a)[-grep("c6",unique(T15_ns@meta.data$clusters_16a))],
	min.pct=0.1, logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "LR", latent.vars = "donor")

#Combine tables
comparisonVec<-c("Ion_vs_rare","Tuft_vs_rare","PNEC_vs_rare","Ion_vs_nonrare","Tuft_vs_nonrare","PNEC_vs_nonrare")

#Get vector of clusters for calculating pct values from RNA
pctIndexList<-list()
for(i in 1:length(clusterList)){
	currClusterIndex<-which(clusterList == clusterList[i])
	pctIndexList[[length(pctIndexList) + 1]]<-list(clusterList[i],clusterList[-currClusterIndex])
}
pctIndexList[[length(pctIndexList) + 1]]<-list("c6a",unique(T15_ns@meta.data$clusters_16a)[-grep("c6",unique(T15_ns@meta.data$clusters_16a))])
pctIndexList[[length(pctIndexList) + 1]]<-list("c6b",unique(T15_ns@meta.data$clusters_16a)[-grep("c6",unique(T15_ns@meta.data$clusters_16a))])
pctIndexList[[length(pctIndexList) + 1]]<-list("c6c",unique(T15_ns@meta.data$clusters_16a)[-grep("c6",unique(T15_ns@meta.data$clusters_16a))])

#Get RNA counts and cluster assignments
counts<-as.matrix(GetAssayData(T15_ns,assay="RNA",slot="data"))
clusterAssignments<-T15_ns@meta.data$clusters_16a

#Run OrganizeDF_list
Rare3_nonsmokerOnly_DEGs<-OrganizeDF_list(datasetList=markers.list,comparisonVec=comparisonVec,pctIndexList=pctIndexList,
	counts=counts,clusterAssignments=clusterAssignments)

#Write all to single txt file
write.table(Rare3_nonsmokerOnly_DEGs,sep="\t",quote=FALSE,row.names=FALSE,file="DEG_Tables/Rare3_nonsmokerOnly_DEGs.txt")

#Write to excel
writeDEGsToExcel(DEG_table = Rare3_nonsmokerOnly_DEGs, savedFile = "DEG_Tables/Rare3_nonsmokerOnly_DEGs.xlsx")

#Export table for Enrichr
Rare3_nonsmokerOnly_DEGs_forEnrichr<-Rare3_nonsmokerOnly_DEGs[which(Rare3_nonsmokerOnly_DEGs$p_adj_FDR < 0.05),] #Only good DEGs

#Write all to single txt file
write.table(Rare3_nonsmokerOnly_DEGs_forEnrichr,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,
	file="Enrichr/Rare3_nonsmokerOnly_DEGs_forEnrichr.txt")

#Do enrichments
DEG_table<-read.table("Enrichr/Rare3_nonsmokerOnly_DEGs_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"Rare3_nonsmokerOnly_DEGs"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)






#Finally, merge the log fold changes and FDRs from across the six tables above.
temp<-Rare3_nonsmokerOnly_DEGs_forEnrichr[,c(1:2,4,9)]
colnames(temp)<-c("gene","LFC","FDR","comparison")

#Make list of the six tables    
DEGList<-list(temp[which(temp$comparison==unique(temp$comparison)[1]),],
	temp[which(temp$comparison==unique(temp$comparison)[2]),],temp[which(temp$comparison==unique(temp$comparison)[3]),],
	temp[which(temp$comparison==unique(temp$comparison)[4]),],temp[which(temp$comparison==unique(temp$comparison)[5]),],
	temp[which(temp$comparison==unique(temp$comparison)[6]),])	

#Loop through merging
rare_nonsmokerOnlly_merged<-DEGList[[1]]
for(i in 2:length(unique(temp$comparison))){
	rare_nonsmokerOnlly_merged<-merge(rare_nonsmokerOnlly_merged,DEGList[[i]],by="gene",all=T,
		suffixes=c(paste(".",unique(as.character(temp$comparison))[i-1],sep=""),
		paste(".",unique(as.character(temp$comparison))[i],sep="")))
}

#Add a column that lists all the comparisons significant for a gene
rare_nonsmokerOnlly_merged<-cbind(rare_nonsmokerOnlly_merged,"Which_DEG"=apply(rare_nonsmokerOnlly_merged[,grep("comparison",colnames(rare_nonsmokerOnlly_merged))],1,
	function(x)paste(x[which(!is.na(x))],collapse=",")))

#Filter out the comparison columns
rare_nonsmokerOnlly_merged<-rare_nonsmokerOnlly_merged[,-grep("comparison",colnames(rare_nonsmokerOnlly_merged))]

#Write to xcel
write.xlsx(rare_nonsmokerOnlly_merged,colWidths="auto",firstRow=T,file="Combined_rareCell_nonsmokerOnly_DEG_table.xlsx")
	