######################################################SMOKING EFFECTS##############################################################
library(Seurat, lib.loc="/Users/jacksonna/R-dev")
library(plotrix)
library(scales)
library(ggplot2)
library(cowplot)
library(openxlsx)
library(venn)
library(heatmap3)
source("CommonFunctions.r")

#Download the processed Seurat R object, "Processed_invitro_seurat.Rdata" from GEO,
#obtainable from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134174
#Then load the object it contains ("T15_int") into the session.
#T15_int is a Seurat dataset that contains the processed, integrated dataset
load("Processed_invitro_seurat.Rdata")


############################# GET SMOKING DEGS and GET CORE AND UNIQUE
#bsub -e test.err -o test.out -R "rusage[mem=42000]" "module load R/3.5.1 && env R_MAX_VSIZE=700Gb R CMD BATCH diffexp_clusters10_smoke_up.r"
#bsub -e test.err -o test.out -R "rusage[mem=42000]" "module load R/3.5.1 && env R_MAX_VSIZE=700Gb R CMD BATCH diffexp_clusters10_smoke_down.r"

#Get combined cluster and smoking column for doing DE
treat.col<-data.frame("clusters10_smoke"=paste(T15_int@meta.data$clusters_10,T15_int@meta.data$smoke,sep="_"),row.names=rownames(T15_int@meta.data))
treat.col[,1]<-as.character(treat.col[,1])
#Get rid of T89
treat.col[which(T15_int@meta.data$donor == "T89"),]<-sapply(strsplit(treat.col[which(T15_int@meta.data$donor == "T89"),],"_"),function(x)paste(x[1],"exclude",sep="_"))
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="clusters10_smoke")

#Repeat this with the 16a
treat.col<-data.frame("clusters16a_smoke"=paste(T15_int@meta.data$clusters_16a,T15_int@meta.data$smoke,sep="_"),row.names=rownames(T15_int@meta.data))
treat.col[,1]<-as.character(treat.col[,1])
#Get rid of T89
treat.col[which(T15_int@meta.data$donor == "T89"),]<-sapply(strsplit(treat.col[which(T15_int@meta.data$donor == "T89"),],"_"),function(x)paste(x[1],"exclude",sep="_"))
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="clusters16a_smoke")

#Finally, while we're at it, add a new smoke column that excluded T89 from the never category
treat.col<-data.frame("smoke_noT89"=T15_int@meta.data$smoke,row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col[,1]<-as.character(treat.col[,1])
treat.col[rownames(T15_int@meta.data)[which(T15_int@meta.data$donor == "T89")],]<-"excluded"
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="smoke_noT89")



#Set identity
Idents(T15_int) <- T15_int@meta.data$clusters10_smoke

#Get list of comparisons to do
#UPREGULATED
DE_list_up<-list()
comparisons_up<-sort(as.character(unique(T15_int@meta.data$clusters10_smoke)))[-grep("light|exclude",
	sort(as.character(unique(T15_int@meta.data$clusters10_smoke))))]
for(i in seq(1,length(unique(T15_int@meta.data$integrated_snn_res.0.22)) * 2,by=2)){
	DE_list_up[[length(DE_list_up) + 1]]<-list(comparisons_up[i],comparisons_up[i+1])
}
#DOWNREGULATED
DE_list_down<-list()
comparisons_down<-sort(as.character(unique(T15_int@meta.data$clusters10_smoke)))[-grep("light|exclude",
	sort(as.character(unique(T15_int@meta.data$clusters10_smoke))))]
for(i in seq(1,length(unique(T15_int@meta.data$integrated_snn_res.0.22)) * 2,by=2)){
	DE_list_down[[length(DE_list_down) + 1]]<-list(comparisons_down[i+1],comparisons_down[i])
}




##Now do differential expression between smoker_heavy and nonsmoker for each cluster in cluster 10
#Did on cluster
#UPREGULATED
clusters10_smoke_markers_up<-list()
for(i in 1:length(DE_list_up)){
clusters10_smoke_markers_up[[length(clusters10_smoke_markers_up) + 1]]<-FindMarkers(T15_int,ident.1 = unlist(DE_list_up[[i]][1]),ident.2 = unlist(DE_list_up[[i]][2]), 
	min.pct = 0.1, only.pos=T, logfc.threshold = 0.25, assay = "SCT", test.use = "wilcox")
}
#DOWNREGULATED
clusters10_smoke_markers_down<-list()
for(i in 1:length(DE_list_down)){
clusters10_smoke_markers_down[[length(clusters10_smoke_markers_down) + 1]]<-FindMarkers(T15_int,ident.1 = unlist(DE_list_down[[i]][1]),ident.2 = unlist(DE_list_down[[i]][2]), 
	min.pct = 0.1, only.pos=T, logfc.threshold = 0.25, assay = "SCT", test.use = "wilcox")
}



#######PROCESS DEGS
#If bringing in from rda files:
load("diffexp_rda/clusters10_smoke_markers_up.rda")
load("diffexp_rda/clusters10_smoke_markers_down.rda")

#UPREGULATED
#Get comparison descriptions
comparisonVec<-sapply(strsplit(sapply(DE_list_up,function(x)paste(x,collapse="_")),"_"),function(x)x[1])

#Get vector of clusters for calculating pct values
pctIndexList<-DE_list_up

#Get counts and cluster assignments
counts<-as.matrix(GetAssayData(T15_int,assay="RNA",slot="data"))
clusterAssignments<-T15_int@meta.data$clusters10_smoke

#Run OrganizeDF_list
T15_clusters10_smoke_up_DEGs<-OrganizeDF_list(datasetList=clusters10_smoke_markers_up,comparisonVec=comparisonVec,
	pctIndexList=pctIndexList,counts=counts,clusterAssignments=clusterAssignments)

#Write all to single txt file
write.table(T15_clusters10_smoke_up_DEGs,sep="\t",quote=FALSE,row.names=FALSE,file="DEG_Tables/T15_clusters10_smoke_up_DEGs.txt")

#Write to excel
wb<-createWorkbook()
for(i in 1:length(unique(T15_clusters10_smoke_up_DEGs$comparison))){
 	addWorksheet(wb,as.character(unique(T15_clusters10_smoke_up_DEGs$comparison)[i]))
 	writeData(wb,i,T15_clusters10_smoke_up_DEGs[which(T15_clusters10_smoke_up_DEGs$comparison == unique(T15_clusters10_smoke_up_DEGs$comparison)[i]),])
 	setColWidths(wb, sheet = i, cols = 1:ncol(currData), widths = "auto")
 	freezePane(wb,i,firstRow=T)
 	saveWorkbook(wb,"DEG_Tables/T15_clusters10_smoke_up_DEGs.xlsx", overwrite = TRUE)  	
}

#Export table for Enrichr
T15_clusters10_smoke_up_DEGs_forEnrichr<-T15_clusters10_smoke_up_DEGs[which(T15_clusters10_smoke_up_DEGs$p_adj_FDR < 0.05),] #Only good DEGs

#Write all to single txt file
write.table(T15_clusters10_smoke_up_DEGs_forEnrichr,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,file="Enrichr/T15_clusters10_smoke_up_DEGs_forEnrichr.txt")

#Do enrichments
DEG_table<-read.table("Enrichr/T15_clusters10_smoke_up_DEGs_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"T15_clusters10_smoke_up_DEGs"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)







#DOWNREGULATED
#Get comparison descriptions
comparisonVec<-sapply(strsplit(sapply(DE_list_down,function(x)paste(x,collapse="_")),"_"),function(x)x[1])

#Get vector of clusters for calculating pct values
pctIndexList<-DE_list_down

#Get counts and cluster assignments
counts<-as.matrix(GetAssayData(T15_int,assay="RNA",slot="data"))
clusterAssignments<-T15_int@meta.data$clusters10_smoke

#Run OrganizeDF_list
T15_clusters10_smoke_down_DEGs<-OrganizeDF_list(datasetList=clusters10_smoke_markers_down,comparisonVec=comparisonVec,
	pctIndexList=pctIndexList,counts=counts,clusterAssignments=clusterAssignments)

#Write all to single txt file
write.table(T15_clusters10_smoke_down_DEGs,sep="\t",quote=FALSE,row.names=FALSE,file="DEG_Tables/T15_clusters10_smoke_down_DEGs.txt")

#Write to excel
wb<-createWorkbook()
for(i in 1:length(unique(T15_clusters10_smoke_down_DEGs$comparison))){
 	addWorksheet(wb,as.character(unique(T15_clusters10_smoke_down_DEGs$comparison)[i]))
 	writeData(wb,i,T15_clusters10_smoke_down_DEGs[which(T15_clusters10_smoke_down_DEGs$comparison == unique(T15_clusters10_smoke_down_DEGs$comparison)[i]),])
 	setColWidths(wb, sheet = i, cols = 1:ncol(currData), widths = "auto")
 	freezePane(wb,i,firstRow=T)
 	saveWorkbook(wb,"DEG_Tables/T15_clusters10_smoke_down_DEGs.xlsx", overwrite = TRUE)  	
}

#Export table for Enrichr
T15_clusters10_smoke_down_DEGs_forEnrichr<-T15_clusters10_smoke_down_DEGs[which(T15_clusters10_smoke_down_DEGs$p_adj_FDR < 0.05),] #Only good DEGs

#Write all to single txt file
write.table(T15_clusters10_smoke_down_DEGs_forEnrichr,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,file="Enrichr/T15_clusters10_smoke_down_DEGs_forEnrichr.txt")

#Do enrichments
DEG_table<-read.table("Enrichr/T15_clusters10_smoke_down_DEGs_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"T15_clusters10_smoke_down_DEGs"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)




















###################################GET CORE SMOKING DEGs
############# Also, get DEGs that are shared among all, 7, 6, 5abc, 4, 3, 2, and 1 clusters
#Isolate the upregulated comparisons
##NOTE: for getting the "core" response, we're going to ignore the rare cell populations
###UP
currData<-T15_clusters10_smoke_up_DEGs_forEnrichr[-grep("c6",T15_clusters10_smoke_up_DEGs_forEnrichr$comparison),]
currData$comparison<-droplevels(currData$comparison)
###DOWN
currData<-T15_clusters10_smoke_down_DEGs_forEnrichr[-grep("c6",T15_clusters10_smoke_down_DEGs_forEnrichr$comparison),]
currData$comparison<-droplevels(currData$comparison)

#Store the genes
sharedDEGsList<-list()

#######Run the following for either UP or DOWN genes (takes ~10 min)
#Count the genes for each replicate
numDEGs<-c()
listCount<-1
for(j in 10:1){
	sharedDEGs<-c()
	for(i in 1:length(unique(currData$gene))){
		if(length(currData$gene[which(currData$gene == unique(currData$gene)[i])]) >= j){
			sharedDEGs<-append(sharedDEGs,unique(currData$gene)[i])
		}
	}
	if(!is.null(sharedDEGs)){
		sharedDEGsList[[length(sharedDEGsList) + 1]]<-sharedDEGs
		names(sharedDEGsList)[listCount]<-paste(j,"clusters",sep="_")
		numDEGs<-append(numDEGs,length(sharedDEGs))
		names(numDEGs)[listCount]<-paste(j,"clusters",sep="_")
		listCount<-listCount + 1
	}
}

##Barplot how DEG sharing is distributed across clusters
#pdf("T15_clusters10_smoke_up_distributionOfDEGsharing_barplot.pdf")
dev.new(height=4,width=6.5)
barplot(numDEGs,col=c("white"),las=1,ylab="Number of DEGs",names=as.character(c(8:1)),xlab="Number of clusters sharing a DEG")
#pdf("T15_clusters10_smoke_down_distributionOfDEGsharing_barplot.pdf")
dev.new(height=4,width=6.5)
barplot(numDEGs,col=c("white"),las=1,ylab="Number of DEGs",names=as.character(c(8:1)),xlab="Number of clusters sharing a DEG")


####Finally, make a table for just the "core" genes (with 5 or more clusters being significant),
####itemizing which clusters were involved
core_clusters<-data.frame("core_gene"=sharedDEGsList$'5_clusters',"clusters"=NA)
core_clusters<-core_clusters[order(core_clusters$core_gene),]
for(i in 1:nrow(core_clusters)){
	currVec<-as.character(currData[which(currData$gene %in% core_clusters$core_gene[i]),]$comparison)
	core_clusters$clusters[i]<-paste(currVec,collapse="_")
}
write.xlsx(core_clusters, file = "T15_clusters10_smoke_up_coreGenes.xlsx")
write.xlsx(core_clusters, file = "T15_clusters10_smoke_down_coreGenes.xlsx")

##Also, print these DEGs to excel file
write.xlsx(sharedDEGsList, file = "T15_clusters10_smoke_up_distributionOfDEGsharing.xlsx")
write.xlsx(sharedDEGsList, file = "T15_clusters10_smoke_down_distributionOfDEGsharing.xlsx")



#What's the overlap of these "core" genes with the old set 
old_up<-read.xlsx("Gene_lists/T0_all_twoSmokingLevels_coreGenes_c5abcVersion_up.xlsx")
old_down<-read.xlsx("Gene_lists/T0_all_twoSmokingLevels_coreGenes_c5abcVersion_down.xlsx")

new_up_5<-read.xlsx("IMAGES_smoking_effects/T15_clusters10_smoke_up_distributionOfDEGsharing.xlsx",colNames=F,sheet = 4)
new_down_5<-read.xlsx("IMAGES_smoking_effects/T15_clusters10_smoke_down_distributionOfDEGsharing.xlsx",colNames=F,sheet = 4)

#Overlapping core genes between old and new
length(intersect(new_up_4[,1],old_up[,1])) #14 overlapping
length(intersect(new_down_5[,1],old_down[,1])) #7 overlapping



#Bring core genes in and run enricher on them
core_up<-read.xlsx("IMAGES_smoking_effects/T15_clusters10_smoke_up_coreGenes.xlsx")$core_gene
core_down<-read.xlsx("IMAGES_smoking_effects/T15_clusters10_smoke_down_coreGenes.xlsx")$core_gene
#Build table for Enrichr
core_genes<-data.frame("gene"=c(core_up,core_down),"comparison"=c(rep("up",length(core_up)),rep("down",length(core_down))),stringsAsFactors=F)
write.table(core_genes,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,file="Enrichr/coreSmokingDEGs_forEnrichr.txt")
#Do enrichments
DEG_table<-read.table("Enrichr/coreSmokingDEGs_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"coreSmokingDEGs"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)






















######################################################################################
########################################## Get tables of unique and core genes for each cluster
######################################################################################
#Read in previously produced smoking DEGs
smoke_up<-read.table("DEG_Tables/T15_clusters10_smoke_up_DEGs.txt",header=T,stringsAsFactors=F)
smoke_down<-read.table("DEG_Tables/T15_clusters10_smoke_down_DEGs.txt",header=T,stringsAsFactors=F)

#Read in previously produced core DEGs
core_up<-read.xlsx("IMAGES_smoking_effects/T15_clusters10_smoke_up_coreGenes.xlsx")$core_gene
core_down<-read.xlsx("IMAGES_smoking_effects/T15_clusters10_smoke_down_coreGenes.xlsx")$core_gene


#Also, get rare cell smoking DEGs
#First, define our 16a clusters
#Bring in the monocle-derived mature secretory cell group and add as metadata
matureSecCells<-scan("monocle/ANALYSIS_lumSec_reintegratedData/Mature_secretory_cells_monocle.txt",what="character")
treat.col<-data.frame("Mature_secretory_cells_monocle"=rep(FALSE,nrow(T15_int@meta.data)),row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col[which(rownames(treat.col) %in% matureSecCells),]<-TRUE
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="Mature_secretory_cells_monocle")
#Also, bring in smgb, ciliated, and rare cell subclusters
load("Subclustering_rda/T15_cil.rda")
load("Subclustering_rda/T15_smgb.rda")
load("Subclustering_rda/T15_rare3.rda")
#Now combine subclusters and clusters I care about together 
treat.col<-data.frame("clusters_16a"=T15_int@meta.data$clusters_10,row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col[,1]<-as.character(treat.col[,1])
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_smgb@meta.data)[which(T15_smgb@meta.data$clusters_3 == "c3a")])],]<-"c3a"
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_smgb@meta.data)[which(T15_smgb@meta.data$clusters_3 == "c3b")])],]<-"c3b"
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_smgb@meta.data)[which(T15_smgb@meta.data$clusters_3 == "c3c")])],]<-"c3c"
treat.col[rownames(T15_int@meta.data)[which(T15_int@meta.data$Mature_secretory_cells_monocle==TRUE)],]<-"c5M"
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_rare3@meta.data)[which(T15_rare3@meta.data$clusters_3 == "c6a")])],]<-"c6a"
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_rare3@meta.data)[which(T15_rare3@meta.data$clusters_3 == "c6b")])],]<-"c6b"
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_rare3@meta.data)[which(T15_rare3@meta.data$clusters_3 == "c6c")])],]<-"c6c"
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_cil@meta.data)[which(T15_cil@meta.data$clusters_3 == "c8a")])],]<-"c8a"
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_cil@meta.data)[which(T15_cil@meta.data$clusters_3 == "c8b")])],]<-"c8b"
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_cil@meta.data)[which(T15_cil@meta.data$clusters_3 == "c8c")])],]<-"c8c"
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="clusters_16a")
#Now, set up smoking status column for each of the 16a clusters
treat.col<-data.frame("clusters16a_smoke"=paste(T15_int@meta.data$clusters_16a,T15_int@meta.data$smoke,sep="_"),row.names=rownames(T15_int@meta.data))
treat.col[,1]<-as.character(treat.col[,1])
#Get rid of T89, who is too young
treat.col[which(T15_int@meta.data$donor == "T89"),]<-sapply(strsplit(treat.col[which(T15_int@meta.data$donor == "T89"),],"_"),function(x)paste(x[1],"exclude",sep="_"))
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="clusters16a_smoke")
#Now get up smoking DEGs in rare3 groups
Idents(T15_int)<-T15_int@meta.data$clusters16a_smoke
c6a_smoke_up<-FindMarkers(T15_int,ident.1 = c("c6a_heavy"),ident.2="c6a_never",min.pct = 0.1, only.pos=T, logfc.threshold = 0.25, assay = "SCT", test.use = "wilcox")
c6b_smoke_up<-FindMarkers(T15_int,ident.1 = c("c6b_heavy"),ident.2="c6b_never",min.pct = 0.1, only.pos=T, logfc.threshold = 0.25, assay = "SCT", test.use = "wilcox")
c6c_smoke_up<-FindMarkers(T15_int,ident.1 = c("c6c_heavy"),ident.2="c6c_never",min.pct = 0.1, only.pos=T, logfc.threshold = 0.25, assay = "SCT", test.use = "wilcox")
#Now get down smoking DEGs in rare3 groups
c6a_smoke_down<-FindMarkers(T15_int,ident.1 = c("c6a_never"),ident.2="c6a_heavy",min.pct = 0.1, only.pos=T, logfc.threshold = 0.25, assay = "SCT", test.use = "wilcox")
c6b_smoke_down<-FindMarkers(T15_int,ident.1 = c("c6b_never"),ident.2="c6b_heavy",min.pct = 0.1, only.pos=T, logfc.threshold = 0.25, assay = "SCT", test.use = "wilcox")
c6c_smoke_down<-FindMarkers(T15_int,ident.1 = c("c6c_never"),ident.2="c6c_heavy",min.pct = 0.1, only.pos=T, logfc.threshold = 0.25, assay = "SCT", test.use = "wilcox")

#Now organize into tables
comparisonVec<-c("c6a_smoke_up","c6b_smoke_up","c6c_smoke_up","c6a_smoke_down","c6b_smoke_down","c6c_smoke_down")
#Run OrganizeDF_list
c6_smoke_stratified_DEGs<-OrganizeDF_list(datasetList=list(c6a_smoke_up,c6b_smoke_up,c6c_smoke_up,c6a_smoke_down,
	c6b_smoke_down,c6c_smoke_down),comparisonVec=comparisonVec)

#Write all to single txt file
write.table(c6_smoke_stratified_DEGs,sep="\t",quote=FALSE,row.names=FALSE,file="DEG_Tables/c6_smoke_stratified_DEGs.txt")

#Write to excel
writeDEGsToExcel(c6_smoke_stratified_DEGs,"DEG_Tables/c6_smoke_stratified_DEGs.xlsx")







########################## And get smoking DEGs in ciliated cells, stratified into hybrid and non-hybrid groups
#Now get up smoking DEGs in non-hybrid cells 
Idents(T15_int)<-T15_int@meta.data$clusters16a_smoke
c8bc_smoke_up<-FindMarkers(T15_int,ident.1 = c("c8b_heavy","c8c_heavy"),ident.2 = c("c8b_never","c8c_never"),min.pct = 0.1, 
	only.pos=T, logfc.threshold = 0.25, assay = "SCT", test.use = "wilcox")
c8bc_smoke_down<-FindMarkers(T15_int,ident.1 = c("c8b_never","c8c_never"),ident.2 = c("c8b_heavy","c8c_heavy"),min.pct = 0.1, 
	only.pos=T, logfc.threshold = 0.25, assay = "SCT", test.use = "wilcox")
#Now get up smoking DEGs in hybrid cells 
c8a_smoke_up<-FindMarkers(T15_int,ident.1 = "c8a_heavy",ident.2 = "c8a_never",min.pct = 0.1, 
	only.pos=T, logfc.threshold = 0.25, assay = "SCT", test.use = "wilcox")
c8a_smoke_down<-FindMarkers(T15_int,ident.1 = "c8a_never",ident.2 = "c8a_heavy",min.pct = 0.1, 
	only.pos=T, logfc.threshold = 0.25, assay = "SCT", test.use = "wilcox")

#Now organize into tables
comparisonVec<-c("c8bc_smoke_up","c8bc_smoke_down","c8a_smoke_up","c8a_smoke_down")
#Run OrganizeDF_list
c8_smoke_stratified_DEGs<-OrganizeDF_list(datasetList=list(c8bc_smoke_up,c8bc_smoke_down,c8a_smoke_up,c8a_smoke_down),
	comparisonVec=comparisonVec)

#Write all to single txt file
write.table(c8_smoke_stratified_DEGs,sep="\t",quote=FALSE,row.names=FALSE,file="DEG_Tables/c8_smoke_stratified_DEGs.txt")

#Write to excel
writeDEGsToExcel(c8_smoke_stratified_DEGs,"DEG_Tables/c8_smoke_stratified_DEGs.xlsx")








####### Now combine together all the smoking DEGs
smoke_up_withSubclusters<-rbind(smoke_up[,-grep("counts",colnames(smoke_up))],c6_smoke_stratified_DEGs[grep("up",c6_smoke_stratified_DEGs$comparison),],
	c8_smoke_stratified_DEGs[grep("up",c8_smoke_stratified_DEGs$comparison),])
smoke_down_withSubclusters<-rbind(smoke_down[,-grep("counts",colnames(smoke_down))],c6_smoke_stratified_DEGs[grep("down",c6_smoke_stratified_DEGs$comparison),],
	c8_smoke_stratified_DEGs[grep("down",c8_smoke_stratified_DEGs$comparison),])






########## Now add columns to define the core/unique status of each gene under four different collections of populations:
#1. main - Only main pops (no rare)
#2. main_rare - Main with rare broken out
#3. main_cil + Main with cilia hybrids broken out
#5. main_rare_cil - Main with ciliated broken out + rare broken out
smoke_up_withSubclusters<-cbind(smoke_up_withSubclusters,"main"=NA,"main_rare"=NA,"main_cil"=NA,"main_rare_cil"=NA)
smoke_down_withSubclusters<-cbind(smoke_down_withSubclusters,"main"=NA,"main_rare"=NA,"main_cil"=NA,"main_rare_cil"=NA)





########## Get rid of c6ab and c6c since these play no role in these comparisons
smoke_up_withSubclusters<-smoke_up_withSubclusters[-grep("^c6ab$|^c6c$",smoke_up_withSubclusters$comparison),]
smoke_down_withSubclusters<-smoke_down_withSubclusters[-grep("^c6ab$|^c6c$",smoke_down_withSubclusters$comparison),]

########## Also, simplify the comparison names
smoke_up_withSubclusters$comparison<-sapply(strsplit(smoke_up_withSubclusters$comparison,"_"),function(x)x[1])
smoke_down_withSubclusters$comparison<-sapply(strsplit(smoke_down_withSubclusters$comparison,"_"),function(x)x[1])

######### And then sort the table by comparison
smoke_up_withSubclusters<-smoke_up_withSubclusters[order(smoke_up_withSubclusters$comparison),]
smoke_down_withSubclusters<-smoke_down_withSubclusters[order(smoke_down_withSubclusters$comparison),]

######## Then reset all the rownames to access later on
rownames(smoke_up_withSubclusters)<-seq(nrow(smoke_up_withSubclusters))
rownames(smoke_down_withSubclusters)<-seq(nrow(smoke_down_withSubclusters))



######## Now create four different datasets that are culled to fit the criteria above
smoke_up_main<-smoke_up_withSubclusters[-grep("^c6|c8a|c8bc",smoke_up_withSubclusters$comparison),]
smoke_up_main_rare<-smoke_up_withSubclusters[-grep("c8a|c8bc",smoke_up_withSubclusters$comparison),]
smoke_up_main_cil<-smoke_up_withSubclusters[-grep("^c6|c8$",smoke_up_withSubclusters$comparison),]
smoke_up_main_rare_cil<-smoke_up_withSubclusters[-grep("^c8$",smoke_up_withSubclusters$comparison),]

smoke_down_main<-smoke_down_withSubclusters[-grep("^c6|c8a|c8bc",smoke_down_withSubclusters$comparison),]
smoke_down_main_rare<-smoke_down_withSubclusters[-grep("c8a|c8bc",smoke_down_withSubclusters$comparison),]
smoke_down_main_cil<-smoke_down_withSubclusters[-grep("^c6|c8$",smoke_down_withSubclusters$comparison),]
smoke_down_main_rare_cil<-smoke_down_withSubclusters[-grep("^c8$",smoke_down_withSubclusters$comparison),]





######## Now, for each of the four datasets, get all the unique DEGs (FDR < 0.05) that are > 0.2 FDR in all other comparisons
#First lists of DEGs with FDR < 0.2 (these will be treated as DEGs for our purposes of determining uniqueness)
smoke_up_main_0.2<-smoke_up_main[which(smoke_up_main$p_adj_FDR < 0.2),]
smoke_up_main_rare_0.2<-smoke_up_main_rare[which(smoke_up_main_rare$p_adj_FDR < 0.2),]
smoke_up_main_cil_0.2<-smoke_up_main_cil[which(smoke_up_main_cil$p_adj_FDR < 0.2),]
smoke_up_main_rare_cil_0.2<-smoke_up_main_rare_cil[which(smoke_up_main_rare_cil$p_adj_FDR < 0.2),]

smoke_down_main_0.2<-smoke_down_main[which(smoke_down_main$p_adj_FDR < 0.2),]
smoke_down_main_rare_0.2<-smoke_down_main_rare[which(smoke_down_main_rare$p_adj_FDR < 0.2),]
smoke_down_main_cil_0.2<-smoke_down_main_cil[which(smoke_down_main_cil$p_adj_FDR < 0.2),]
smoke_down_main_rare_cil_0.2<-smoke_down_main_rare_cil[which(smoke_down_main_rare_cil$p_adj_FDR < 0.2),]


#Now get unique genes (those showing up only once)
smoke_up_main_unique<-smoke_up_main_0.2[which(!(smoke_up_main_0.2$gene %in% 
	smoke_up_main_0.2$gene[which(duplicated(smoke_up_main_0.2$gene))]) & 
	smoke_up_main_0.2$p_adj_FDR < 0.05),]
smoke_up_main_rare_unique<-smoke_up_main_rare_0.2[which(!(smoke_up_main_rare_0.2$gene %in% 
	smoke_up_main_rare_0.2$gene[which(duplicated(smoke_up_main_rare_0.2$gene))]) & 
	smoke_up_main_rare_0.2$p_adj_FDR < 0.05),]
smoke_up_main_cil_unique<-smoke_up_main_cil_0.2[which(!(smoke_up_main_cil_0.2$gene %in% 
	smoke_up_main_cil_0.2$gene[which(duplicated(smoke_up_main_cil_0.2$gene))]) & 
	smoke_up_main_cil_0.2$p_adj_FDR < 0.05),]
smoke_up_main_rare_cil_unique<-smoke_up_main_rare_cil_0.2[which(!(smoke_up_main_rare_cil_0.2$gene %in% 
	smoke_up_main_rare_cil_0.2$gene[which(duplicated(smoke_up_main_rare_cil_0.2$gene))]) & 
	smoke_up_main_rare_cil_0.2$p_adj_FDR < 0.05),]

smoke_down_main_unique<-smoke_down_main_0.2[which(!(smoke_down_main_0.2$gene %in% 
	smoke_down_main_0.2$gene[which(duplicated(smoke_down_main_0.2$gene))]) & 
	smoke_down_main_0.2$p_adj_FDR < 0.05),]
smoke_down_main_rare_unique<-smoke_down_main_rare_0.2[which(!(smoke_down_main_rare_0.2$gene %in% 
	smoke_down_main_rare_0.2$gene[which(duplicated(smoke_down_main_rare_0.2$gene))]) & 
	smoke_down_main_rare_0.2$p_adj_FDR < 0.05),]
smoke_down_main_cil_unique<-smoke_down_main_cil_0.2[which(!(smoke_down_main_cil_0.2$gene %in% 
	smoke_down_main_cil_0.2$gene[which(duplicated(smoke_down_main_cil_0.2$gene))]) & 
	smoke_down_main_cil_0.2$p_adj_FDR < 0.05),]
smoke_down_main_rare_cil_unique<-smoke_down_main_rare_cil_0.2[which(!(smoke_down_main_rare_cil_0.2$gene %in% 
	smoke_down_main_rare_cil_0.2$gene[which(duplicated(smoke_down_main_rare_cil_0.2$gene))]) & 
	smoke_down_main_rare_cil_0.2$p_adj_FDR < 0.05),]






###### Now fill in the table using the original rowname positions
#Unique for unique genes
smoke_up_withSubclusters$main[as.numeric(rownames(smoke_up_main_unique))]<-"Unique"
smoke_up_withSubclusters$main_rare[as.numeric(rownames(smoke_up_main_rare_unique))]<-"Unique"
smoke_up_withSubclusters$main_cil[as.numeric(rownames(smoke_up_main_cil_unique))]<-"Unique"
smoke_up_withSubclusters$main_rare_cil[as.numeric(rownames(smoke_up_main_rare_cil_unique))]<-"Unique"

smoke_down_withSubclusters$main[as.numeric(rownames(smoke_down_main_unique))]<-"Unique"
smoke_down_withSubclusters$main_rare[as.numeric(rownames(smoke_down_main_rare_unique))]<-"Unique"
smoke_down_withSubclusters$main_cil[as.numeric(rownames(smoke_down_main_cil_unique))]<-"Unique"
smoke_down_withSubclusters$main_rare_cil[as.numeric(rownames(smoke_down_main_rare_cil_unique))]<-"Unique"


#Core for core genes
smoke_up_withSubclusters$main[which(smoke_up_withSubclusters$gene %in% core_up)]<-"Core"
smoke_up_withSubclusters$main_rare[which(smoke_up_withSubclusters$gene %in% core_up)]<-"Core"
smoke_up_withSubclusters$main_cil[which(smoke_up_withSubclusters$gene %in% core_up)]<-"Core"
smoke_up_withSubclusters$main_rare_cil[which(smoke_up_withSubclusters$gene %in% core_up)]<-"Core"

smoke_down_withSubclusters$main[which(smoke_down_withSubclusters$gene %in% core_down)]<-"Core"
smoke_down_withSubclusters$main_rare[which(smoke_down_withSubclusters$gene %in% core_down)]<-"Core"
smoke_down_withSubclusters$main_cil[which(smoke_down_withSubclusters$gene %in% core_down)]<-"Core"
smoke_down_withSubclusters$main_rare_cil[which(smoke_down_withSubclusters$gene %in% core_down)]<-"Core"


#Semi-unique for anything else
smoke_up_withSubclusters$main[which(is.na(smoke_up_withSubclusters$main))]<-"Semiunique"
smoke_up_withSubclusters$main_rare[which(is.na(smoke_up_withSubclusters$main_rare))]<-"Semiunique"
smoke_up_withSubclusters$main_cil[which(is.na(smoke_up_withSubclusters$main_cil))]<-"Semiunique"
smoke_up_withSubclusters$main_rare_cil[which(is.na(smoke_up_withSubclusters$main_rare_cil))]<-"Semiunique"

smoke_down_withSubclusters$main[which(is.na(smoke_down_withSubclusters$main))]<-"Semiunique"
smoke_down_withSubclusters$main_rare[which(is.na(smoke_down_withSubclusters$main_rare))]<-"Semiunique"
smoke_down_withSubclusters$main_cil[which(is.na(smoke_down_withSubclusters$main_cil))]<-"Semiunique"
smoke_down_withSubclusters$main_rare_cil[which(is.na(smoke_down_withSubclusters$main_rare_cil))]<-"Semiunique"







######### Return genes to NA in comparisons that aren't relevant to a given dataset
smoke_up_withSubclusters$main[grep("^c6|c8a|c8bc",smoke_up_withSubclusters$comparison)]<-NA
smoke_up_withSubclusters$main_rare[grep("c8a|c8bc",smoke_up_withSubclusters$comparison)]<-NA
smoke_up_withSubclusters$main_cil[grep("^c6|c8$",smoke_up_withSubclusters$comparison)]<-NA
smoke_up_withSubclusters$main_rare_cil[grep("^c8$",smoke_up_withSubclusters$comparison)]<-NA

smoke_down_withSubclusters$main[grep("^c6|c8a|c8bc",smoke_down_withSubclusters$comparison)]<-NA
smoke_down_withSubclusters$main_rare[grep("c8a|c8bc",smoke_down_withSubclusters$comparison)]<-NA
smoke_down_withSubclusters$main_cil[grep("^c6|c8$",smoke_down_withSubclusters$comparison)]<-NA
smoke_down_withSubclusters$main_rare_cil[grep("^c8$",smoke_down_withSubclusters$comparison)]<-NA





######## Cull out any genes that aren't significant DEGs at FDR < 0.05
smoke_up_withSubclusters<-smoke_up_withSubclusters[which(smoke_up_withSubclusters$p_adj_FDR < 0.05),]
smoke_down_withSubclusters<-smoke_down_withSubclusters[which(smoke_down_withSubclusters$p_adj_FDR < 0.05),]





######## Finally export tables
write.table(smoke_up_withSubclusters,quote=F,row.names=F,file="DEG_Tables/Smoking_up_DEGs_coreUniqueDesignations.txt")
writeDEGsToExcel(smoke_up_withSubclusters,"DEG_Tables/Smoking_up_DEGs_coreUniqueDesignations.xlsx")

write.table(smoke_down_withSubclusters,quote=F,row.names=F,file="DEG_Tables/Smoking_down_DEGs_coreUniqueDesignations.txt")
writeDEGsToExcel(smoke_down_withSubclusters,"DEG_Tables/Smoking_down_DEGs_coreUniqueDesignations.xlsx")














##########Post-processing
####UP DEGs
############# Now tally up number of unique DEGs for each of the four datasets
#Main
table(smoke_up_withSubclusters$comparison[which(smoke_up_withSubclusters$main == "Unique")])
# c1 c1c  c2  c3  c4  c5  c7  c8 
# 23 233   6  13  33  30  13 114 
#Main_rare
table(smoke_up_withSubclusters$comparison[which(smoke_up_withSubclusters$main_rare == "Unique")])
# c1 c1c  c2  c3  c4  c5 c6a c6b c6c  c7  c8 
# 20 171   4  12  13  22  16 594   8  11  59
#Main_cil
table(smoke_up_withSubclusters$comparison[which(smoke_up_withSubclusters$main_cil == "Unique")])
#  c1  c1c   c2   c3   c4   c5   c7 c8bc 
#  23  229    5   13   33   25   11  135
#Main_rare_cil
table(smoke_up_withSubclusters$comparison[which(smoke_up_withSubclusters$main_rare_cil == "Unique")])
#  c1  c1c   c2   c3   c4   c5  c6a  c6b  c6c   c7 c8bc 
#  20  168    4   12   13   19   16  585    8   10   76 


######### Now tally up the proportion of DEGs that are unique based on the four datasets
#Main
round(table(smoke_up_withSubclusters$comparison[which(smoke_up_withSubclusters$main == "Unique")])/
	table(smoke_up_withSubclusters$comparison)[-c(7,8,9,12,13)],3)
#   c1   c1c    c2    c3    c4    c5    c7    c8 
#0.095 0.476 0.054 0.080 0.163 0.203 0.232 0.509 
#Main_rare
round(table(smoke_up_withSubclusters$comparison[which(smoke_up_withSubclusters$main_rare == "Unique")])/
	table(smoke_up_withSubclusters$comparison)[-c(12,13)],3)
#   c1   c1c    c2    c3    c4    c5   c6a   c6b   c6c    c7    c8 
#0.083 0.350 0.036 0.074 0.064 0.149 0.327 0.774 0.471 0.196 0.263 
#Main_cil
round(table(smoke_up_withSubclusters$comparison[which(smoke_up_withSubclusters$main_cil == "Unique")])/
	table(smoke_up_withSubclusters$comparison)[-c(7,8,9,11,12)],3)
#   c1   c1c    c2    c3    c4    c5    c7  c8bc 
#0.095 0.468 0.045 0.080 0.163 0.169 0.196 0.540 
#Main_rare_cil
round(table(smoke_up_withSubclusters$comparison[which(smoke_up_withSubclusters$main_rare_cil == "Unique")])/
	table(smoke_up_withSubclusters$comparison)[-c(11,12)],3)
#   c1   c1c    c2    c3    c4    c5   c6a   c6b   c6c    c7  c8bc 
#0.083 0.344 0.036 0.074 0.064 0.128 0.327 0.763 0.471 0.179 0.304 


######### Now tally up the number and proportion of DEGs that are core (this will be the same across all four datasets)
#Main
round(table(smoke_up_withSubclusters$comparison[which(smoke_up_withSubclusters$main_rare == "Core" | 
	smoke_up_withSubclusters$main_cil == "Core")]))
#	  c1  c1c   c2   c3   c4   c5  c6a  c6b  c6c   c7   c8  c8a c8bc 
#  78   81   71   72   77   56   15   30    3   26   52    1   49
round(table(smoke_up_withSubclusters$comparison[which(smoke_up_withSubclusters$main_rare == "Core" | 
	smoke_up_withSubclusters$main_cil == "Core")]) / table(smoke_up_withSubclusters$comparison),3)
#   c1   c1c    c2    c3    c4    c5   c6a   c6b   c6c    c7    c8   c8a  c8bc 
#0.324 0.166 0.634 0.444 0.379 0.378 0.306 0.039 0.176 0.464 0.232 0.333 0.196 





####DOWN DEGs
############# Now tally down number of unique DEGs for each of the four datasets
#Main
table(smoke_down_withSubclusters$comparison[which(smoke_down_withSubclusters$main == "Unique")])
#  c1 c1c  c2  c3  c4  c5  c7  c8 
#   9  13  20  11   6  31  12  24 
#Main_rare
table(smoke_down_withSubclusters$comparison[which(smoke_down_withSubclusters$main_rare == "Unique")])
#  c1 c1c  c2  c3  c4  c5 c6a c6c  c7  c8  
#   9   9  13  11   3  26 757 194  10  18 
#Main_cil
table(smoke_down_withSubclusters$comparison[which(smoke_down_withSubclusters$main_cil == "Unique")])
#    c1  c1c   c2   c3   c4   c5   c7  c8a c8bc 
#     9   13   19   10    6   26   10  267   12 
#Main_rare_cil
table(smoke_down_withSubclusters$comparison[which(smoke_down_withSubclusters$main_rare_cil == "Unique")])
#  c1  c1c   c2   c3   c4   c5  c6a  c6c   c7  c8a c8bc  
#   9    9   13   10    3   22  736  187    8  245    8  


######### Now tally down the proportion of DEGs that are unique based on the four datasets
#Main
round(table(smoke_down_withSubclusters$comparison[which(smoke_down_withSubclusters$main == "Unique")])/
	table(smoke_down_withSubclusters$comparison)[-c(7,8,11,12)],3)
#   c1   c1c    c2    c3    c4    c5    c7    c8 
#0.153 0.255 0.274 0.239 0.102 0.326 0.255 0.545 
#Main_rare
round(table(smoke_down_withSubclusters$comparison[which(smoke_down_withSubclusters$main_rare == "Unique")])/
	table(smoke_down_withSubclusters$comparison)[-c(11,12)],3)
#   c1   c1c    c2    c3    c4    c5   c6a   c6c    c7    c8 
#0.153 0.176 0.178 0.239 0.051 0.274 0.907 0.840 0.213 0.409 
#Main_cil
round(table(smoke_down_withSubclusters$comparison[which(smoke_down_withSubclusters$main_cil == "Unique")])/
	table(smoke_down_withSubclusters$comparison)[-c(7,8,10)],3)
#   c1   c1c    c2    c3    c4    c5    c7   c8a  c8bc 
#0.153 0.255 0.260 0.217 0.102 0.274 0.213 0.924 0.308 
#Main_rare_cil
round(table(smoke_down_withSubclusters$comparison[which(smoke_down_withSubclusters$main_rare_cil == "Unique")])/
	table(smoke_down_withSubclusters$comparison)[-c(10)],3)
#   c1   c1c    c2    c3    c4    c5   c6a   c6c    c7   c8a  c8bc 
#0.153 0.176 0.178 0.217 0.051 0.232 0.881 0.810 0.170 0.848 0.205 


######### Now tally down the number and proportion of DEGs that are core (this will be the same across all four datasets)
#Main
round(table(smoke_down_withSubclusters$comparison[which(smoke_down_withSubclusters$main_rare == "Core" | 
	smoke_down_withSubclusters$main_cil == "Core")]))
#  c1  c1c   c2   c3   c4   c5  c6a  c6c   c7   c8  c8a c8bc 
#  26   20   25   18   28   23    9    1   15   11    6   11
round(table(smoke_down_withSubclusters$comparison[which(smoke_down_withSubclusters$main_rare == "Core" | 
	smoke_down_withSubclusters$main_cil == "Core")]) / table(smoke_down_withSubclusters$comparison),3)
#   c1   c1c    c2    c3    c4    c5   c6a   c6c    c7    c8   c8a  c8bc 
#0.441 0.392 0.342 0.391 0.475 0.242 0.011 0.004 0.319 0.250 0.021 0.282










##############Finally, compare c8 DEGs to c8bc DEGs - does including the hybrids matter?
#89% (199 out of 224) of c8 (with hybrids) up DEGs are still there when hybrids are removed (c8bc DEGs) (only 51 novel DEGs with c8bc)
length(intersect(smoke_up_withSubclusters$gene[which(smoke_up_withSubclusters$comparison == "c8")],
	smoke_up_withSubclusters$gene[which(smoke_up_withSubclusters$comparison == "c8bc")]))/
	length(smoke_up_withSubclusters$gene[which(smoke_up_withSubclusters$comparison == "c8bc")])

#84% (37 out of 44) of c8 (with hybrids) down DEGs are still there when hybrids are removed (c8bc DEGs) (only two novel DEGs with c8bc)
length(intersect(smoke_down_withSubclusters$gene[which(smoke_down_withSubclusters$comparison == "c8")],
	smoke_down_withSubclusters$gene[which(smoke_down_withSubclusters$comparison == "c8bc")]))/
	length(smoke_down_withSubclusters$gene[which(smoke_down_withSubclusters$comparison == "c8")])











##############Finally, as we have decided to focus on the Main_cil set of populations, let's do enrichments on
##############just those unique genes and populations
#First, isolate the unique genes and the appropriate set of populations
smoke_up_mainCil_unique_forEnrichr<-smoke_up_withSubclusters[which(smoke_up_withSubclusters$main_cil == "Unique"),]
smoke_down_mainCil_unique_forEnrichr<-smoke_down_withSubclusters[which(smoke_down_withSubclusters$main_cil == "Unique"),]

#Write all to single txt file
write.table(smoke_up_mainCil_unique_forEnrichr,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,file="Enrichr/Smoke_up_mainCil_uniqueDEGs_forEnrichr.txt")
write.table(smoke_down_mainCil_unique_forEnrichr,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,file="Enrichr/Smoke_down_mainCil_uniqueDEGs_forEnrichr.txt")

#Do enrichments
DEG_table<-read.table("Enrichr/Smoke_up_mainCil_uniqueDEGs_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"Smoke_up_mainCil_uniqueDEGs"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)

DEG_table<-read.table("Enrichr/Smoke_down_mainCil_uniqueDEGs_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"Smoke_down_mainCil_uniqueDEGs"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)



##############And one more set for just the unique and semiunique together, excluding the core
#Bring in datasets if necessary
smoke_up_withSubclusters<-read.table("DEG_Tables/Smoking_up_DEGs_coreUniqueDesignations.txt",header=T,stringsAsFactors=F)
smoke_down_withSubclusters<-read.table("DEG_Tables/Smoking_down_DEGs_coreUniqueDesignations.txt",header=T,stringsAsFactors=F)

#First, isolate the unique genes and the appropriate set of populations
smoke_up_mainCil_semiunique_forEnrichr<-smoke_up_withSubclusters[which(smoke_up_withSubclusters$main_cil == "Unique" |
	smoke_up_withSubclusters$main_cil == "Semiunique"),]
smoke_down_mainCil_semiunique_forEnrichr<-smoke_down_withSubclusters[which(smoke_down_withSubclusters$main_cil == "Unique" |
	smoke_down_withSubclusters$main_cil == "Semiunique"),]

#Write all to single txt file
write.table(smoke_up_mainCil_semiunique_forEnrichr,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,file="Enrichr/smoke_up_mainCil_semiuniqueDEGs_forEnrichr.txt")
write.table(smoke_down_mainCil_semiunique_forEnrichr,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,file="Enrichr/smoke_down_mainCil_semiuniqueDEGs_forEnrichr.txt")

#Do enrichments
DEG_table<-read.table("Enrichr/smoke_up_mainCil_semiuniqueDEGs_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"smoke_up_mainCil_semiuniqueDEGs"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)

DEG_table<-read.table("Enrichr/smoke_down_mainCil_semiuniqueDEGs_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"smoke_down_mainCil_semiuniqueDEGs"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)




######## Finally export tables
write.table(smoke_up_withSubclusters,quote=F,row.names=F,file="DEG_Tables/Smoking_up_DEGs_coreUniqueDesignations.txt")
writeDEGsToExcel(smoke_up_withSubclusters,"DEG_Tables/Smoking_up_DEGs_coreUniqueDesignations.xlsx")

write.table(smoke_down_withSubclusters,quote=F,row.names=F,file="DEG_Tables/Smoking_down_DEGs_coreUniqueDesignations.txt")
writeDEGsToExcel(smoke_down_withSubclusters,"DEG_Tables/Smoking_down_DEGs_coreUniqueDesignations.xlsx")
































################################ Make pie charts that describe proportion of core, semiunique, and unique DEGs per cluster

#Proportion of smoking DEGs were unique or core
#First, bring in smoking DEGs for clusters and subclusters
genes_smokeUp<-read.table("DEG_Tables/Smoking_up_DEGs_coreUniqueDesignations.txt",header=T,stringsAsFactors=F)
genes_smokeDown<-read.table("DEG_Tables/Smoking_down_DEGs_coreUniqueDesignations.txt",header=T,stringsAsFactors=F)

#######################UP
#Get proportions
#UNIQUE
#Non-rare
propUnique<-round(c(table(genes_smokeUp$comparison[which(genes_smokeUp$main_cil == "Unique")])[-8],"c8a"=0,table(genes_smokeUp$comparison[which(genes_smokeUp$main_cil == "Unique")])[8]) / 
	table(genes_smokeUp$comparison[which(!is.na(genes_smokeUp$main_cil))]),2)
#  c1  c1c   c2   c3   c4   c5   c7  c8a c8bc 
#0.10 0.47 0.04 0.08 0.16 0.17 0.20 0.00 0.54
#Rare
propUnique_rare<-round(table(genes_smokeUp$comparison[which(genes_smokeUp$main_rare == "Unique")])[7:9] / 
	table(genes_smokeUp$comparison[which(!is.na(genes_smokeUp$main_rare))])[7:9],2)
# c6a  c6b  c6c 
#0.33 0.76 0.47 
#Combine non-rare and rare
propUnique<-append(propUnique,propUnique_rare)


#CORE
#Non-rare
propCore<-round(table(genes_smokeUp$comparison[which(genes_smokeUp$main_cil == "Core")]) / 
	table(genes_smokeUp$comparison[which(!is.na(genes_smokeUp$main_cil))]),2)
#  c1  c1c   c2   c3   c4   c5   c7  c8a c8bc 
#0.32 0.17 0.63 0.44 0.38 0.38 0.46 0.33 0.20 
#Rare
propCore_rare<-round(table(genes_smokeUp$comparison[which(genes_smokeUp$main_rare == "Core")])[7:9] / 
	table(genes_smokeUp$comparison[which(!is.na(genes_smokeUp$main_rare))])[7:9],2)
# c6a  c6b  c6c 
#0.31 0.04 0.18 
propCore<-append(propCore,propCore_rare)


#SEMI-UNIQUE
#Non-rare
propSemiUnique<-round(table(genes_smokeUp$comparison[which(genes_smokeUp$main_cil == "Semiunique")]) / 
	table(genes_smokeUp$comparison[which(!is.na(genes_smokeUp$main_cil))]),2)
#  c1  c1c   c2   c3   c4   c5   c7  c8a c8bc 
#0.58 0.37 0.32 0.48 0.46 0.45 0.34 0.67 0.26 
#Rare
propSemiUnique_rare<-round(table(genes_smokeUp$comparison[which(genes_smokeUp$main_rare == "Semiunique")])[7:9] / 
	table(genes_smokeUp$comparison[which(!is.na(genes_smokeUp$main_rare))])[7:9],2)
# c6a  c6b  c6c 
#0.37 0.20 0.35 
propSemiUnique<-append(propSemiUnique,propSemiUnique_rare)



#Get list of proportions
propList<-list()
for(i in 1:length(propSemiUnique)){
	propList[[length(propList)+1]]<-c(propCore[i],propUnique[i],propSemiUnique[i])
	names(propList)[[length(propList)]]<-names(propUnique)[i]
}

#Plot pie
pdf("T15_pieCharts_smokerUniqueCoreProportions_up.pdf")
par(mfrow=(c(3,3)))
for(i in 1:length(propList)){
	pct <- propList[[i]]*100
	lbls <- paste(pct,"%",sep="") # ad % to labels 
	colors<-c("darkolivegreen2","darkorchid4","thistle3")
	pie(propList[[i]],labels = lbls, col=colors, cex=2,main = names(propList)[[i]])
}
dev.off()







#######################Down
#Get proportions
#UNIQUE
#Non-rare
propUnique<-propCore<-round(table(genes_smokeDown$comparison[which(genes_smokeDown$main_cil == "Unique")]) / 
	table(genes_smokeDown$comparison[which(!is.na(genes_smokeDown$main_cil))]),2)
#  c1  c1c   c2   c3   c4   c5   c7  c8a c8bc 
#0.15 0.25 0.26 0.22 0.10 0.27 0.21 0.92 0.31 
#Rare
propUnique_rare<-round(table(genes_smokeDown$comparison[which(genes_smokeDown$main_rare == "Unique")])[7:8] / 
	table(genes_smokeDown$comparison[which(!is.na(genes_smokeDown$main_rare))])[7:8],2)
# c6a  c6c
#0.88 0.81
#Combine non-rare and rare
propUnique<-append(propUnique,propUnique_rare)


#CORE
#Non-rare
propCore<-round(table(genes_smokeDown$comparison[which(genes_smokeDown$main_cil == "Core")]) / 
	table(genes_smokeDown$comparison[which(!is.na(genes_smokeDown$main_cil))]),2)
#  c1  c1c   c2   c3   c4   c5   c7  c8a c8bc 
#0.44 0.39 0.34 0.39 0.47 0.24 0.32 0.02 0.28 
#Rare
propCore_rare<-round(table(genes_smokeDown$comparison[which(genes_smokeDown$main_rare == "Core")])[7:8] / 
	table(genes_smokeDown$comparison[which(!is.na(genes_smokeDown$main_rare))])[7:8],2)
# c6a  c6c 
#0.01 0.00 
propCore<-append(propCore,propCore_rare)


#SEMI-UNIQUE
#Non-rare
propSemiUnique<-round(table(genes_smokeDown$comparison[which(genes_smokeDown$main_cil == "Semiunique")]) / 
	table(genes_smokeDown$comparison[which(!is.na(genes_smokeDown$main_cil))]),2)
#  c1  c1c   c2   c3   c4   c5   c7  c8a c8bc 
#0.41 0.35 0.40 0.39 0.42 0.48 0.47 0.06 0.41 
#Rare
propSemiUnique_rare<-round(table(genes_smokeDown$comparison[which(genes_smokeDown$main_rare == "Semiunique")])[7:8] / 
	table(genes_smokeDown$comparison[which(!is.na(genes_smokeDown$main_rare))])[7:8],2)
# c6a  c6c 
#0.11 0.19 
propSemiUnique<-append(propSemiUnique,propSemiUnique_rare)



#Get list of proportions
propList<-list()
for(i in 1:length(propSemiUnique)){
	propList[[length(propList)+1]]<-c(propCore[i],propUnique[i],propSemiUnique[i])
	names(propList)[[length(propList)]]<-names(propUnique)[i]
}

#Plot pie
pdf("T15_pieCharts_smokerUniqueCoreProportions_down.pdf")
par(mfrow=(c(3,3)))
for(i in 1:length(propList)){
	pct <- propList[[i]]*100
	lbls <- paste(pct,"%",sep="") # ad % to labels 
	colors<-c("darkolivegreen2","darkorchid4","thistle3")
	pie(propList[[i]],labels = lbls, col=colors, cex=2,main=names(propList)[[i]])
}
dev.off()






























########################################## BEANE ET AL GENES OVERLAY
########### For the paper, do one more version with the beane up and down genes
dataSCT<-as.matrix(GetAssayData(T15_int,slot="data",assay="SCT"))
beane_up<-unique(scan("Gene_lists/Beane.etal_currSmokingDEGs_up.txt",what="character"))
beane_down<-unique(scan("Gene_lists/Beane.etal_currSmokingDEGs_down.txt",what="character"))

#Also, have object that merges these altogether in whatever way you want
smoking_genes<-unique(c(beane_up)) #up
#smoking_genes<-unique(c(beane_down)) #down

#Only those genes present in the dataset
smoking_genes<-smoking_genes[which(smoking_genes %in% rownames(dataSCT))]

##Add arithmetic mean of smoking genes to metadata
#smoking_genes_mean<-data.frame("smoking_genes"=colMeans(dataSCT[smoking_genes,]),row.names=rownames(T15_int@meta.data))
#T15_int <- AddMetaData(T15_int,smoking_genes_mean)

#Add geometric mean of smoking genes to metadata
smoking_genes_mean<-data.frame("smoking_genes"=exp(colMeans(log(dataSCT[smoking_genes,] + 1))) - 1,row.names=rownames(T15_int@meta.data))
T15_int <- AddMetaData(T15_int,smoking_genes_mean)

#Make dataframe with mean expression of smoking genes across each cluster and smoking status
smoking_genesDF<-data.frame()
for(i in sort(unique(T15_int@meta.data$clusters_10))){
	smoking_genesDF<-rbind(smoking_genesDF,data.frame("norm_exp"=c(T15_int@meta.data$smoking_genes[which(T15_int@meta.data$clusters_10 == i & 
		T15_int@meta.data$smoke == "never")],
		T15_int@meta.data$smoking_genes[which(T15_int@meta.data$clusters_10 == i & 
		T15_int@meta.data$smoke == "heavy")]),
		"group"=c(rep(paste(i,"_nonsmoker",sep=""),length(T15_int@meta.data$smoking_genes[which(T15_int@meta.data$clusters_10 == i & 
		T15_int@meta.data$smoke == "never")])),
		rep(paste(i,"_smoker",sep=""),length(T15_int@meta.data$smoking_genes[which(T15_int@meta.data$clusters_10 == i & 
		T15_int@meta.data$smoke == "heavy")])))))
}

#Test differences between heavy and non smokers (for UP)
pvals<-data.frame("pvals"=rep(NA,length(sort(unique(T15_int@meta.data$clusters_10)))),"cluster"=sort(unique(T15_int@meta.data$clusters_10)))
count<-1
for(i in seq(1,length(levels(smoking_genesDF$group)),by=2)){
	pvals$pvals[count]<-t.test(smoking_genesDF$norm_exp[which(smoking_genesDF$group == levels(smoking_genesDF$group)[i + 1])],
		smoking_genesDF$norm_exp[which(smoking_genesDF$group == levels(smoking_genesDF$group)[i])],
		alternative="greater")$p.value
	count<-count + 1
}

#Test differences between heavy and non smokers (for DOWN)
pvals<-data.frame("pvals"=rep(NA,length(sort(unique(T15_int@meta.data$clusters_10)))),"cluster"=sort(unique(T15_int@meta.data$clusters_10)))
count<-1
for(i in seq(1,length(levels(smoking_genesDF$group)),by=2)){
	pvals$pvals[count]<-t.test(smoking_genesDF$norm_exp[which(smoking_genesDF$group == levels(smoking_genesDF$group)[i + 1])],
		smoking_genesDF$norm_exp[which(smoking_genesDF$group == levels(smoking_genesDF$group)[i])],
		alternative="less")$p.value
	count<-count + 1
}
	
#Make box plots
#pdf("SmokingGenes_beane_up_boxplots.pdf")`
#Get colorVec based on pvalues
colorVec<-rep("black",length(sort(unique(T15_int@meta.data$clusters_10))))
colorVec[which(pvals$pvals <= 0.05)]<-"red"
dev.new(height=4,width=6)
par(bty="l")
boxplot(smoking_genesDF$norm_exp~smoking_genesDF$group,las=1,col=c("black","red"),pch=16,cex=0.5,ylab="Smoking genes (activated)",
	names=rep("",2*10),medcol="white",outcex=0.3,medlwd=1,lwd=0.8)
text(x=seq(1.5,20.5,by=2),y=par()$usr[3]-0.12*(par()$usr[4]-par()$usr[3]),labels=sort(unique(T15_int@meta.data$clusters_10)),
	srt=45,adj=1,xpd=T,cex=0.8,col=colorVec)








####### Do Fisher test of enrichment for each cluster. Are DEGs enriched for Beane genes?
########### Do the test
#Bring in the significant up smoking DEGs
T15_clusters10_smoke_up_DEGs_forEnrichr<-read.table("Enrichr/T15_clusters10_smoke_down_DEGs_forEnrichr.txt",header=T,stringsAsFactors=F)

#Bring in the Beane genes
beane_up<-unique(scan("Gene_lists/Beane.etal_currSmokingDEGs_down.txt",what="character"))

#Only those genes present in the dataset
smoking_genes<-smoking_genes[which(smoking_genes %in% rownames(dataSCT))]

#Bring in data
genes<-smoking_genes

#Bring in modules
mods<-data.frame("Gene"=T15_clusters10_smoke_up_DEGs_forEnrichr$gene,"Module"=T15_clusters10_smoke_up_DEGs_forEnrichr$comparison,stringsAsFactors=F)

#Do enrichments
upEnrich<-DoEnrichment(genes=genes,gene2module=mods,background=nrow(GetAssayData(T15_int,assay="SCT")))








#######Also, look at the proportion of up Beane genes that fall within DEGs of the 24 populations

#Get proportions of detected Beane genes that are activated in each of the cell populations
proportionVec<-data.frame("cluster"=sort(unique(T15_int@meta.data$clusters_10)),
	"proportion"=rep(NA,length(sort(unique(T15_int@meta.data$clusters_10)))),
	"Num_beane"=rep(NA,length(sort(unique(T15_int@meta.data$clusters_10)))),
	"Num_degs"=rep(NA,length(sort(unique(T15_int@meta.data$clusters_10)))))
count<-1
for(i in sort(unique(T15_int@meta.data$clusters_10))){
	proportionVec$proportion[count]<-round(length(which(smoking_genes %in% 
		T15_clusters10_smoke_up_DEGs_forEnrichr$gene[which(T15_clusters10_smoke_up_DEGs_forEnrichr$comparison == i)])) / 
		length(smoking_genes),2)
	proportionVec$Num_beane[count]<-length(which(smoking_genes %in% 
		T15_clusters10_smoke_up_DEGs_forEnrichr$gene[which(T15_clusters10_smoke_up_DEGs_forEnrichr$comparison == i)]))
	proportionVec$Num_degs[count]<-length(T15_clusters10_smoke_up_DEGs_forEnrichr$gene[which(T15_clusters10_smoke_up_DEGs_forEnrichr$comparison == i)])
	count<-count + 1
}

#Add in Fisher test padjust to table
proportionVec<-cbind(proportionVec,"padj"=upEnrich$padj)

#Print proportion output table
write.table(proportionVec,row.names=F,sep="\t",quote=F,file="OverlapOfBeaneGenesInClusterDEGs.txt")

#Color points based on whether they are significantly enriched (FDR < 0.05)
colorVec<-rep("black",length(sort(unique(T15_int@meta.data$clusters_10))))
colorVec[which(upEnrich$padj <= 0.05)]<-"red"

#Make line plots of proportions
#pdf("LinePlot_overlapOfBeaneGenesInClusters.pdf")
dev.new(width=3.5,height=3)
barplot(proportionVec$proportion, col = NA, border=NA,las=1,bty="l",ylim=c(0,0.1),ylab="Proportion of Beane et al. smoking genes")
points(proportionVec$proportion,pch=16,col=colorVec)
lines(proportionVec$proportion)
text(x=seq(1,10),y=par()$usr[3]-0.12*(par()$usr[4]-par()$usr[3]),labels=sort(unique(T15_int@meta.data$clusters_10)),
	srt=45,adj=1,xpd=T,cex=0.9)
	
	
	
	
	


























#################################### SMOKING EFFECTS - TO WHAT EXTENT ARE THEY ENRICHED FOR BASELINE CELL TYPE SIGNATURES
############ Here, we want to ask of core, total, and unique smoking signatures, to what extent are these similar to 
############ cell type signatures, as estimated from nonsmoker cells only?

#First, bring in smoking DEGs for clusters and subclusters
genes_smokeUp<-read.table("DEG_Tables/Smoking_up_DEGs_coreUniqueDesignations.txt",header=T,stringsAsFactors=F)
genes_smokeDown<-read.table("DEG_Tables/Smoking_down_DEGs_coreUniqueDesignations.txt",header=T,stringsAsFactors=F)

#Also, break these into core, non-core, and unique (for the core, there's just one set)
#CORE
genes_smokeUp_core<-unique(genes_smokeUp[which(genes_smokeUp$main_cil == "Core"),]$gene)
genes_smokeDown_core<-unique(genes_smokeDown[which(genes_smokeDown$main_cil == "Core"),]$gene)
#NONCORE
genes_smokeUp_noncore<-genes_smokeUp[-which(genes_smokeUp$main_cil == "Core"),]
genes_smokeDown_noncore<-genes_smokeDown[-which(genes_smokeDown$main_cil == "Core"),]
#UNIQUE
genes_smokeUp_unique<-genes_smokeUp[which(genes_smokeUp$main_cil == "Unique"),]
genes_smokeDown_unique<-genes_smokeDown[which(genes_smokeDown$main_cil == "Unique"),]



#Next, bring in baseline cluster DEGs
genes_cellTypes<-read.table("DEG_Tables/T15_1againstAll_nonsmokers_DEGs.txt",header=T,stringsAsFactors=F)
genes_cellTypes<-genes_cellTypes[which(genes_cellTypes$p_adj_FDR < 0.05),]



#Finally, calculate the background for these tests
background=nrow(GetAssayData(T15_int,assay="SCT"))


#########Now run enrichments
#### ALL UP SMOKING DEGs
#Bring in smoking response "modules"
mods<-data.frame("Gene"=genes_smokeUp$gene,"Module"=genes_smokeUp$comparison,stringsAsFactors=F)
#The cycle enrichments through each of the clusters
enrichments_up_allSmokeDEGs<-list()
for(i in 1:length(unique(genes_cellTypes$comparison))){
	genes<-genes_cellTypes[which(genes_cellTypes$comparison == unique(genes_cellTypes$comparison)[i]),]$gene
	enrichments_up_allSmokeDEGs[[length(enrichments_up_allSmokeDEGs) + 1]]<-DoEnrichment(genes=genes,gene2module=mods,background=background)
	colnames(enrichments_up_allSmokeDEGs[[length(enrichments_up_allSmokeDEGs)]])[1]<-"smoke_DEGs"
	names(enrichments_up_allSmokeDEGs)[[length(enrichments_up_allSmokeDEGs)]]<-unique(genes_cellTypes$comparison)[i]
	#Also, add in the number of smoking DEGs
	Nsmoke<-as.vector(table(mods$Module))
	#And the number of cell type signature DEGs
	NcellType<-length(genes)
	#And the number of smoking DEGs that overlap with the cell type signature DEGs
	overlap<-unlist(lapply(split(mods,f=mods$Module),function(x)sum(genes %in% x$Gene)))
	#And the proportion of smoking DEGs that overlap with the given cell type signature (get rid of this for now)
	#prop_smoke<-round(unlist(lapply(split(mods,f=mods$Module),function(x)sum(x$Gene %in% genes) / length(x$Gene))),2)
	#And the proportion of cell type signature genes that fall into the list of smoking DEGs
	prop_cellType<-round(unlist(lapply(split(mods,f=mods$Module),function(x)sum(x$Gene %in% genes) / length(genes))),2)
	#Merge these into a single table
	enrichments_up_allSmokeDEGs[[length(enrichments_up_allSmokeDEGs)]]<-cbind(enrichments_up_allSmokeDEGs[[length(enrichments_up_allSmokeDEGs)]],
		"N_smoke_DEGs"=Nsmoke,"N_cluster_DEGs"=NcellType,"N_overlap"=overlap,"Prop_overlap_cluster"=prop_cellType)
	#Finally, order by FDR
	enrichments_up_allSmokeDEGs[[length(enrichments_up_allSmokeDEGs)]]<-enrichments_up_allSmokeDEGs[[length(enrichments_up_allSmokeDEGs)]][order(
		enrichments_up_allSmokeDEGs[[length(enrichments_up_allSmokeDEGs)]]$padj),]
}

#Write to excel
#write.xlsx(enrichments_up_allSmokeDEGs,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_upAll_byClusterDEGs.xlsx")

#Also print these in a single table so they can be sorted by smoke_DEGs if desired
enrichments_up_allSmokeDEGs_oneTab<-cbind(enrichments_up_allSmokeDEGs[[1]],"cluster_DEGs"=names(enrichments_up_allSmokeDEGs)[[1]])[,c(1,8,2:7)]
for(i in 2:length(enrichments_up_allSmokeDEGs)){
	enrichments_up_allSmokeDEGs_oneTab<-rbind(enrichments_up_allSmokeDEGs_oneTab,cbind(enrichments_up_allSmokeDEGs[[i]],
		"cluster_DEGs"=names(enrichments_up_allSmokeDEGs)[[i]])[,c(1,8,2:7)])
}
#write.xlsx(enrichments_up_allSmokeDEGs_oneTab,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_upAll_oneTab.xlsx")

#Finally, flip, so there is one tab for each smoke_DEG cluster (this makes the most sense to me)
enrichments_up_allSmokeDEGs_bySmokeDEGs<-split(enrichments_up_allSmokeDEGs_oneTab,f=enrichments_up_allSmokeDEGs_oneTab$smoke_DEGs)
enrichments_up_allSmokeDEGs_bySmokeDEGs<-lapply(enrichments_up_allSmokeDEGs_bySmokeDEGs,function(x)x[order(x$padj),])
write.xlsx(enrichments_up_allSmokeDEGs_bySmokeDEGs,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_upAll.xlsx")










#### ALL DOWN SMOKING DEGs
#Bring in smoking response "modules"
mods<-data.frame("Gene"=genes_smokeDown$gene,"Module"=genes_smokeDown$comparison,stringsAsFactors=F)
#The cycle enrichments through each of the clusters
enrichments_down_allSmokeDEGs<-list()
for(i in 1:length(unique(genes_cellTypes$comparison))){
	genes<-genes_cellTypes[which(genes_cellTypes$comparison == unique(genes_cellTypes$comparison)[i]),]$gene
	enrichments_down_allSmokeDEGs[[length(enrichments_down_allSmokeDEGs) + 1]]<-DoEnrichment(genes=genes,gene2module=mods,background=background)
	colnames(enrichments_down_allSmokeDEGs[[length(enrichments_down_allSmokeDEGs)]])[1]<-"smoke_DEGs"
	names(enrichments_down_allSmokeDEGs)[[length(enrichments_down_allSmokeDEGs)]]<-unique(genes_cellTypes$comparison)[i]
	#Also, add in the number of smoking DEGs
	Nsmoke<-as.vector(table(mods$Module))
	#And the number of cell type signature DEGs
	NcellType<-length(genes)
	#And the number of smoking DEGs that overlap with the cell type signature DEGs
	overlap<-unlist(lapply(split(mods,f=mods$Module),function(x)sum(genes %in% x$Gene)))
	#And the proportion of smoking DEGs that overlap with the given cell type signature (get rid of this for now)
	#prop_smoke<-round(unlist(lapply(split(mods,f=mods$Module),function(x)sum(x$Gene %in% genes) / length(x$Gene))),2)
	#And the proportion of cell type signature genes that fall into the list of smoking DEGs
	prop_cellType<-round(unlist(lapply(split(mods,f=mods$Module),function(x)sum(x$Gene %in% genes) / length(genes))),2)
	#Merge these into a single table
	enrichments_down_allSmokeDEGs[[length(enrichments_down_allSmokeDEGs)]]<-cbind(enrichments_down_allSmokeDEGs[[length(enrichments_down_allSmokeDEGs)]],
		"N_smoke_DEGs"=Nsmoke,"N_cluster_DEGs"=NcellType,"N_overlap"=overlap,"Prop_overlap_cluster"=prop_cellType)
	#Finally, order by FDR
	enrichments_down_allSmokeDEGs[[length(enrichments_down_allSmokeDEGs)]]<-enrichments_down_allSmokeDEGs[[length(enrichments_down_allSmokeDEGs)]][order(
		enrichments_down_allSmokeDEGs[[length(enrichments_down_allSmokeDEGs)]]$padj),]
}

#Write to excel
#write.xlsx(enrichments_down_allSmokeDEGs,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_downAll_byClusterDEGs.xlsx")

#Also print these in a single table so they can be sorted by smoke_DEGs if desired
enrichments_down_allSmokeDEGs_oneTab<-cbind(enrichments_down_allSmokeDEGs[[1]],"cluster_DEGs"=names(enrichments_down_allSmokeDEGs)[[1]])[,c(1,8,2:7)]
for(i in 2:length(enrichments_down_allSmokeDEGs)){
	enrichments_down_allSmokeDEGs_oneTab<-rbind(enrichments_down_allSmokeDEGs_oneTab,cbind(enrichments_down_allSmokeDEGs[[i]],
		"cluster_DEGs"=names(enrichments_down_allSmokeDEGs)[[i]])[,c(1,8,2:7)])
}
#write.xlsx(enrichments_down_allSmokeDEGs_oneTab,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_downAll_oneTab.xlsx")

#Finally, flip, so there is one tab for each smoke_DEG cluster (this makes the most sense to me)
enrichments_down_allSmokeDEGs_bySmokeDEGs<-split(enrichments_down_allSmokeDEGs_oneTab,f=enrichments_down_allSmokeDEGs_oneTab$smoke_DEGs)
enrichments_down_allSmokeDEGs_bySmokeDEGs<-lapply(enrichments_down_allSmokeDEGs_bySmokeDEGs,function(x)x[order(x$padj),])
write.xlsx(enrichments_down_allSmokeDEGs_bySmokeDEGs,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_downAll.xlsx")














#### UNIQUE UP SMOKING DEGs
#Bring in smoking response "modules"
mods<-data.frame("Gene"=genes_smokeUp_unique$gene,"Module"=genes_smokeUp_unique$comparison,stringsAsFactors=F)
#The cycle enrichments through each of the clusters
enrichments_up_uniqueSmokeDEGs<-list()
for(i in 1:length(unique(genes_cellTypes$comparison))){
	genes<-genes_cellTypes[which(genes_cellTypes$comparison == unique(genes_cellTypes$comparison)[i]),]$gene
	enrichments_up_uniqueSmokeDEGs[[length(enrichments_up_uniqueSmokeDEGs) + 1]]<-DoEnrichment(genes=genes,gene2module=mods,background=background)
	colnames(enrichments_up_uniqueSmokeDEGs[[length(enrichments_up_uniqueSmokeDEGs)]])[1]<-"smoke_DEGs"
	names(enrichments_up_uniqueSmokeDEGs)[[length(enrichments_up_uniqueSmokeDEGs)]]<-unique(genes_cellTypes$comparison)[i]
	#Also, add in the number of smoking DEGs
	Nsmoke<-as.vector(table(mods$Module))
	#And the number of cell type signature DEGs
	NcellType<-length(genes)
	#And the number of smoking DEGs that overlap with the cell type signature DEGs
	overlap<-unlist(lapply(split(mods,f=mods$Module),function(x)sum(genes %in% x$Gene)))
	#And the proportion of smoking DEGs that overlap with the given cell type signature (get rid of this for now)
	#prop_smoke<-round(unlist(lapply(split(mods,f=mods$Module),function(x)sum(x$Gene %in% genes) / length(x$Gene))),2)
	#And the proportion of cell type signature genes that funique into the list of smoking DEGs
	prop_cellType<-round(unlist(lapply(split(mods,f=mods$Module),function(x)sum(x$Gene %in% genes) / length(genes))),2)
	#Merge these into a single table
	enrichments_up_uniqueSmokeDEGs[[length(enrichments_up_uniqueSmokeDEGs)]]<-cbind(enrichments_up_uniqueSmokeDEGs[[length(enrichments_up_uniqueSmokeDEGs)]],
		"N_smoke_DEGs"=Nsmoke,"N_cluster_DEGs"=NcellType,"N_overlap"=overlap,"Prop_overlap_cluster"=prop_cellType)
	#Finuniquey, order by FDR
	enrichments_up_uniqueSmokeDEGs[[length(enrichments_up_uniqueSmokeDEGs)]]<-enrichments_up_uniqueSmokeDEGs[[length(enrichments_up_uniqueSmokeDEGs)]][order(
		enrichments_up_uniqueSmokeDEGs[[length(enrichments_up_uniqueSmokeDEGs)]]$padj),]
}

#Write to excel
#write.xlsx(enrichments_up_uniqueSmokeDEGs,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_upUnique_byClusterDEGs.xlsx")

#Also print these in a single table so they can be sorted by smoke_DEGs if desired
enrichments_up_uniqueSmokeDEGs_oneTab<-cbind(enrichments_up_uniqueSmokeDEGs[[1]],"cluster_DEGs"=names(enrichments_up_uniqueSmokeDEGs)[[1]])[,c(1,8,2:7)]
for(i in 2:length(enrichments_up_uniqueSmokeDEGs)){
	enrichments_up_uniqueSmokeDEGs_oneTab<-rbind(enrichments_up_uniqueSmokeDEGs_oneTab,cbind(enrichments_up_uniqueSmokeDEGs[[i]],
		"cluster_DEGs"=names(enrichments_up_uniqueSmokeDEGs)[[i]])[,c(1,8,2:7)])
}
#write.xlsx(enrichments_up_uniqueSmokeDEGs_oneTab,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_upUnique_oneTab.xlsx")

#Finally, flip, so there is one tab for each smoke_DEG cluster (this makes the most sense to me)
enrichments_up_uniqueSmokeDEGs_bySmokeDEGs<-split(enrichments_up_uniqueSmokeDEGs_oneTab,f=enrichments_up_uniqueSmokeDEGs_oneTab$smoke_DEGs)
enrichments_up_uniqueSmokeDEGs_bySmokeDEGs<-lapply(enrichments_up_uniqueSmokeDEGs_bySmokeDEGs,function(x)x[order(x$padj),])
write.xlsx(enrichments_up_uniqueSmokeDEGs_bySmokeDEGs,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_upUnique.xlsx")














#### UNIQUE DOWN SMOKING DEGs
#Bring in smoking response "modules"
mods<-data.frame("Gene"=genes_smokeDown_unique$gene,"Module"=genes_smokeDown_unique$comparison,stringsAsFactors=F)
#The cycle enrichments through each of the clusters
enrichments_down_uniqueSmokeDEGs<-list()
for(i in 1:length(unique(genes_cellTypes$comparison))){
	genes<-genes_cellTypes[which(genes_cellTypes$comparison == unique(genes_cellTypes$comparison)[i]),]$gene
	enrichments_down_uniqueSmokeDEGs[[length(enrichments_down_uniqueSmokeDEGs) + 1]]<-DoEnrichment(genes=genes,gene2module=mods,background=background)
	colnames(enrichments_down_uniqueSmokeDEGs[[length(enrichments_down_uniqueSmokeDEGs)]])[1]<-"smoke_DEGs"
	names(enrichments_down_uniqueSmokeDEGs)[[length(enrichments_down_uniqueSmokeDEGs)]]<-unique(genes_cellTypes$comparison)[i]
	#Also, add in the number of smoking DEGs
	Nsmoke<-as.vector(table(mods$Module))
	#And the number of cell type signature DEGs
	NcellType<-length(genes)
	#And the number of smoking DEGs that overlap with the cell type signature DEGs
	overlap<-unlist(lapply(split(mods,f=mods$Module),function(x)sum(genes %in% x$Gene)))
	#And the proportion of smoking DEGs that overlap with the given cell type signature (get rid of this for now)
	#prop_smoke<-round(unlist(lapply(split(mods,f=mods$Module),function(x)sum(x$Gene %in% genes) / length(x$Gene))),2)
	#And the proportion of cell type signature genes that funique into the list of smoking DEGs
	prop_cellType<-round(unlist(lapply(split(mods,f=mods$Module),function(x)sum(x$Gene %in% genes) / length(genes))),2)
	#Merge these into a single table
	enrichments_down_uniqueSmokeDEGs[[length(enrichments_down_uniqueSmokeDEGs)]]<-cbind(enrichments_down_uniqueSmokeDEGs[[length(enrichments_down_uniqueSmokeDEGs)]],
		"N_smoke_DEGs"=Nsmoke,"N_cluster_DEGs"=NcellType,"N_overlap"=overlap,"Prop_overlap_cluster"=prop_cellType)
	#Finuniquey, order by FDR
	enrichments_down_uniqueSmokeDEGs[[length(enrichments_down_uniqueSmokeDEGs)]]<-enrichments_down_uniqueSmokeDEGs[[length(enrichments_down_uniqueSmokeDEGs)]][order(
		enrichments_down_uniqueSmokeDEGs[[length(enrichments_down_uniqueSmokeDEGs)]]$padj),]
}

#Write to excel
#write.xlsx(enrichments_down_uniqueSmokeDEGs,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_downUnique_byClusterDEGs.xlsx")

#Also print these in a single table so they can be sorted by smoke_DEGs if desired
enrichments_down_uniqueSmokeDEGs_oneTab<-cbind(enrichments_down_uniqueSmokeDEGs[[1]],"cluster_DEGs"=names(enrichments_down_uniqueSmokeDEGs)[[1]])[,c(1,8,2:7)]
for(i in 2:length(enrichments_down_uniqueSmokeDEGs)){
	enrichments_down_uniqueSmokeDEGs_oneTab<-rbind(enrichments_down_uniqueSmokeDEGs_oneTab,cbind(enrichments_down_uniqueSmokeDEGs[[i]],
		"cluster_DEGs"=names(enrichments_down_uniqueSmokeDEGs)[[i]])[,c(1,8,2:7)])
}
#write.xlsx(enrichments_down_uniqueSmokeDEGs_oneTab,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_downUnique_oneTab.xlsx")

#Finally, flip, so there is one tab for each smoke_DEG cluster (this makes the most sense to me)
enrichments_down_uniqueSmokeDEGs_bySmokeDEGs<-split(enrichments_down_uniqueSmokeDEGs_oneTab,f=enrichments_down_uniqueSmokeDEGs_oneTab$smoke_DEGs)
enrichments_down_uniqueSmokeDEGs_bySmokeDEGs<-lapply(enrichments_down_uniqueSmokeDEGs_bySmokeDEGs,function(x)x[order(x$padj),])
write.xlsx(enrichments_down_uniqueSmokeDEGs_bySmokeDEGs,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_downUnique.xlsx")














#### CORE UP SMOKING DEGs
#Bring in smoking response "modules"
mods<-data.frame("Gene"=genes_smokeUp_core,"Module"="core",stringsAsFactors=F)
#The cycle enrichments through each of the clusters
enrichments_up_coreSmokeDEGs<-list()
for(i in 1:length(unique(genes_cellTypes$comparison))){
	genes<-genes_cellTypes[which(genes_cellTypes$comparison == unique(genes_cellTypes$comparison)[i]),]$gene
	enrichments_up_coreSmokeDEGs[[length(enrichments_up_coreSmokeDEGs) + 1]]<-DoEnrichment(genes=genes,gene2module=mods,background=background)
	colnames(enrichments_up_coreSmokeDEGs[[length(enrichments_up_coreSmokeDEGs)]])[1]<-"smoke_DEGs"
	names(enrichments_up_coreSmokeDEGs)[[length(enrichments_up_coreSmokeDEGs)]]<-unique(genes_cellTypes$comparison)[i]
	#Also, add in the number of smoking DEGs
	Nsmoke<-as.vector(table(mods$Module))
	#And the number of cell type signature DEGs
	NcellType<-length(genes)
	#And the number of smoking DEGs that overlap with the cell type signature DEGs
	overlap<-unlist(lapply(split(mods,f=mods$Module),function(x)sum(genes %in% x$Gene)))
	#And the proportion of smoking DEGs that overlap with the given cell type signature (get rid of this for now)
	#prop_smoke<-round(unlist(lapply(split(mods,f=mods$Module),function(x)sum(x$Gene %in% genes) / length(x$Gene))),2)
	#And the proportion of cell type signature genes that fcore into the list of smoking DEGs
	prop_cellType<-round(unlist(lapply(split(mods,f=mods$Module),function(x)sum(x$Gene %in% genes) / length(genes))),2)
	#Merge these into a single table
	enrichments_up_coreSmokeDEGs[[length(enrichments_up_coreSmokeDEGs)]]<-cbind(enrichments_up_coreSmokeDEGs[[length(enrichments_up_coreSmokeDEGs)]],
		"N_smoke_DEGs"=Nsmoke,"N_cluster_DEGs"=NcellType,"N_overlap"=overlap,"Prop_overlap_cluster"=prop_cellType)
	#Fincorey, order by FDR
	enrichments_up_coreSmokeDEGs[[length(enrichments_up_coreSmokeDEGs)]]<-enrichments_up_coreSmokeDEGs[[length(enrichments_up_coreSmokeDEGs)]][order(
		enrichments_up_coreSmokeDEGs[[length(enrichments_up_coreSmokeDEGs)]]$padj),]
}

#Write to excel
#write.xlsx(enrichments_up_coreSmokeDEGs,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_upCore_byClusterDEGs.xlsx")

#Also print these in a single table so they can be sorted by smoke_DEGs if desired
enrichments_up_coreSmokeDEGs_oneTab<-cbind(enrichments_up_coreSmokeDEGs[[1]],"cluster_DEGs"=names(enrichments_up_coreSmokeDEGs)[[1]])[,c(1,8,2:7)]
for(i in 2:length(enrichments_up_coreSmokeDEGs)){
	enrichments_up_coreSmokeDEGs_oneTab<-rbind(enrichments_up_coreSmokeDEGs_oneTab,cbind(enrichments_up_coreSmokeDEGs[[i]],
		"cluster_DEGs"=names(enrichments_up_coreSmokeDEGs)[[i]])[,c(1,8,2:7)])
}
#write.xlsx(enrichments_up_coreSmokeDEGs_oneTab,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_upCore_oneTab.xlsx")

#Finally, flip, so there is one tab for each smoke_DEG cluster (this makes the most sense to me)
enrichments_up_coreSmokeDEGs_bySmokeDEGs<-split(enrichments_up_coreSmokeDEGs_oneTab,f=enrichments_up_coreSmokeDEGs_oneTab$smoke_DEGs)
enrichments_up_coreSmokeDEGs_bySmokeDEGs<-lapply(enrichments_up_coreSmokeDEGs_bySmokeDEGs,function(x)x[order(x$padj),])
write.xlsx(enrichments_up_coreSmokeDEGs_bySmokeDEGs,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_upCore.xlsx")















#### CORE DOWN SMOKING DEGs
#Bring in smoking response "modules"
mods<-data.frame("Gene"=genes_smokeDown_core,"Module"="core",stringsAsFactors=F)
#The cycle enrichments through each of the clusters
enrichments_down_coreSmokeDEGs<-list()
for(i in 1:length(unique(genes_cellTypes$comparison))){
	genes<-genes_cellTypes[which(genes_cellTypes$comparison == unique(genes_cellTypes$comparison)[i]),]$gene
	enrichments_down_coreSmokeDEGs[[length(enrichments_down_coreSmokeDEGs) + 1]]<-DoEnrichment(genes=genes,gene2module=mods,background=background)
	colnames(enrichments_down_coreSmokeDEGs[[length(enrichments_down_coreSmokeDEGs)]])[1]<-"smoke_DEGs"
	names(enrichments_down_coreSmokeDEGs)[[length(enrichments_down_coreSmokeDEGs)]]<-unique(genes_cellTypes$comparison)[i]
	#Also, add in the number of smoking DEGs
	Nsmoke<-as.vector(table(mods$Module))
	#And the number of cell type signature DEGs
	NcellType<-length(genes)
	#And the number of smoking DEGs that overlap with the cell type signature DEGs
	overlap<-unlist(lapply(split(mods,f=mods$Module),function(x)sum(genes %in% x$Gene)))
	#And the proportion of smoking DEGs that overlap with the given cell type signature (get rid of this for now)
	#prop_smoke<-round(unlist(lapply(split(mods,f=mods$Module),function(x)sum(x$Gene %in% genes) / length(x$Gene))),2)
	#And the proportion of cell type signature genes that fcore into the list of smoking DEGs
	prop_cellType<-round(unlist(lapply(split(mods,f=mods$Module),function(x)sum(x$Gene %in% genes) / length(genes))),2)
	#Merge these into a single table
	enrichments_down_coreSmokeDEGs[[length(enrichments_down_coreSmokeDEGs)]]<-cbind(enrichments_down_coreSmokeDEGs[[length(enrichments_down_coreSmokeDEGs)]],
		"N_smoke_DEGs"=Nsmoke,"N_cluster_DEGs"=NcellType,"N_overlap"=overlap,"Prop_overlap_cluster"=prop_cellType)
	#Fincorey, order by FDR
	enrichments_down_coreSmokeDEGs[[length(enrichments_down_coreSmokeDEGs)]]<-enrichments_down_coreSmokeDEGs[[length(enrichments_down_coreSmokeDEGs)]][order(
		enrichments_down_coreSmokeDEGs[[length(enrichments_down_coreSmokeDEGs)]]$padj),]
}

#Write to excel
#write.xlsx(enrichments_down_coreSmokeDEGs,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_downCore_byClusterDEGs.xlsx")

#Also print these in a single table so they can be sorted by smoke_DEGs if desired
enrichments_down_coreSmokeDEGs_oneTab<-cbind(enrichments_down_coreSmokeDEGs[[1]],"cluster_DEGs"=names(enrichments_down_coreSmokeDEGs)[[1]])[,c(1,8,2:7)]
for(i in 2:length(enrichments_down_coreSmokeDEGs)){
	enrichments_down_coreSmokeDEGs_oneTab<-rbind(enrichments_down_coreSmokeDEGs_oneTab,cbind(enrichments_down_coreSmokeDEGs[[i]],
		"cluster_DEGs"=names(enrichments_down_coreSmokeDEGs)[[i]])[,c(1,8,2:7)])
}
#write.xlsx(enrichments_down_coreSmokeDEGs_oneTab,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_downCore_oneTab.xlsx")

#Finally, flip, so there is one tab for each smoke_DEG cluster (this makes the most sense to me)
enrichments_down_coreSmokeDEGs_bySmokeDEGs<-split(enrichments_down_coreSmokeDEGs_oneTab,f=enrichments_down_coreSmokeDEGs_oneTab$smoke_DEGs)
enrichments_down_coreSmokeDEGs_bySmokeDEGs<-lapply(enrichments_down_coreSmokeDEGs_bySmokeDEGs,function(x)x[order(x$padj),])
write.xlsx(enrichments_down_coreSmokeDEGs_bySmokeDEGs,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_downCore.xlsx")


















#### Non-core UP SMOKING DEGs
#Bring in smoking response "modules"
mods<-data.frame("Gene"=genes_smokeUp_noncore$gene,"Module"=genes_smokeUp_noncore$comparison,stringsAsFactors=F)
#The cycle enrichments through each of the clusters
enrichments_up_noncoreSmokeDEGs<-list()
for(i in 1:length(unique(genes_cellTypes$comparison))){
	genes<-genes_cellTypes[which(genes_cellTypes$comparison == unique(genes_cellTypes$comparison)[i]),]$gene
	enrichments_up_noncoreSmokeDEGs[[length(enrichments_up_noncoreSmokeDEGs) + 1]]<-DoEnrichment(genes=genes,gene2module=mods,background=background)
	colnames(enrichments_up_noncoreSmokeDEGs[[length(enrichments_up_noncoreSmokeDEGs)]])[1]<-"smoke_DEGs"
	names(enrichments_up_noncoreSmokeDEGs)[[length(enrichments_up_noncoreSmokeDEGs)]]<-unique(genes_cellTypes$comparison)[i]
	#Also, add in the number of smoking DEGs
	Nsmoke<-as.vector(table(mods$Module))
	#And the number of cell type signature DEGs
	NcellType<-length(genes)
	#And the number of smoking DEGs that overlap with the cell type signature DEGs
	overlap<-unlist(lapply(split(mods,f=mods$Module),function(x)sum(genes %in% x$Gene)))
	#And the proportion of smoking DEGs that overlap with the given cell type signature (get rid of this for now)
	#prop_smoke<-round(unlist(lapply(split(mods,f=mods$Module),function(x)sum(x$Gene %in% genes) / length(x$Gene))),2)
	#And the proportion of cell type signature genes that fall into the list of smoking DEGs
	prop_cellType<-round(unlist(lapply(split(mods,f=mods$Module),function(x)sum(x$Gene %in% genes) / length(genes))),2)
	#Merge these into a single table
	enrichments_up_noncoreSmokeDEGs[[length(enrichments_up_noncoreSmokeDEGs)]]<-cbind(enrichments_up_noncoreSmokeDEGs[[length(enrichments_up_noncoreSmokeDEGs)]],
		"N_smoke_DEGs"=Nsmoke,"N_cluster_DEGs"=NcellType,"N_overlap"=overlap,"Prop_overlap_cluster"=prop_cellType)
	#Finally, order by FDR
	enrichments_up_noncoreSmokeDEGs[[length(enrichments_up_noncoreSmokeDEGs)]]<-enrichments_up_noncoreSmokeDEGs[[length(enrichments_up_noncoreSmokeDEGs)]][order(
		enrichments_up_noncoreSmokeDEGs[[length(enrichments_up_noncoreSmokeDEGs)]]$padj),]
}

#Write to excel
#write.xlsx(enrichments_up_noncoreSmokeDEGs,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_upNoncore_byClusterDEGs.xlsx")

#Also print these in a single table so they can be sorted by smoke_DEGs if desired
enrichments_up_noncoreSmokeDEGs_oneTab<-cbind(enrichments_up_noncoreSmokeDEGs[[1]],"cluster_DEGs"=names(enrichments_up_noncoreSmokeDEGs)[[1]])[,c(1,8,2:7)]
for(i in 2:length(enrichments_up_noncoreSmokeDEGs)){
	enrichments_up_noncoreSmokeDEGs_oneTab<-rbind(enrichments_up_noncoreSmokeDEGs_oneTab,cbind(enrichments_up_noncoreSmokeDEGs[[i]],
		"cluster_DEGs"=names(enrichments_up_noncoreSmokeDEGs)[[i]])[,c(1,8,2:7)])
}
#write.xlsx(enrichments_up_noncoreSmokeDEGs_oneTab,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_upNoncore_oneTab.xlsx")

#Finally, flip, so there is one tab for each smoke_DEG cluster (this makes the most sense to me)
enrichments_up_noncoreSmokeDEGs_bySmokeDEGs<-split(enrichments_up_noncoreSmokeDEGs_oneTab,f=enrichments_up_noncoreSmokeDEGs_oneTab$smoke_DEGs)
enrichments_up_noncoreSmokeDEGs_bySmokeDEGs<-lapply(enrichments_up_noncoreSmokeDEGs_bySmokeDEGs,function(x)x[order(x$padj),])
write.xlsx(enrichments_up_noncoreSmokeDEGs_bySmokeDEGs,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_upNoncore.xlsx")










#### Non-core DOWN SMOKING DEGs
#Bring in smoking response "modules"
mods<-data.frame("Gene"=genes_smokeDown_noncore$gene,"Module"=genes_smokeDown_noncore$comparison,stringsAsFactors=F)
#The cycle enrichments through each of the clusters
enrichments_down_noncoreSmokeDEGs<-list()
for(i in 1:length(unique(genes_cellTypes$comparison))){
	genes<-genes_cellTypes[which(genes_cellTypes$comparison == unique(genes_cellTypes$comparison)[i]),]$gene
	enrichments_down_noncoreSmokeDEGs[[length(enrichments_down_noncoreSmokeDEGs) + 1]]<-DoEnrichment(genes=genes,gene2module=mods,background=background)
	colnames(enrichments_down_noncoreSmokeDEGs[[length(enrichments_down_noncoreSmokeDEGs)]])[1]<-"smoke_DEGs"
	names(enrichments_down_noncoreSmokeDEGs)[[length(enrichments_down_noncoreSmokeDEGs)]]<-unique(genes_cellTypes$comparison)[i]
	#Also, add in the number of smoking DEGs
	Nsmoke<-as.vector(table(mods$Module))
	#And the number of cell type signature DEGs
	NcellType<-length(genes)
	#And the number of smoking DEGs that overlap with the cell type signature DEGs
	overlap<-unlist(lapply(split(mods,f=mods$Module),function(x)sum(genes %in% x$Gene)))
	#And the proportion of smoking DEGs that overlap with the given cell type signature (get rid of this for now)
	#prop_smoke<-round(unlist(lapply(split(mods,f=mods$Module),function(x)sum(x$Gene %in% genes) / length(x$Gene))),2)
	#And the proportion of cell type signature genes that fall into the list of smoking DEGs
	prop_cellType<-round(unlist(lapply(split(mods,f=mods$Module),function(x)sum(x$Gene %in% genes) / length(genes))),2)
	#Merge these into a single table
	enrichments_down_noncoreSmokeDEGs[[length(enrichments_down_noncoreSmokeDEGs)]]<-cbind(enrichments_down_noncoreSmokeDEGs[[length(enrichments_down_noncoreSmokeDEGs)]],
		"N_smoke_DEGs"=Nsmoke,"N_cluster_DEGs"=NcellType,"N_overlap"=overlap,"Prop_overlap_cluster"=prop_cellType)
	#Finally, order by FDR
	enrichments_down_noncoreSmokeDEGs[[length(enrichments_down_noncoreSmokeDEGs)]]<-enrichments_down_noncoreSmokeDEGs[[length(enrichments_down_noncoreSmokeDEGs)]][order(
		enrichments_down_noncoreSmokeDEGs[[length(enrichments_down_noncoreSmokeDEGs)]]$padj),]
}

#Write to excel
#write.xlsx(enrichments_down_noncoreSmokeDEGs,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_downNoncore_byClusterDEGs.xlsx")

#Also print these in a single table so they can be sorted by smoke_DEGs if desired
enrichments_down_noncoreSmokeDEGs_oneTab<-cbind(enrichments_down_noncoreSmokeDEGs[[1]],"cluster_DEGs"=names(enrichments_down_noncoreSmokeDEGs)[[1]])[,c(1,8,2:7)]
for(i in 2:length(enrichments_down_noncoreSmokeDEGs)){
	enrichments_down_noncoreSmokeDEGs_oneTab<-rbind(enrichments_down_noncoreSmokeDEGs_oneTab,cbind(enrichments_down_noncoreSmokeDEGs[[i]],
		"cluster_DEGs"=names(enrichments_down_noncoreSmokeDEGs)[[i]])[,c(1,8,2:7)])
}
#write.xlsx(enrichments_down_noncoreSmokeDEGs_oneTab,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_downNoncore_oneTab.xlsx")

#Finally, flip, so there is one tab for each smoke_DEG cluster (this makes the most sense to me)
enrichments_down_noncoreSmokeDEGs_bySmokeDEGs<-split(enrichments_down_noncoreSmokeDEGs_oneTab,f=enrichments_down_noncoreSmokeDEGs_oneTab$smoke_DEGs)
enrichments_down_noncoreSmokeDEGs_bySmokeDEGs<-lapply(enrichments_down_noncoreSmokeDEGs_bySmokeDEGs,function(x)x[order(x$padj),])
write.xlsx(enrichments_down_noncoreSmokeDEGs_bySmokeDEGs,firstRow=T,colWidths="auto",file="SmokeDEG_clusterDEG_enrichments_downNoncore.xlsx")




















################################# MAKE ENRICHMENT HEAT MAP FOR SUBSETTED POPULTIONS FOR THE PAPER FIGURE, 
################################# USE NONCORE, AS CORE IS BROKEN OUT

#First, need to get a matrix of -log10 FDR values with clusters on the rows, and smoking DEGs on the columns
#First, toss c8, c8b, and c8c and then sort so in desired order
#Also, build a matrix of padj values for these
makePadjMat<-function(currData=NULL){
	for(i in 1:length(currData)){
		currData[[i]]<-currData[[i]][-which(as.character(currData[[i]]$cluster_DEGs) == "c8" |
			as.character(currData[[i]]$cluster_DEGs) == "c8b" | as.character(currData[[i]]$cluster_DEGs) == "c8c" |
			as.character(currData[[i]]$cluster_DEGs) == "c8a" | as.character(currData[[i]]$cluster_DEGs) == "c5M" |
			as.character(currData[[i]]$cluster_DEGs) == "c6a" | as.character(currData[[i]]$cluster_DEGs) == "c6b"|
			as.character(currData[[i]]$cluster_DEGs) == "c6c"),]
		currData[[i]]$cluster_DEGs<-factor(currData[[i]]$cluster_DEGs,
			levels=c("c1","c1c","c2","c4","c5","c8bc","c3","c7"))
		currData[[i]]<-currData[[i]][order(currData[[i]]$cluster_DEGs),]
		if(i == 1){
			currMat<-data.frame(matrix(nrow=nrow(currData[[i]]),ncol=length(currData)))
			rownames(currMat)<-currData[[i]]$cluster_DEGs
		}
		currMat[,i]<-as.numeric(format((-log10(currData[[i]]$padj)),scientific=F,digits=3))
		colnames(currMat)[i]<-names(currData)[[i]]
	}
	return(currMat)
}

#Now do core up, core down, all up, and all down matrices
coreUpMat<-makePadjMat(enrichments_up_coreSmokeDEGs_bySmokeDEGs)
coreDownMat<-makePadjMat(enrichments_down_coreSmokeDEGs_bySmokeDEGs)
allUpMat<-makePadjMat(enrichments_up_noncoreSmokeDEGs_bySmokeDEGs)
allDownMat<-makePadjMat(enrichments_down_noncoreSmokeDEGs_bySmokeDEGs)

#Add core up and down to the beginning of the all matrices so plotted together
allUpMat<-cbind(coreUpMat,allUpMat)
allDownMat<-cbind(coreDownMat,allDownMat)

#Now sort the columns in the order we want them (and exclude c8, c8a, c6a-c)
allUpMat<-allUpMat[,c(1:4,6,7,14,5,11)]
allDownMat<-allDownMat[,c(1:4,6,7,13,5,10)]

#Forgot to make the rows backwards
allUpMat<-allUpMat[rev(seq(nrow(allUpMat))),]
allDownMat<-allDownMat[rev(seq(nrow(allDownMat))),]

#For the down, there is 1 infinity. Just replace this with a large number
#allDownMat[4,10]<-150





##############################Now make enrichment heat map

#Choose up or down to do
mat<-allUpMat
#mat<-allDownMat

#Column colors
moduleOrdering_col<-c("black","#1619BA","mediumpurple3","#3281FF","tomato","orange","green3","darkturquoise","saddlebrown")

#Set module orderings
moduleOrdering_col<-moduleOrdering_col
#Get row colors for modules
modules<-color.scale(as.numeric(as.factor(moduleOrdering_col)),extremes=sort(unique(moduleOrdering_col)), color.spec="rgb")
ColSideColors<-data.frame(modules=modules,row.names=colnames(mat))

moduleOrdering_row<-rev(c("#1619BA","mediumpurple3","#3281FF","tomato","orange","green3","darkturquoise","saddlebrown"))
#moduleOrdering_row<-rev(c("#1619BA","mediumpurple3","#3281FF","tomato","orange","greenyellow","green3","wheat","gold4","yellow3",
#	"darkturquoise","saddlebrown"))
#Get row colors for modules
modules<-color.scale(as.numeric(as.factor(moduleOrdering_row)),extremes=sort(unique(moduleOrdering_row)), color.spec="rgb")
RowSideColors<-data.frame(modules=modules,row.names=rownames(mat))

####Make vector for cells that should be highlighted (comparisons with FDR < 0.05)
row.index<-c()
col.index<-c()
for(i in 1:ncol(mat)){
	row.index<-append(row.index,which(mat[,i] > 1.3))
	col.index<-append(col.index,rep(i,length(which(mat[,i] > 1.3))))
}
col.border<-rep("black",length(row.index))
border.width=rep(3,length(row.index))
highlightCellVec<-data.frame(row.index,col.index,col.border,border.width)




######### Scale values (one for now used for the paper
mycols = colorRampPalette(c("blue","white","red"))(1000)
breakscale = c(-8,seq(-1.2,1.2, length.out=length(mycols)-1), 8)
#distance.col = dist(t(as.matrix(mat)), method = 'canberra')
#cluster.col = hclust(distance.col, method = 'ward.D2')

#pdf("Heatmap_noncoreSmokeDEG_clusterDEG_enrichments_UP.pdf")
heatmap3(as.matrix(mat),scale="column",col=mycols,ColSideColors=as.matrix(ColSideColors),
	RowSideColors=as.matrix(RowSideColors),cexRow=1,cexCol=1,Rowv=NA,Colv=NA, breaks=breakscale,
	highlightCell=highlightCellVec)
#dev.off()



########## Breakscaling unscaled values (this actually looks pretty good)
mycols = colorRampPalette(c("blue","white","red"))(1000)
breakscale = c(0,seq(0,20, length.out=length(mycols)-1), 200)
#distance.col = dist(t(as.matrix(mat)), method = 'canberra')
#cluster.col = hclust(distance.col, method = 'ward.D2')

#pdf("Heatmap_SmokeDEG_clusterDEG_enrichments_UP.pdf")
heatmap3(as.matrix(mat),scale="none",col=mycols,ColSideColors=as.matrix(ColSideColors),
	RowSideColors=as.matrix(RowSideColors),cexRow=1,cexCol=1,Rowv=NA,Colv=NA, breaks=breakscale,
	highlightCell=highlightCellVec)
#dev.off()


























################################# FINALLY, REPETE FOR RARE CELL FIGURE

#First, need to get a matrix of -log10 FDR values with clusters on the rows, and smoking DEGs on the columns
#First, toss c8, c8b, and c8c and then sort so in desired order
#Also, build a matrix of padj values for these
makePadjMat<-function(currData=NULL){
	for(i in 1:length(currData)){
		currData[[i]]<-currData[[i]][-which(as.character(currData[[i]]$cluster_DEGs) == "c8" |
			as.character(currData[[i]]$cluster_DEGs) == "c8b" | as.character(currData[[i]]$cluster_DEGs) == "c8c" |
			as.character(currData[[i]]$cluster_DEGs) == "c8a" | as.character(currData[[i]]$cluster_DEGs) == "c5M" |
			as.character(currData[[i]]$cluster_DEGs) == "c1" |
			as.character(currData[[i]]$cluster_DEGs) == "c2" | as.character(currData[[i]]$cluster_DEGs) == "c3" |
			as.character(currData[[i]]$cluster_DEGs) == "c7" | as.character(currData[[i]]$cluster_DEGs) == "c8bc"),]
		currData[[i]]$cluster_DEGs<-factor(currData[[i]]$cluster_DEGs,
			levels=c("c1c","c4","c5","c6b","c6a","c6c"))
		currData[[i]]<-currData[[i]][order(currData[[i]]$cluster_DEGs),]
		if(i == 1){
			currMat<-data.frame(matrix(nrow=nrow(currData[[i]]),ncol=length(currData)))
			rownames(currMat)<-currData[[i]]$cluster_DEGs
		}
		currMat[,i]<-as.numeric(format((-log10(currData[[i]]$padj)),scientific=F,digits=3))
		colnames(currMat)[i]<-names(currData)[[i]]
	}
	return(currMat)
}

#Now do core up, core down, all up, and all down matrices
coreUpMat<-makePadjMat(enrichments_up_coreSmokeDEGs_bySmokeDEGs)
coreDownMat<-makePadjMat(enrichments_down_coreSmokeDEGs_bySmokeDEGs)
allUpMat<-makePadjMat(enrichments_up_allSmokeDEGs_bySmokeDEGs)
allDownMat<-makePadjMat(enrichments_down_allSmokeDEGs_bySmokeDEGs)

#Add core up and down to the beginning of the all matrices so plotted together
allUpMat<-cbind(coreUpMat,allUpMat)
allDownMat<-cbind(coreDownMat,allDownMat)

#Now sort the columns in the order we want them (and exclude everything except rare cells)
allUpMat<-allUpMat[,c(9,8,10)]
allDownMat<-allDownMat[,c(8,9)]

#Forgot to make the rows backwards
allUpMat<-allUpMat[rev(seq(nrow(allUpMat))),]
allDownMat<-allDownMat[rev(seq(nrow(allDownMat))),]

#For the down, there is 1 infinity. Just replace this with a large number
allDownMat[2,2]<-140





##############################Now make enrichment heat map

#Choose up or down to do
mat<-allUpMat
#mat<-allDownMat

#Column colors
#UP
moduleOrdering_col<-c("yellow3","gold4","wheat")
#DOWN
#moduleOrdering_col<-c("gold4","wheat")

#Set module orderings
moduleOrdering_col<-moduleOrdering_col
#Get row colors for modules
modules<-color.scale(as.numeric(as.factor(moduleOrdering_col)),extremes=sort(unique(moduleOrdering_col)), color.spec="rgb")
ColSideColors<-data.frame(modules=modules,row.names=colnames(mat))

moduleOrdering_row<-rev(c("mediumpurple3","tomato","orange","yellow3","gold4","wheat"))
modules<-color.scale(as.numeric(as.factor(moduleOrdering_row)),extremes=sort(unique(moduleOrdering_row)), color.spec="rgb")
RowSideColors<-data.frame(modules=modules,row.names=rownames(mat))

####Make vector for cells that should be highlighted (in this case, ALL cells)
row.index<-c()
col.index<-c()
for(i in 1:ncol(mat)){
	row.index<-append(row.index,which(mat[,i] >= 0))
	col.index<-append(col.index,rep(i,length(which(mat[,i] >= 0))))
}
col.border<-rep("black",length(row.index))
border.width=rep(3,length(row.index))
highlightCellVec<-data.frame(row.index,col.index,col.border,border.width)



#Scale so that it's white if p-value < 0.05 (-log10 of 1.3) and then becomes blazing red at 12
mycols = colorRampPalette(c("white","red"))(1000)
breakscale = c(0,seq(1.3,12, length.out=length(mycols)-1), 150)
#distance.col = dist(t(as.matrix(mat)), method = 'canberra')
#cluster.col = hclust(distance.col, method = 'ward.D2')

#pdf("Heatmap_SmokeDEG_clusterDEG_enrichments_rareOnly_UP_noScale.pdf")
heatmap3(as.matrix(mat),scale="none",col=mycols,ColSideColors=as.matrix(ColSideColors),
	RowSideColors=as.matrix(RowSideColors),cexRow=1,cexCol=1,Rowv=NA,Colv=NA, breaks=breakscale,
	highlightCell=highlightCellVec)
#dev.off()














