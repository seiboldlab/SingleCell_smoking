library(Seurat, lib.loc="/Users/jacksonna/R-dev")
library(plotrix)
library(scales)
library(ggplot2)
library(cowplot)
library(openxlsx)
source("CommonFunctions.r")

#Download the processed Seurat R object, "Processed_invitro_seurat.Rdata" from GEO,
#obtainable from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134174
#Then load the object it contains ("T15_int") into the session.
#T15_int is a Seurat dataset that contains the processed, integrated dataset
load("Processed_invitro_seurat.Rdata")



################### Differential expression and Enrichments for the major clusters
######################################################Find markers and do Enrichments

##########################First, do each against one another
#bsub -e test.err -o test.out -R "rusage[mem=42000]" "module load R/3.5.1 && env R_MAX_VSIZE=700Gb R CMD BATCH diffexp1.r"
#bsub -e test.err -o test.out -R "rusage[mem=42000]" "module load R/3.5.1 && env R_MAX_VSIZE=700Gb R CMD BATCH diffexp1_LR.r"

#Test with logistic regression, with smoking as a latent variable
Idents(T15_int) <- T15_int@meta.data$clusters_10
clusterList<-sort(as.character(unique(Idents(T15_int))))
markers.list <- list()
for(i in 1:length(clusterList)) {
    markers.list[[i]] <- FindMarkers(T15_int, ident.1=clusterList[i], min.pct=0.1, logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "LR", latent.vars = "smoke")
}


#Get comparison descriptions
Idents(T15_int) <- T15_int@meta.data$clusters_10
clusterList<-sort(as.character(unique(Idents(T15_int))))
comparisonVec<-clusterList

#Get vector of clusters for calculating pct values
pctIndexList<-list()
for(i in 1:length(clusterList)){
	currClusterIndex<-which(clusterList == clusterList[i])
	pctIndexList[[length(pctIndexList) + 1]]<-list(clusterList[i],clusterList[-currClusterIndex])
}

#Get counts and cluster assignments
counts<-as.matrix(GetAssayData(T15_int,assay="RNA",slot="data"))
clusterAssignments<-T15_int@meta.data$clusters_10

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
write.table(T15_1againstAll_DEGs_forEnrichr,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,file="Enrichr/T15_1againstAll_DEGs_forEnrichr.txt")

#Do enrichments
DEG_table<-read.table("Enrichr/T15_1againstAll_DEGs_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"T15_1againstAll_DEGs"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)





















################################ Dot plot of marker list
##Read in Kate's biomarkers
Kate_markers<-scan("Gene_lists/Kate_biomarker_list7.txt",what="character")
Kate_markers<-Kate_markers[which(Kate_markers %in% rownames(GetAssayData(T15_int,assay="SCT",slot="data")))]

##CLUSTERS
#Set desired order of clusters to ident
temp<-data.frame(T15_int@meta.data$clusters_10,row.names=rownames(T15_int@meta.data))
colnames(temp)<-"clusters_10"
Idents(T15_int)<-factor(temp[,1],levels=c("c1","c1c","c2","c4","c5","c8","c6c","c6ab","c3","c7"))

#pdf("T15_biomarkers_dotPlot_clusters_10_ordered.pdf")
dev.new(width=6.65,height=9.4)
dotPlots<-DotPlot(T15_int,features=Kate_markers,plot.legend=T,col.min=-1,col.max=1,dot.scale=9.5,
	cols = c("lightgrey", "blue"),do.return=T,assay="SCT")
dotPlots + theme(axis.text.x=element_text(angle = -45, hjust = 0, size = 12)) + coord_flip()
















################################### Get distributions of cells across clusters ##################################

########## Export table with counts for each cluster, broken out by donor
#Create a table with all the values I want to count up
T15_metadata<-T15_int@meta.data[,c(4:5,50)]

#Generate count table for each combination
library(plyr)
T15_metadata_counts<-count(T15_metadata)
colnames(T15_metadata_counts)[4]<-"cellCounts"

#Export
write.table(T15_metadata_counts,sep="\t",quote=F,row.names=F,file="T15_metadata_cellFreq_16aVersion.txt")




#Barplot with total cells per donor
#pdf("T15_barplots_donorOverall.pdf")
dev.new(height=5,width=5)
colorVec<-c("red","black","red","black","black","black","red","red")
barplot(NumRemainVec,las=1,names.arg="",col=colorVec,ylim=c(0,4000))
text(x=c(0.8,2,3.2,4.4,5.6,6.8,8,9.2),y=par()$usr[3]-0.08*(par()$usr[4]-par()$usr[3]),labels=names(table(T15_int@meta.data$donor)),srt=45,adj=0.8,xpd=T)
box(which = "plot", bty = "l")













###########Calculate out and plot differences in clusters between smokers and non-smokers, stratified by SMG and epithelium

#Here, for each donor, we want to calculate the proportion of all its cells in each cluster. Then for each cluster, I'll plot these proportions, grouping 
#smokers and non-smoker together
#As for the clusters, I will do this for all the clusters together, but then also by stratifying clusters into SMG and surface

#Cull out the metadata to exclude light smokers and T89
T15_metadata_counts_culled<-T15_metadata_counts[which(T15_metadata_counts$smoke != "light" & T15_metadata_counts$donor != "T89"),]

##Add a row for T137 c6c, which is missing because there were zero cells (for when using clusters_10)
#T15_metadata_counts_culled<-rbind(T15_metadata_counts_culled[1:37,],data.frame("donor"="T137","smoke"="never","clusters_10"="c6c","cellCounts"=0),
#	T15_metadata_counts_culled[38:nrow(T15_metadata_counts_culled),])
#Add rows for T137 c6a and c6c and T154 c6b abd c8a, which is missing because there were zero cells (for when using clusters_12)
T15_metadata_counts_culled<-rbind(T15_metadata_counts_culled[1:42,],data.frame("donor"="T137","smoke"="never","clusters_12"="c6a","cellCounts"=0),
	T15_metadata_counts_culled[43,],data.frame("donor"="T137","smoke"="never","clusters_12"="c6c","cellCounts"=0),
	T15_metadata_counts_culled[44:65,],data.frame("donor"="T154","smoke"="heavy","clusters_12"="c6b","cellCounts"=0),
	T15_metadata_counts_culled[66:67,],data.frame("donor"="T154","smoke"="heavy","clusters_12"="c8a","cellCounts"=0),
	T15_metadata_counts_culled[68:nrow(T15_metadata_counts_culled),])
#Sort table by smoke, then by donor, then by cluster
T15_metadata_counts_culled$smoke<-factor(T15_metadata_counts_culled$smoke,levels=c("never","heavy"))
T15_metadata_counts_culled<-T15_metadata_counts_culled[order(T15_metadata_counts_culled$smoke,T15_metadata_counts_culled$donor,T15_metadata_counts_culled$clusters_12),]

#Vectors of populations to run
surfVec<-c("c1","c1c","c2","c4","c5","c6a","c6b","c6c","c8a","c8bc")
smgVec<-c("c3","c7")
allVec<-sort(c(surfVec,smgVec))

##Get versions of the table for surface and SMG-only groups
#T15_metadata_counts_culled_surf<-T15_metadata_counts_culled[which(T15_metadata_counts_culled$clusters_12 %in% surfVec),]
#T15_metadata_counts_culled_smg<-T15_metadata_counts_culled[which(T15_metadata_counts_culled$clusters_12 %in% smgVec),]

#####Get vectors of total cell numbers for each donor, based on different groups of clusters
getTotals<-function(dataset){
	totals<-c()
	for(i in 1:length(unique(dataset$donor))){
		totals<-append(totals,sum(dataset$cellCounts[which(dataset$donor == unique(dataset$donor)[i])]))
	}
	names(totals)<-unique(dataset$donor)
	return(totals)
}
totals_all<-getTotals(T15_metadata_counts_culled)

#totals_surf<-getTotals(T15_metadata_counts_culled_surf)
#totals_smg<-getTotals(T15_metadata_counts_culled_smg)

#Get donor colors
colorVec<-c("pink","tan","greenyellow","yellow","lightblue","grey","midnightblue","black","orangered4","darkgreen","purple4","slateblue")

#Now, for one of the three datasets above, plot proportions for each cluster
dataset<-T15_metadata_counts_culled
#dataset<-T15_metadata_counts_culled_surf
#dataset<-T15_metadata_counts_culled_smg

totals<-totals_all
#totals<-totals_surf
#totals<-totals_smg

#If doing one-sided...
altVec<-c("greater","greater","greater","less","less","greater","less","greater","greater","less")
#altVec<-c("greater","greater","greater","less","greater","less","greater","less")
#altVec<-c("less","greater")

#Now make box plots
#pdf("T15_boxplots_cellProportionsBySmoking_clusters12.pdf")
dev.new(width=6,height=6.5)
par(mfrow=c(3,4),bty="l")
for(i in 1:length(unique(dataset$clusters_12))){
	currData<-dataset[which(dataset$clusters_12 == unique(dataset$clusters_12)[i]),]
	currData$cellProportions<-currData$cellCounts / totals
	boxplot(currData$cellProportions~currData$smoke,las=1,names=rep("",2),main=unique(dataset$clusters_12)[i],
		border=c("black","red"),medcol=c("black","red"),yaxt="n",
		ylim=c(min(currData$cellProportions),max(currData$cellProportions)*1.1))
	points(x=c(rep(1,6),rep(2,6)),y=currData$cellProportions,col=colorVec,pch=16,cex=1.5)
	text(x=1:2,y=par()$usr[3]-0.12*(par()$usr[4]-par()$usr[3]),labels=c("Nonsmokers","Smokers"),
		srt=45,adj=1,xpd=T,cex=1)
	mtext(paste("p = ",formatC(t.test(currData$cellProportions[which(currData$smoke=="heavy")],
	currData$cellProportions[which(currData$smoke=="never")],alternative="two.sided")$p.value,digits=2),sep=""),
	side=3,line=0,at=2.3,adj=1,cex=0.7)
	axis(2,at=round(seq(from=min(currData$cellProportions),to=max(currData$cellProportions),length.out=5),digits=6),las=1,
		labels=formatC(seq(from=min(currData$cellProportions),to=max(currData$cellProportions),length.out=5),format="f",digits=3))
	title(ylab="Proportion of cells", line=3.3, cex.lab=1)
}



#Now also make bar plots
#First get raw shift in cell proportions for each cell type and donor
#For each heavy value, get shift from the mean never values
props<-lapply(split(dataset,f=dataset$clusters_12),function(x)x$cellCounts / totals)
props_shift<-t(apply(sapply(props,function(x)x[unique(dataset$donor[which(dataset$smoke == "heavy")])]),1,
	function(x)x - apply(sapply(props,function(x)x[unique(dataset$donor[which(dataset$smoke == "never")])]),2,mean)))

#Then get means and standard errors of this
props_shift_mean<-apply(props_shift,2,mean)
props_shift_se<-apply(props_shift,2,function(x)sd(x)/sqrt(length(x)))
#Finally, select and order the clusters as I want
props_shift_mean<-props_shift_mean[c(6,2,10,1,9,7,8,5,4,3,12)]
props_shift_se<-props_shift_se[c(6,2,10,1,9,7,8,5,4,3,12)]
colorVec<-c("orange","mediumpurple3","saddlebrown","#1619BA","wheat","yellow4","yellow3","tomato",
	"darkturquoise","#3281FF","green3")

library(gplots)
#pdf("T15_barplots_cellProportionsBySmoking_clusters12.pdf",height=8,width=5)
dev.new(height=4,width=6)
barplot2(rev(props_shift_mean),las=1,col=rev(colorVec),ylim=c(1,11),names.arg="",xaxt="n",
	xlab="Shift in proportion of Cells",horiz=T, xlim=c(-0.05,0.05),ci.l = rev(props_shift_mean - props_shift_se), 
	ci.u = rev(props_shift_mean + props_shift_se), plot.ci=T, ci.color="black",space=c(0.1,1),width=0.7)
abline(v=0)
axis(1,c(-0.05,-0.025,0,0.025,0.05),labels=c(-0.05,-0.025,0,0.025,0.05))





#Here is a version of bar plots that looks at the % shift in the proportions, rather than a raw shift in the proportions
#First get raw shift in cell proportions for each cell type and donor
#For each heavy value, get shift from the mean never values
props<-lapply(split(dataset,f=dataset$clusters_12),function(x)x$cellCounts / totals)
props_shift<-t(apply(sapply(props,function(x)x[unique(dataset$donor[which(dataset$smoke == "heavy")])]),1,
	function(x)x - apply(sapply(props,function(x)x[unique(dataset$donor[which(dataset$smoke == "never")])]),2,mean)))

#If you want to convert this to a % shift in the proportion
#For this, divide each proportional shift by the mean proportion for the never smokers
props_shift_temp<-props_shift
for(i in 1:ncol(props_shift)){
	props_shift[,i]<-(props_shift[,i]/apply(sapply(props,function(x)x[unique(dataset$donor[which(dataset$smoke == "never")])]),2,mean)[i]) * 100
}

#Then get means and standard errors of this
props_shift_mean<-apply(props_shift,2,mean)
props_shift_se<-apply(props_shift,2,function(x)sd(x)/sqrt(length(x)))
#Finally, select and order the clusters as I want
props_shift_mean<-props_shift_mean[c(6,10,2,9,7,1,5,4,3,8,12)]
props_shift_se<-props_shift_se[c(6,10,2,9,7,1,5,4,3,8,12)]
colorVec<-c("orange","saddlebrown","mediumpurple3","wheat","yellow4","#1619BA","tomato",
	"darkturquoise","#3281FF","yellow3","green3")

library(gplots)
#pdf("T15_barplots_percentCellProportionsBySmoking_clusters12.pdf",height=8,width=5)
dev.new(height=4,width=6)
barplot2(rev(props_shift_mean),las=1,col=rev(colorVec),ylim=c(1,11),names.arg="",xaxt="n",
	xlab="Percent shift in proportion of Cells",horiz=T, xlim=c(-75,100),ci.l = rev(props_shift_mean - props_shift_se), 
	ci.u = rev(props_shift_mean + props_shift_se), plot.ci=T, ci.color="black",space=c(0.1,1),width=0.7)
abline(v=0)
axis(1,c(-75,-50,-25,0,25,50,75),labels=c(-75,-50,-25,0,25,50,75))


















############################## Get markers for the main clusters

dataset<-T15_int
cluster.meta.data<-T15_int@meta.data$clusters_10
downsample<-"median" #2,021 cells
min.pct=0.1
logfc.threshold=0.25
test.use="LR"
latent.vars = "smoke"

DEGsCons<-GetConservedMarkers(dataset=dataset,cluster.meta.data=cluster.meta.data,min.pct=min.pct,
	logfc.threshold=logfc.threshold,test.use=test.use,latent.vars=latent.vars,
	downsample=downsample)







########################### Repeat by adding in three rare cell clusters

#Bring in rare cell subclusters
load("Subclustering_rda/T15_rare3.rda")

#Now combine subclusters and clusters I care about together 
treat.col<-data.frame("clusters_10_plusRare3"=T15_int@meta.data$clusters_10,row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col[,1]<-as.character(treat.col[,1])
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_rare3@meta.data)[which(T15_rare3@meta.data$clusters_3 == "c6a")])],]<-"c6a"
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_rare3@meta.data)[which(T15_rare3@meta.data$clusters_3 == "c6b")])],]<-"c6b"
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_rare3@meta.data)[which(T15_rare3@meta.data$clusters_3 == "c6c")])],]<-"c6c"
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="clusters_10_plusRare3")

dataset<-T15_int
cluster.meta.data<-T15_int@meta.data$clusters_10_plusRare3
downsample<-"median" #1685 cells
min.pct=0.1
logfc.threshold=0.25
test.use="LR"
latent.vars = "smoke"

DEGsCons<-GetConservedMarkers(dataset=dataset,cluster.meta.data=cluster.meta.data,min.pct=min.pct,
	logfc.threshold=logfc.threshold,test.use=test.use,latent.vars=latent.vars,
	downsample=downsample)




########################### Repeat by adding in three SMG basal subclusters

#Bring in SMG clusters
load("Subclustering_rda/T15_smgb.rda")

#Now combine subclusters and clusters I care about together 
treat.col<-data.frame("clusters_10_plusSMGB"=T15_int@meta.data$clusters_10,row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col[,1]<-as.character(treat.col[,1])
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_smgb@meta.data)[which(T15_smgb@meta.data$clusters_3 == "c3a")])],]<-"c3a"
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_smgb@meta.data)[which(T15_smgb@meta.data$clusters_3 == "c3b")])],]<-"c3b"
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_smgb@meta.data)[which(T15_smgb@meta.data$clusters_3 == "c3c")])],]<-"c3c"
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="clusters_10_plusSMGB")

dataset<-T15_int
cluster.meta.data<-T15_int@meta.data$clusters_10_plusSMGB
downsample<-"median" #1,270 cells
min.pct=0.1
logfc.threshold=0.25
test.use="LR"
latent.vars = "smoke"

DEGsCons<-GetConservedMarkers(dataset=dataset,cluster.meta.data=cluster.meta.data,min.pct=min.pct,
	logfc.threshold=logfc.threshold,test.use=test.use,latent.vars=latent.vars,
	downsample=downsample)



########################### Repeat by adding in three ciliated subclusters

#Bring in ciliated subclusters
load("Subclustering_rda/T15_cil.rda")

#Now combine subclusters and clusters I care about together 
treat.col<-data.frame("clusters_10_plusCil"=T15_int@meta.data$clusters_10,row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col[,1]<-as.character(treat.col[,1])
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_cil@meta.data)[which(T15_cil@meta.data$clusters_3 == "c8a")])],]<-"c8a"
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_cil@meta.data)[which(T15_cil@meta.data$clusters_3 == "c8b")])],]<-"c8b"
treat.col[rownames(T15_int@meta.data)[which(rownames(T15_int@meta.data) %in% rownames(T15_cil@meta.data)[which(T15_cil@meta.data$clusters_3 == "c8c")])],]<-"c8c"
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="clusters_10_plusCil")

dataset<-T15_int
cluster.meta.data<-T15_int@meta.data$clusters_10_plusCil
downsample<-"median" #1,475 cells
min.pct=0.1
logfc.threshold=0.25
test.use="LR"
latent.vars = "smoke"

DEGsCons<-GetConservedMarkers(dataset=dataset,cluster.meta.data=cluster.meta.data,min.pct=min.pct,
	logfc.threshold=logfc.threshold,test.use=test.use,latent.vars=latent.vars,
	downsample=downsample)




########################### Repeat by adding in the monocle mature secretory cluster

#First, bring in the monocle-derived mature secretory cell group and add as metadata
matureSecCells<-scan("monocle/ANALYSIS_lumSec_reintegratedData/Mature_secretory_cells_monocle.txt",what="character")
treat.col<-data.frame("Mature_secretory_cells_monocle"=rep(FALSE,nrow(T15_int@meta.data)),row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col[which(rownames(treat.col) %in% matureSecCells),]<-TRUE
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="Mature_secretory_cells_monocle")

#Now combine subclusters and clusters I care about together 
treat.col<-data.frame("clusters_20_plusMatureSec"=T15_int@meta.data$clusters_10,row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col[,1]<-as.character(treat.col[,1])
treat.col[rownames(T15_int@meta.data)[which(T15_int@meta.data$Mature_secretory_cells_monocle==TRUE)],]<-"c5M"
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="clusters_20_plusMatureSec")

dataset<-T15_int
cluster.meta.data<-T15_int@meta.data$clusters_20_plusMatureSec
downsample<-"median" #1,685 cells
min.pct=0.1
logfc.threshold=0.25
test.use="LR"
latent.vars = "smoke"

DEGsCons<-GetConservedMarkers(dataset=dataset,cluster.meta.data=cluster.meta.data,min.pct=min.pct,
	logfc.threshold=logfc.threshold,test.use=test.use,latent.vars=latent.vars,
	downsample=downsample)














######################################  CULLING THE MARKER LIST AND MAKING HEAT MAPS FOR FIGURE S1

#Bring in the conserved markers list from above
DEGsCons<-read.table("ConservedBiomarkers/Biomarkers_main/Conserved_DEGs.txt",header=T,stringsAsFactors=F)

#Also bring in the 1againstAll DEGs, in case we need to grab some of these
oaa_DEGs<-read.table("DEG_Tables/T15_1againstAll_DEGs.txt",header=T,stringsAsFactors=F)

#Get culled conserved markers gene list for plotting
p_adj_FDR=1e-5
avg_logFC=0.25
pct.1=0.2
pct_ratio="q25"
minGene=20
maxGene=200
rank.by="pct_ratio"
genesList<-GetCulledConservedMarkersList(DEGsCons=DEGsCons,p_adj_FDR=p_adj_FDR,avg_logFC=avg_logFC,pct.1=pct.1,pct_ratio=pct_ratio,
	minGene=minGene,maxGene=maxGene,rank.by=rank.by)
	
#For c1c, c2, c3, and c4, add in extra genes from oaa_DEGs
#Do c1c
oaa_DEGs_subset<-oaa_DEGs[which(oaa_DEGs$comparison == "c1c" & oaa_DEGs$p_adj_FDR < 0.05),]
newGeneVec<-c()
for(i in 1:nrow(oaa_DEGs_subset)){
	#For each gene, if Log fold change in c4 is the max seen anywhere, and if the FDR in c4 is the
	#min seen anywhere, and if the gene is not already among the conserved markers anywhere, then keep the gene
	if(oaa_DEGs_subset$avg_logFC[which(oaa_DEGs_subset$gene == oaa_DEGs_subset$gene[i])] == 
		max(oaa_DEGs$avg_logFC[which(oaa_DEGs$gene == oaa_DEGs_subset$gene[i])]) &
		oaa_DEGs_subset$p_adj_FDR[which(oaa_DEGs_subset$gene == oaa_DEGs_subset$gene[i])] == 
		min(oaa_DEGs$p_adj_FDR[which(oaa_DEGs$gene == oaa_DEGs_subset$gene[i])]) &
		!(oaa_DEGs_subset$gene[i] %in% unlist(genesList))){
		newGeneVec<-append(newGeneVec,oaa_DEGs_subset$gene[i])
	}
}
genesList$markers_c1c<-c(genesList$markers_c1c,newGeneVec[1:11])

#Do c2
oaa_DEGs_subset<-oaa_DEGs[which(oaa_DEGs$comparison == "c2" & oaa_DEGs$p_adj_FDR < 0.05),]
newGeneVec<-c()
for(i in 1:nrow(oaa_DEGs_subset)){
	#For each gene, if Log fold change in c4 is the max seen anywhere, and if the FDR in c4 is the
	#min seen anywhere, and if the gene is not already among the conserved markers anywhere, then keep the gene
	if(oaa_DEGs_subset$avg_logFC[which(oaa_DEGs_subset$gene == oaa_DEGs_subset$gene[i])] == 
		max(oaa_DEGs$avg_logFC[which(oaa_DEGs$gene == oaa_DEGs_subset$gene[i])]) &
		oaa_DEGs_subset$p_adj_FDR[which(oaa_DEGs_subset$gene == oaa_DEGs_subset$gene[i])] == 
		min(oaa_DEGs$p_adj_FDR[which(oaa_DEGs$gene == oaa_DEGs_subset$gene[i])]) &
		!(oaa_DEGs_subset$gene[i] %in% unlist(genesList))){
		newGeneVec<-append(newGeneVec,oaa_DEGs_subset$gene[i])
	}
}
genesList$markers_c2<-c(genesList$markers_c2,newGeneVec[1:2])

#Do c3, only keep genes with high pct.2
oaa_DEGs_subset<-oaa_DEGs[which(oaa_DEGs$comparison == "c3" & oaa_DEGs$p_adj_FDR < 0.05),]
newGeneVec<-c()
for(i in 1:nrow(oaa_DEGs_subset)){
	#For each gene, if Log fold change in c4 is the max seen anywhere, and if the FDR in c4 is the
	#min seen anywhere, and if the gene is not already among the conserved markers anywhere, then keep the gene
	if(oaa_DEGs_subset$avg_logFC[which(oaa_DEGs_subset$gene == oaa_DEGs_subset$gene[i])] == 
		max(oaa_DEGs$avg_logFC[which(oaa_DEGs$gene == oaa_DEGs_subset$gene[i])]) &
		oaa_DEGs_subset$p_adj_FDR[which(oaa_DEGs_subset$gene == oaa_DEGs_subset$gene[i])] == 
		min(oaa_DEGs$p_adj_FDR[which(oaa_DEGs$gene == oaa_DEGs_subset$gene[i])]) &
		!(oaa_DEGs_subset$gene[i] %in% unlist(genesList))){
		newGeneVec<-append(newGeneVec,oaa_DEGs_subset$gene[i])
	}
}
genesList$markers_c3<-c(genesList$markers_c3,newGeneVec[1:50])
genesList$markers_c3<-genesList$markers_c3[which(genesList$markers_c3 %in% oaa_DEGs_subset[which(oaa_DEGs_subset$pct.2 > 0.15),]$gene)][1:20]	

#Do c4
oaa_DEGs_subset<-oaa_DEGs[which(oaa_DEGs$comparison == "c4" & oaa_DEGs$p_adj_FDR < 0.05),]
newGeneVec<-c()
for(i in 1:nrow(oaa_DEGs_subset)){
	#For each gene, if Log fold change in c4 is the max seen anywhere, and if the FDR in c4 is the
	#min seen anywhere, and if the gene is not already among the conserved markers anywhere, then keep the gene
	if(oaa_DEGs_subset$avg_logFC[which(oaa_DEGs_subset$gene == oaa_DEGs_subset$gene[i])] == 
		max(oaa_DEGs$avg_logFC[which(oaa_DEGs$gene == oaa_DEGs_subset$gene[i])]) &
		oaa_DEGs_subset$p_adj_FDR[which(oaa_DEGs_subset$gene == oaa_DEGs_subset$gene[i])] == 
		min(oaa_DEGs$p_adj_FDR[which(oaa_DEGs$gene == oaa_DEGs_subset$gene[i])]) &
		!(oaa_DEGs_subset$gene[i] %in% unlist(genesList))){
		newGeneVec<-append(newGeneVec,oaa_DEGs_subset$gene[i])
	}
}
genesList$markers_c4<-c(genesList$markers_c4,newGeneVec[1:16])


#For c6bc, only keep genes with low pct.2
oaa_DEGs_subset<-oaa_DEGs[which(oaa_DEGs$comparison == "c6ab" & oaa_DEGs$p_adj_FDR < 0.05),]
newGeneVec<-c()
for(i in 1:nrow(oaa_DEGs_subset)){
	#For each gene, if Log fold change in c4 is the max seen anywhere, and if the FDR in c4 is the
	#min seen anywhere, and if the gene is not already among the conserved markers anywhere, then keep the gene
	if(oaa_DEGs_subset$avg_logFC[which(oaa_DEGs_subset$gene == oaa_DEGs_subset$gene[i])] == 
		max(oaa_DEGs$avg_logFC[which(oaa_DEGs$gene == oaa_DEGs_subset$gene[i])]) &
		oaa_DEGs_subset$p_adj_FDR[which(oaa_DEGs_subset$gene == oaa_DEGs_subset$gene[i])] == 
		min(oaa_DEGs$p_adj_FDR[which(oaa_DEGs$gene == oaa_DEGs_subset$gene[i])]) &
		!(oaa_DEGs_subset$gene[i] %in% unlist(genesList))){
		newGeneVec<-append(newGeneVec,oaa_DEGs_subset$gene[i])
	}
}
genesList$markers_c6ab<-c(genesList$markers_c6ab,newGeneVec[1:50])
genesList$markers_c6ab<-genesList$markers_c6ab[which(genesList$markers_c6ab %in% oaa_DEGs_subset[which(oaa_DEGs_subset$pct.2 < 0.1),]$gene)][1:20]
	
#For c7, only keep genes with high pct.2
oaa_DEGs_subset<-oaa_DEGs[which(oaa_DEGs$comparison == "c7" & oaa_DEGs$p_adj_FDR < 0.05),]
genesList$markers_c7<-genesList$markers_c7[which(genesList$markers_c7 %in% oaa_DEGs_subset[which(oaa_DEGs_subset$pct.2 > 0.12),]$gene)][1:20]	

#For c6c, only keep genes with low pct.2
oaa_DEGs_subset<-oaa_DEGs[which(oaa_DEGs$comparison == "c6c" & oaa_DEGs$p_adj_FDR < 0.05),]
genesList$markers_c6c<-genesList$markers_c6c[which(genesList$markers_c6c %in% oaa_DEGs_subset[which(oaa_DEGs_subset$pct.2 < 0.05),]$gene)][1:20]

#Cull rest	
genesList$markers_c1<-genesList$markers_c1[1:20]
genesList$markers_c5<-genesList$markers_c5[1:20]
genesList$markers_c8<-genesList$markers_c8[1:20]

























#################################################### BIOMARKER HEATMAP (by cluster)
sct_exp<-as.matrix(GetAssayData(T15_int,assay="SCT"))

cols_clusters<-c("#1619BA","mediumpurple3","#3281FF","darkturquoise","tomato","orange","gold4","tan","saddlebrown","green3")

#First, need to isolate those genes we want and cells we want in the order that we want them
############### For su
cellOrdering<-c(
	sample(which(T15_int@meta.data$clusters_10 == "c1"),
		length(which(T15_int@meta.data$clusters_10 == "c1")),replace=FALSE),
	sample(which(T15_int@meta.data$clusters_10 == "c1c"),
		length(which(T15_int@meta.data$clusters_10 == "c1c")),replace=FALSE),
	sample(which(T15_int@meta.data$clusters_10 == "c2"),
		length(which(T15_int@meta.data$clusters_10 == "c2")),replace=FALSE),
	sample(which(T15_int@meta.data$clusters_10 == "c4"),
		length(which(T15_int@meta.data$clusters_10 == "c4")),replace=FALSE),
	sample(which(T15_int@meta.data$clusters_10 == "c5"),
		length(which(T15_int@meta.data$clusters_10 == "c5")),replace=FALSE),
	sample(which(T15_int@meta.data$clusters_10 == "c8"),
		length(which(T15_int@meta.data$clusters_10 == "c8")),replace=FALSE),
	sample(which(T15_int@meta.data$clusters_10 == "c6c"),
		length(which(T15_int@meta.data$clusters_10 == "c6c")),replace=FALSE),
	sample(which(T15_int@meta.data$clusters_10 == "c6ab"),
		length(which(T15_int@meta.data$clusters_10 == "c6ab")),replace=FALSE),
	sample(which(T15_int@meta.data$clusters_10 == "c3"),
		length(which(T15_int@meta.data$clusters_10 == "c3")),replace=FALSE),
	sample(which(T15_int@meta.data$clusters_10 == "c7"),
		length(which(T15_int@meta.data$clusters_10 == "c7")),replace=FALSE))
cellOrdering<-cellOrdering

#Get gene ordering based on conserved markers (first 20 markers from each)
geneOrdering<-c()
for(i in c(4,7,1,3,5,6,10,8,2,9)){
	geneOrdering<-append(geneOrdering, sample(genesList[[i]],length(genesList[[i]]),replace=F))
}
geneOrdering<-rev(geneOrdering)


#Get module ordering (cluster for targeted DEGs)
moduleOrdering<-rep(c("#1619BA","mediumpurple3","#3281FF","tomato","orange","green3","tan","gold4","darkturquoise","saddlebrown"),each=20)
moduleOrdering<-rev(moduleOrdering)

#Now make a heatmap using heatmap3 of expression across genes and cells
mat<-sct_exp[geneOrdering,cellOrdering]

#Get side colors for population
#Note that the colors need to be in numerical order, not in the order of cellOrdering
clusters<-color.scale(as.numeric(as.factor(T15_int@meta.data[cellOrdering,]$clusters_10)),extremes=cols_clusters,color.spec="rgb")
ColSideColors<-data.frame(clusters=clusters)

#Get row colors for modules
modules<-color.scale(as.numeric(as.factor(moduleOrdering)),extremes=sort(unique(moduleOrdering)), color.spec="rgb")
RowSideColors<-data.frame(modules=modules)

mycols = colorRampPalette(c("blue","white","red"))(1000)
breakscale = c(-8,seq(-1.3,1.3, length.out=length(mycols)-1), 20)

pdf('T15_markers_heatmap_clusters.pdf',height=24,width=18)
heatmap3(mat, scale="row",Rowv=NA,Colv=NA, 
	col=mycols,ColSideColors=as.matrix(ColSideColors),RowSideColors=as.matrix(RowSideColors),
	cexRow=0.2,cexCol=1,breaks=breakscale,labCol="",useRaster=F)
dev.off()

#Export as png, since useRaster = F causes part of the heat map to be chopped off when converting to pdf later
png('T15_markers_heatmap_clusters.png',height=24,width=18,units="in",res=400)
heatmap3(mat, scale="row",Rowv=NA,Colv=NA, 
	col=mycols,
	cexRow=0.5,cexCol=1,breaks=breakscale,labCol="",useRaster=T)
dev.off()




























####################  Plot Fischer genes onto the UMAP
sct_exp<-as.matrix(GetAssayData(T15_int,assay="SCT"))

####################### SURFACE GENES
surface_fischer  = c("S100P", "PLAC8", "GSTA1", "SERPINB3", "CXCL1", "PRSS23", "FLJ32626",
                      "CYP4B1", "S100A2", "UNG2", "LOC222171", "ADH1C", "DHRS9",
                      "TSPAN1", "CEACAM5", "AGR2", "MYB", "RDH10", "RIPK4", "ANXA8L1",
                      "DNCL1", "MAL2", "TSGA2", "TUBA3", "TSPAN3", "CSTB", "NQO1", 
                      "HSPA1B", "CGI-38", "TNFRSF21", "FLJ34512", "KRT5", "ORF1-FL49",
                      "ALDH2", "DNAI2", "CLCA2", "SLC27A2", "FLJ20160", "LRP11",
                      "CCDC17", "BCMP11", "ALDH1A1", "GSTP1", "LOC165186", 
                      "CIB1", "AKR1C2", "ALDH3A1", "ADPRHL2", "IL8", "C5AR1")

#Which not in dataset?
surface_fischer[-which(surface_fischer %in% rownames(sct_exp))]

# in the genes for the surface samples, there are 15 genes missing. I looked them up
# to confirm that they are not under different names

# genes that do not have aliases/are the same names: FLJ32626, FLJ34512, FLJ20160
alt_surface = list(c("UNG2", "UNG"), c("LOC222171", "PRR15"), c("DNCL1", "DYNLL1"), 
                      c("TSGA2", "RSPH1"), c("TUBA3", "TUBA1A"),c("CGI-38", "TPPP3"),
                      c("ORF1-FL49", "CYSTM1"), c("BCMP11", "AGR3"), 
                      c("LOC165186", "TOGARAM2"), c("IL8", "CXCL8"))

old_surface_names = sapply(alt_surface, function(x)x[1])
new_surface_names = sapply(alt_surface, function(x)x[2])

surface_fischer = surface_fischer[which(! surface_fischer %in% old_surface_names)]
surface_fischer = c(surface_fischer, new_surface_names)

#Filter out the missing three genes
surface_fischer<-surface_fischer[which(surface_fischer %in% rownames(sct_exp))]






############################# SMG GENES
smg_fischer = c("PRR4", "AZGP1", "LFT", "SERPINA3", "DMBT1", "IGHG1", "IGLC2",                                 
				"IGJ", "C7", "PIP", "IGL_LOCUS", "TIMP3", "ACTA2", "COL6A3",
                "NNMT", "DMN", "LYZ", "STATH", "A2M", "LOC124220", "AQP1",
                "LHFP", "DAPK1", "HRTA1", "PMP22", "JAM3", "APIN", "MYLK", "DPYSL2",
                "CAV1", "VIM", "SOX9", "PTGIS", "KCNN4", "TRAF5", "NOSTRIN",
                "THBS2", "SPARCL1", "CDR1", "APOD", "SAA1", "SAA1/2", "TAGLN",
                "ST3GAL6", "LOH3CR2A", "PLGLB1", "PRSS11", "NDRG2", "PTRF", "IGL LOCUS")     

#Which not in dataset?
smg_fischer[-which(smg_fischer %in% rownames(sct_exp))]

# genes not in the dataset: "IGL_LOCUS" "LHFP"      "HRTA1"     "CDR1"      "PTRF"      "IGL LOCUS" "LIX1"      "LINC00312"
alt_smg = list(c("LFT", "LIX1"), c("IGJ", "JCHAIN"), c("DMN", "SYNM"), 
               c("LOC124220", "ZG16B"),c("PRSS11", "HTRA1"),
               c("APIN", "ODAM"), c("SAA1/2", "SAA1"), c("LOH3CR2A", "LINC00312"))

old_smg_names = sapply(alt_smg, function(x)x[1])
new_smg_names = sapply(alt_smg, function(x)x[2])

smg_fischer = smg_fischer[which(! smg_fischer %in% old_smg_names)]
smg_fischer = c(smg_fischer, new_smg_names)

#Filter out the missing three genes
smg_fischer<-smg_fischer[which(smg_fischer %in% rownames(sct_exp))]








######################Now make the UMAP based on mean expression across these genes
QUANTS_SMG = c(0.5,0.75, 0.85, 0.95)
QUANTS_SURFACE = c(0.5,0.25,0.15,0.05)
for (i in 1:length(QUANTS)){
	classificationDF = data.frame(isSurface = rep(NA, nrow(T15_int@meta.data)),
	                              isSMG     = rep(NA, nrow(T15_int@meta.data)),
	                              row.names = colnames(sct_exp))
	surfaceMeans       = colMeans(sct_exp[surface_fischer,])
	threshold_surface  = quantile(surfaceMeans, QUANTS_SURFACE[i])
	
	smgMeans = colMeans(sct_exp[smg_fischer,])
	threshold_smg = quantile(smgMeans, QUANTS_SMG[i])
	
	classificationDF$isSurface[surfaceMeans < threshold_surface] = 0
	classificationDF$isSurface[surfaceMeans >= threshold_surface] = 1
	
	classificationDF$isSMG[smgMeans < threshold_smg] = 0
	classificationDF$isSMG[smgMeans >= threshold_smg] = 1
	
	# red : high in submucosal genes and low in epithelial genes
	colorsDF = data.frame(smg_surface_class= rep("grey", nrow(T15_int@meta.data)), 
	                      stringsAsFactors = FALSE, row.names = colnames(sct_exp))
	
	# Color in just the SMG ones
	colorsDF$smg_surface_class[which(classificationDF$isSMG == 1 & classificationDF$isSurface == 0)] = "red"
	T15_int <- AddMetaData(T15_int, colorsDF)
	
	
	fn = paste0("PanelA_fischer_classification_exploration_q", QUANTS_SMG[i], ".pdf")
	pdf(fn,height=5,width=5.3)
	#dev.new(height=5,width=6)
	DimPlot(T15_int, reduction='umap',group.by="smg_surface_class",
	         cols=rev(c("red", "grey"))) + NoLegend()
	dev.off()
}


