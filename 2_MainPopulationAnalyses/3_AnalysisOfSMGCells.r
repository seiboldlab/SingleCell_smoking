#########################################
#############Seurat analysis#############
#########################################

library(Seurat, lib.loc="/Users/jacksonna/R-dev")
library(plotrix)
library(scales)
library(ggplot2)
library(cowplot)
library(openxlsx)
library(heatmap3)
library(gplots)
source("CommonFunctions.r")

#Download the processed Seurat R object, "Processed_invitro_seurat.Rdata" from GEO,
#obtainable from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134174
#Then load the object it contains ("T15_int") into the session.
#T15_int is a Seurat dataset that contains the processed, integrated dataset
load("Processed_invitro_seurat.Rdata")



################################### First, get dot plot of Markers

#First need to make new metadata that allocates clusters of interest
treat.col<-data.frame("SMGrelatedGroups"=T15_int@meta.data$clusters_18,row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col$SMGrelatedGroups[which(T15_int@meta.data$clusters_18 %in% c("c2"))] <-"Surface differentiating basal"
treat.col$SMGrelatedGroups[which(T15_int@meta.data$clusters_18 %in% c("c1c"))] <-"Surface proteasomal basal"
treat.col$SMGrelatedGroups[which(T15_int@meta.data$clusters_18 %in% c("c3c"))] <-"Myoepithelial"
treat.col$SMGrelatedGroups[which(T15_int@meta.data$clusters_18 %in% c("c3a"))] <-"SMG basal a"
treat.col$SMGrelatedGroups[which(T15_int@meta.data$clusters_18 %in% c("c3b"))] <-"SMG basal b"
treat.col$SMGrelatedGroups[which(T15_int@meta.data$clusters_18 %in% c("c7"))] <-"SMG secretory"
treat.col$SMGrelatedGroups[which(T15_int@meta.data$clusters_18 %in% c("c1","c4","c5a","c5b","c5c","c5d","c6a","c6b","c6c","c8a","c8b","c8c"))] <-"Surface other"

T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="SMGrelatedGroups")
T15_int@meta.data$SMGrelatedGroups<-factor(T15_int@meta.data$SMGrelatedGroups,levels=c("Surface differentiating basal",
	"Surface proteasomal basal","Myoepithelial","SMG basal a", "SMG basal b","SMG secretory","Surface other"))

markers<-c("KRT5","MKI67","IL33","KRT14","ACTA2","SOX9","CAV2","ANXA5","IFITM3")
markers<-c("KRT5","MKI67","IL33","KRT14","IFITM3","CAV1","CAV2","ANXA5","ACTN1","ACTA2","MIA","VIM","LGALS1","FOLR1","BARX2","DMBT1")

#Unite c1c and c2
markers<-c("IL33","TP63","TGFB1","IFITM1","COL7A1","TNC","LIMA1","SGK1","KRT14")
#Unite SMG
markers<-c("FOXC1","CFI","SCARA3","SOX9","RUNX3","CYP4X1","HTRA1","TOMM7","ARL4C","CTNNB1","NINJ1","SCPEP1")
#Unite c3b and c7
markers<-c("TF","BARX2","MANSC1","TMEM63A","DMBT1","CEACAM1","ARHGEF10L","CA2","HOMER2","TMEM213")
#Only up in c3a
markers<-c("ID1","G0S2","DUSP1","SPRY1","SNHG8","LOXL4","PLAT","TMEM213")
#Only up in c3b
markers<-c("PLAT","TMEM213","SHISA2","STC1")
#Only up in c3c
markers<-c("MYH11","ACTN2","ACTA2","CNN1","BGN","SFRP1")
#Unite c3a with myo
markers<-c("MMP2","COL14A1","PTK2B")
#Unite c3a with c3b
markers<-c("CREG1","DDIT3","TWIST1","CTSH","BCL11A","IL1R1")
#Final markers
markers<-c("IL33","COL7A1","TNC","IFITM3","KRT14","CAV2","CTNNB1","HTRA1","ARL4C","SOX9","SCARA3","RUNX3","FOXC1",
	"MYH11","ACTA2","CNN1","BGN","SFRP1","VIM","ACTN1","PTK2B","MMP2","COL14A1","ID1","G0S2","NUAK1",
	"CREG1","BCL11A","DDIT3","PLAT","SHISA2","STC1","BARX2","TMEM63A","CA2","CEACAM1","DMBT1","FOLR1","LTF","STATH")


##CLUSTERS
#Set desired order of clusters to ident
Idents(T15_int)<-T15_int@meta.data$SMGrelatedGroups

#pdf("DotPlot_SMGmarkers.pdf")
#pdf("DotPlot_SMGmarkers_expanded.pdf")
dev.new(width=5.8,height=14.3)
dotPlots<-DotPlot(T15_int,features=markers,plot.legend=T,col.min=-1,col.max=1,dot.scale=9.5,
	cols = c("lightgrey", "blue"),do.return=T,assay="SCT")
dotPlots + theme(axis.text.x=element_text(angle = -45, hjust = 0, size = 12)) + coord_flip()

















########################################## Do DE to determine differences between surface and SMG secretory
############## First, get clusters into the Metadata if not already there
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





#Now do DE for each mature secretory group (c5M and c7) against all non-secretory cells (c5 and c7)
Idents(T15_int) <- T15_int@meta.data$clusters_16a
clusterList<-sort(as.character(unique(Idents(T15_int))))
markers.list <- list()
markers.list[[length(markers.list) + 1]] <- FindMarkers(T15_int, ident.1="c5M", ident.2=clusterList[-c(8,9,13)],
	min.pct=0.1, logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "LR", latent.vars = "smoke")
markers.list[[length(markers.list) + 1]] <- FindMarkers(T15_int, ident.1="c7", ident.2=clusterList[-c(8,9,13)],
	min.pct=0.1, logfc.threshold=0.25, only.pos=T, assay="SCT", test.use = "LR", latent.vars = "smoke")

#Combine tables
comparisonVec<-c("c5M_v_nonsec","c7_v_nonsec")

#Get vector of clusters for calculating pct values from RNA
pctIndexList<-list()
pctIndexList[[length(pctIndexList) + 1]]<-list("c5M",clusterList[-c(8,9,13)])
pctIndexList[[length(pctIndexList) + 1]]<-list("c7",clusterList[-c(8,9,13)])

#Get RNA counts and cluster assignments
counts<-as.matrix(GetAssayData(T15_int,assay="RNA",slot="data"))
clusterAssignments<-T15_int@meta.data$clusters_16a

#Run OrganizeDF_list
SurfSMG_DEGs<-OrganizeDF_list(datasetList=markers.list,comparisonVec=comparisonVec,pctIndexList=pctIndexList,
	counts=counts,clusterAssignments=clusterAssignments)

#####Add in unique and common gene lists between these
#First, get common genes between the two
c5M_DEGs<-SurfSMG_DEGs[which(SurfSMG_DEGs$comparison == "c5M_v_nonsec" & SurfSMG_DEGs$p_adj_FDR < 0.05),]
c7_DEGs<-SurfSMG_DEGs[which(SurfSMG_DEGs$comparison == "c7_v_nonsec" & SurfSMG_DEGs$p_adj_FDR < 0.05),]
c5M_DEGs_common<-c5M_DEGs[which(c5M_DEGs$gene %in% c7_DEGs$gene),]
c7_DEGs_common<-c7_DEGs[which(c7_DEGs$gene %in% c5M_DEGs$gene),]
c5M_DEGs_common<-c5M_DEGs_common[order(c5M_DEGs_common$gene),]
c7_DEGs_common<-c7_DEGs_common[order(c7_DEGs_common$gene),]
all(c5M_DEGs_common$gene == c7_DEGs_common$gene)
c5Mc7_v_nonsec_common<-data.frame("gene.c5M"=c5M_DEGs_common$gene,"avg_logFC.c5M"=c5M_DEGs_common$avg_logFC,"p_adj_FDR.c5M"=c5M_DEGs_common$p_adj_FDR,
	"pct.1.c5M"=c5M_DEGs_common$pct.1,"pct.2.c5M"=c5M_DEGs_common$pct.2,"gene.c7"=c7_DEGs_common$gene,"avg_logFC.c7"=c7_DEGs_common$avg_logFC,
	"p_adj_FDR.c7"=c7_DEGs_common$p_adj_FDR,"pct.1.c7"=c7_DEGs_common$pct.1,"pct.2.c7"=c7_DEGs_common$pct.2)

#Get unique genes
c5M_DEGs_unique<-c5M_DEGs[which(!(c5M_DEGs$gene %in% c7_DEGs$gene)),]
c7_DEGs_unique<-c7_DEGs[which(!(c7_DEGs$gene %in% c5M_DEGs$gene)),]

#Write to excel
write.xlsx(list(c5M_DEGs,c7_DEGs,c5Mc7_v_nonsec_common,c5M_DEGs_unique,c7_DEGs_unique),file="DEG_Tables/c5M_c7_DEGs.xlsx",
	sheetName=c("c5M_v_nonsec","c7_v_nonsec","c5M_c7_v_nonsec_COMMON","c5M_v_nonsec_UNIQUE","c7_v_nonsec_UNIQUE"))

#Make version to export for Enrichr
c5M_c7_DEGs_forEnrichr<-rbind(c5M_DEGs,c7_DEGs,data.frame(c5M_DEGs_common[,-9],"comparison"="c5M_c7_v_nonsec_COMMON"),
	data.frame(c5M_DEGs_unique[,-9],"comparison"="c5M_v_nonsec_UNIQUE"),data.frame(c7_DEGs_unique[,-9],"comparison"="c7_v_nonsec_UNIQUE"))

#Write all to single txt file
write.table(c5M_c7_DEGs_forEnrichr,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,
	file="Enrichr/c5M_c7_DEGs_forEnrichr.txt")

#Do enrichments
DEG_table<-read.table("Enrichr/c5M_c7_DEGs_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"c5M_c7_DEGs"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)




















############### Do heat map of these genes (SCT expression; Supplementary Figure 3b)

#Bring in c5M_c7 genes
c5M_c7_DEGs_forEnrichr<-read.table("Enrichr/c5M_c7_DEGs_forEnrichr.txt",sep="\t",header=T,stringsAsFactors=F)

c5up<-c5M_c7_DEGs_forEnrichr$gene[which(c5M_c7_DEGs_forEnrichr$comparison == "c5M_v_nonsec_UNIQUE")]
c7up<-c5M_c7_DEGs_forEnrichr$gene[which(c5M_c7_DEGs_forEnrichr$comparison == "c7_v_nonsec_UNIQUE")]
c5_7up<-c5M_c7_DEGs_forEnrichr[which(c5M_c7_DEGs_forEnrichr$comparison == "c5M_c7_v_nonsec_COMMON"),]
c5_7up<-c5_7up[order(c5_7up$p_adj_FDR),]$gene

#Randomize order of genes
c5up<-sample(c5up,length(c5up),replace=FALSE)
c7up<-sample(c7up,length(c7up),replace=FALSE)
c5_7up<-sample(c5_7up,length(c5_7up),replace=FALSE)

cols_clusters<-c("orange","saddlebrown","brown")
cols_clusters_cells<-c("orange","saddlebrown")
cols_subclusters_cells<-c("tomato","darkgoldenrod","darkorange2","brown","black","saddlebrown")
cols_donor<-c("midnightblue","black","red","pink","tan","greenyellow","orangered4","yellow","lightblue",
	"grey","darkgreen","blue","purple4","palegreen","slateblue")

#Subset data
T15_c5_c7<-subset(T15_int,cells=rownames(T15_int@meta.data)[which(T15_int@meta.data$clusters_16a == "c5M" | T15_int@meta.data$clusters_16a == "c7")])

#First, need to isolate those genes we want and cells we want in the order that we want them
###############
#cellOrdering<-c(
#	which(T15_c5_c7@meta.data$clusters_18 == "c4"),
#	which(T15_c5_c7@meta.data$clusters_18 == "c5a"),
#	which(T15_c5_c7@meta.data$clusters_18 == "c5b"),
#	which(T15_c5_c7@meta.data$clusters_18 == "c5c"),
#	which(T15_c5_c7@meta.data$clusters_18 == "c5d"),
#	which(T15_c5_c7@meta.data$clusters_18 == "c7"))

#cellOrdering<-c(
#	which(T15_c5_c7@meta.data$clusters_16a == "c5M"),
#	which(T15_c5_c7@meta.data$clusters_16a == "c7"))

#Randomize (to minimize donor differences)
cellOrdering<-c(
	sample(which(T15_c5_c7@meta.data$clusters_16a == "c5M"),
		length(which(T15_c5_c7@meta.data$clusters_16a == "c5M")),replace=FALSE),
	sample(which(T15_c5_c7@meta.data$clusters_16a == "c7"),
		length(which(T15_c5_c7@meta.data$clusters_16a == "c7")),replace=FALSE))

#Get gene ordering based on DEGs
geneOrdering<-rev(c(c5up,c7up,c5_7up))

#Get module ordering (cluster for targeted DEGs)
moduleOrdering<-c(rep("orange",length(c5up)),rep("saddlebrown",length(c7up)),rep("brown",length(c5_7up)))
moduleOrdering<-rev(moduleOrdering)

#Now make a heatmap using heatmap3 of expression across genes and cells
mat<-as.matrix(GetAssayData(T15_c5_c7,assay="SCT"))[geneOrdering,cellOrdering]

#Get side colors for population
#Note that the colors need to be in numerical order, not in the order of cellOrdering
clusters<-color.scale(as.numeric(as.factor(T15_c5_c7@meta.data[cellOrdering,]$clusters_16a)),extremes=cols_clusters_cells,color.spec="rgb")
subclusters<-color.scale(as.numeric(as.factor(T15_c5_c7@meta.data[cellOrdering,]$clusters_18)),extremes=cols_subclusters_cells,color.spec="rgb")
donors<-color.scale(as.numeric(as.factor(T15_c5_c7@meta.data[cellOrdering,]$donor)),extremes=cols_donor,color.spec="rgb")
ColSideColors<-data.frame(clusters=clusters,subclusters=subclusters,donors=donors)

#Get row colors for modules
modules<-color.scale(as.numeric(as.factor(moduleOrdering)),extremes=sort(unique(moduleOrdering)), color.spec="rgb")
RowSideColors<-data.frame(modules=modules)

mycols = colorRampPalette(c("blue","white","red"))(1000)
breakscale = c(-8,seq(-1.2,1.2, length.out=length(mycols)-1), 8)

pdf('T15_c5M_c7_markersSCT_heatmap.pdf',height=24,width=18)
heatmap3(mat, scale="row",Rowv=NA,Colv=NA, 
	col=mycols,ColSideColors=as.matrix(ColSideColors),RowSideColors=as.matrix(RowSideColors),
	cexRow=0.2,cexCol=1,breaks=breakscale,labCol="",useRaster=T)
dev.off()













############################### Make box plots for supplement for each of the three gene sets in the heat map

sct_exp<-as.matrix(GetAssayData(T15_int,assay="SCT"))

#Now get mean expression and make box plots
mean_exp<-CalculateMeanExpression(sct_exp,list(c5up,c7up,c5_7up),method="geometric")
T15_int<-AddMetaData(T15_int,metadata=mean_exp,col.name=c("c5up","c7up","c5_7up"))

#Divide into three groups
treat.col<-data.frame("SMGrelatedGroups"=T15_int@meta.data$clusters_16a,row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col$SMGrelatedGroups[which(T15_int@meta.data$clusters_16a %in% c("c5M"))] <-"Surface"
treat.col$SMGrelatedGroups[which(T15_int@meta.data$clusters_16a %in% c("c7"))] <-"SMG"
treat.col$SMGrelatedGroups[which(T15_int@meta.data$clusters_16a %in% unique(T15_int@meta.data$clusters_16a)[-c(5,12)])] <-"Other"

T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="SMGrelatedGroups")
T15_int@meta.data$SMGrelatedGroups<-factor(T15_int@meta.data$SMGrelatedGroups,levels=c("Surface","SMG","Other"))

#Now make box plots
#pdf('T15_c5M_c7_markersSCT_boxplots_surfaceGenes.pdf',width=2, height=3)
dev.new(height=3,width=1.8)
par(bty="n")
boxplot(T15_int@meta.data$c5up~T15_int@meta.data$SMGrelatedGroups,las=1,ylab="Mean normalized expression",
	main="c5up",cex.main=0.3,col=c("yellow2","saddlebrown","gray"),medcol="white",cex=0.4,names=c("Surface","SMG","Other"),outlwd=0.5,outcex=0.3)

#pdf('T15_c5M_c7_markersSCT_boxplots_smgGenes.pdf',width=2, height=3)
dev.new(height=3,width=1.8)
par(bty="n")
boxplot(T15_int@meta.data$c7up~T15_int@meta.data$SMGrelatedGroups,las=1,ylab="Mean normalized expression",
	main="c7up",cex.main=0.3,col=c("yellow2","saddlebrown","gray"),medcol="white",cex=0.4,names=c("Surface","SMG","Other"),outlwd=0.5,outcex=0.3)
	
#pdf('T15_c5M_c7_markersSCT_boxplots_bothGenes.pdf',width=2, height=3)
dev.new(height=3,width=1.8)
par(bty="n")
boxplot(T15_int@meta.data$c5_7up~T15_int@meta.data$SMGrelatedGroups,las=1,ylab="Mean normalized expression",
	main="c5_7up",cex.main=0.3,col=c("yellow2","saddlebrown","gray"),medcol="white",cex=0.4,names=c("Surface","SMG","Other"),outlwd=0.5,outcex=0.3)


















############### Do heat map of these genes (use subsetted genes for main Figure 3 - SCT data)

#Bring in enrichments and DEGs isolate the c5M unique enrichments
enrich_culled<-readExcelTabs("Enrichr/c5M_c7_DEGs_Enrichr.output.xlsx",returnSingleTable=F)
c5M_c7_DEGs_forEnrichr<-read.table("Enrichr/c5M_c7_DEGs_forEnrichr.txt",header=T,stringsAsFactors=F)
enrich_culled_c5M<-enrich_culled[[4]]


####Want gene sets for the following enrichments for the c5M DEGs
#Protein processing in endoplasmic reticulum
#Transport to the Golgi and subsequent modification
#O−linked glycosylation of mucins
#Asparagine N−linked glycosylation
enrichedGeneLists<-lapply(strsplit(enrich_culled_c5M$Genes,";"),function(x)x)[c(9,16,8,3)]
names(enrichedGeneLists)<-c("Protein processing in endoplasmic reticulum",
	"Golgi vesicle transport","O−linked glycosylation of mucins",
	"Asparagine N−linked glycosylation")

#Select genes
c5_enrich_genes1<-c5M_c7_DEGs_forEnrichr[which(c5M_c7_DEGs_forEnrichr$comparison == "c5M_v_nonsec_UNIQUE" & c5M_c7_DEGs_forEnrichr$gene %in% 
	enrichedGeneLists$"Protein processing in endoplasmic reticulum"),][c(1,3,4,9,10),]$gene
c5_enrich_genes2<-c5M_c7_DEGs_forEnrichr[which(c5M_c7_DEGs_forEnrichr$comparison == "c5M_v_nonsec_UNIQUE" & c5M_c7_DEGs_forEnrichr$gene %in% 
	enrichedGeneLists$"Golgi vesicle transport"),][c(1,2,4,5,9),]$gene
c5_enrich_genes3<-c5M_c7_DEGs_forEnrichr[which(c5M_c7_DEGs_forEnrichr$comparison == "c5M_v_nonsec_UNIQUE" & c5M_c7_DEGs_forEnrichr$gene %in% 
	enrichedGeneLists$"O−linked glycosylation of mucins"),][c(1,2,4,7,11),]$gene
c5_enrich_genes4<-c5M_c7_DEGs_forEnrichr[which(c5M_c7_DEGs_forEnrichr$comparison == "c5M_v_nonsec_UNIQUE" & c5M_c7_DEGs_forEnrichr$gene %in% 
	enrichedGeneLists$"Asparagine N−linked glycosylation"),][c(2,3,4,7,12),]$gene

###Now just get the top 20 DEGs for the c7 unique
c7up<-c5M_c7_DEGs_forEnrichr$gene[which(c5M_c7_DEGs_forEnrichr$comparison == "c7_v_nonsec_UNIQUE")]
c7_genes<-c7up[c(1:9,11:14,18:24)]
#Randomize order
c7_genes<-sample(c7_genes,length(c7_genes),replace=FALSE)


###Finally, plot TFs
### Get unique DEGs that overlap with both the IPA upstream regulators list and with the integrated genes available
TFs<-read.table("Gene_lists/database_genes_IPA.txt",header=T,sep="\t",stringsAsFactors=F)
TFs<-TFs$Gene[which(TFs$Type_IPA=="transcription regulator")]
c5up<-c5M_c7_DEGs_forEnrichr$gene[which(c5M_c7_DEGs_forEnrichr$comparison == "c5M_v_nonsec_UNIQUE")]
c7up<-c5M_c7_DEGs_forEnrichr$gene[which(c5M_c7_DEGs_forEnrichr$comparison == "c7_v_nonsec_UNIQUE")]
c5up_TFs<-c5up[which(c5up %in% TFs)]
c7up_TFs<-c7up[which(c7up %in% TFs)]
TFs_genes<-c("SPDEF","MAGED1","CREB3L1","CREB3L4","PRRX2","BARX2","NFIB","FOXC1","SOX9","TFCP2L1")


#Get colors for bars
cols_clusters<-c("orange","saddlebrown","brown")
cols_clusters_cells<-c("orange","saddlebrown")
cols_subclusters_cells<-c("tomato","darkgoldenrod","darkorange2","brown","black","saddlebrown")
cols_donor<-c("midnightblue","black","red","pink","tan","greenyellow","orangered4","yellow","lightblue",
	"grey","darkgreen","blue","purple4","palegreen","slateblue")

#Subset data
T15_c5_c7<-subset(T15_int,cells=rownames(T15_int@meta.data)[which(T15_int@meta.data$clusters_16a == "c5M" | T15_int@meta.data$clusters_16a == "c7")])

#First, need to isolate those cells we want (in a randomized order to minimize donor differences)
###############
cellOrdering<-c(
	sample(which(T15_c5_c7@meta.data$clusters_16a == "c5M"),
		length(which(T15_c5_c7@meta.data$clusters_16a == "c5M")),replace=FALSE),
	sample(which(T15_c5_c7@meta.data$clusters_16a == "c7"),
		length(which(T15_c5_c7@meta.data$clusters_16a == "c7")),replace=FALSE))

##Order cells in donor-clumps to see batch effects
#cellOrdering<-c(
#	which(T15_c5_c7@meta.data$clusters_16a == "c5M"),
#	which(T15_c5_c7@meta.data$clusters_16a == "c7"))

#Get gene ordering based on DEGs
geneOrdering<-rev(c(c5_enrich_genes1,c5_enrich_genes2,c5_enrich_genes3,c5_enrich_genes4,c7_genes,TFs_genes))

#Get module ordering (cluster for targeted DEGs)
moduleOrdering<-c(rep("orange",length(c(c5_enrich_genes1,c5_enrich_genes2,c5_enrich_genes3,c5_enrich_genes4))),
	rep("saddlebrown",length(c7_genes)),rep("brown",length(TFs_genes)))
moduleOrdering<-rev(moduleOrdering)

#Now make a heatmap using heatmap3 of expression across genes and cells
mat<-as.matrix(GetAssayData(T15_c5_c7,assay="SCT",slot="data"))[geneOrdering,cellOrdering]

#Get side colors for population
#Note that the colors need to be in numerical order, not in the order of cellOrdering
clusters<-color.scale(as.numeric(as.factor(T15_c5_c7@meta.data[cellOrdering,]$clusters_16a)),extremes=cols_clusters_cells,color.spec="rgb")
subclusters<-color.scale(as.numeric(as.factor(T15_c5_c7@meta.data[cellOrdering,]$clusters_18)),extremes=cols_subclusters_cells,color.spec="rgb")
donors<-color.scale(as.numeric(as.factor(T15_c5_c7@meta.data[cellOrdering,]$donor)),extremes=cols_donor,color.spec="rgb")
ColSideColors<-data.frame(clusters=clusters,subclusters=subclusters,donors=donors)

#Get row colors for modules
modules<-color.scale(as.numeric(as.factor(moduleOrdering)),extremes=sort(unique(moduleOrdering)), color.spec="rgb")
RowSideColors<-data.frame(modules=modules)

mycols = colorRampPalette(c("blue","white","red"))(1000)
breakscale = c(-8,seq(-1.2,1.2, length.out=length(mycols)-1),8)

pdf('T15_c5M_c7_markersSCT_culled_heatmap.pdf',height=24,width=18)
heatmap3(mat, scale="row",Rowv=NA,Colv=NA, 
	col=mycols,ColSideColors=as.matrix(ColSideColors),RowSideColors=as.matrix(RowSideColors),
	cexRow=2,cexCol=1,breaks=breakscale,labCol="",useRaster=F)
dev.off()


























#########################Now, make a bar graph of cell proportions of MUCs, breaking out smoking status
#Subset c7 cells
T15_smgSec<-subset(T15_int, cells=rownames(T15_int@meta.data)[which(T15_int@meta.data$clusters_10 == "c7")])

#Define clusters to plot
T15_sct<-as.matrix(GetAssayData(T15_smgSec,assay="SCT"))

donorList_ns<-unique(T15_smgSec@meta.data$donor[which(T15_smgSec@meta.data$smoke_noT89 == "never")])
donorList_hs<-unique(T15_smgSec@meta.data$donor[which(T15_smgSec@meta.data$smoke_noT89 == "heavy")])
pvalList<-list()

#Donor colors
colorVec<-c("pink","tan","greenyellow","yellow","lightblue","grey","midnightblue","black","orangered4","darkgreen","purple4","slateblue")

#Do boxplot or do barplot
do.boxplot <- F
#Make plot
#pdf("T15_smgSec_MUC_smoking_barplots.pdf")
#pdf("T15_smgSec_MUC_smoking_boxplots.pdf")
dev.new(height=3,width=3)
for(i in 1:1){
	currData_ns<-list()
	currData_hs<-list()
	for(j in 1:length(donorList_ns)){
		currData_ns[[length(currData_ns) + 1]]<-T15_sct[c("MUC5AC","MUC5B"),which(T15_smgSec@meta.data$donor == donorList_ns[j])]
		currData_hs[[length(currData_hs) + 1]]<-T15_sct[c("MUC5AC","MUC5B"),which(T15_smgSec@meta.data$donor == donorList_hs[j])]
	}
	
	bp_mat<-data.frame(matrix(NA,nrow=6,ncol=8))
	
	#Proportion of cells expressing MUC5AC
	bp_mat[,1]<-sapply(currData_ns,function(x)length(which(x["MUC5AC",] != 0 & x["MUC5B",] == 0)) / ncol(x))
	bp_mat[,2]<-sapply(currData_hs,function(x)length(which(x["MUC5AC",] != 0 & x["MUC5B",] == 0)) / ncol(x))
	#Proportion of cells expressing MUC5B
	bp_mat[,3]<-sapply(currData_ns,function(x)length(which(x["MUC5AC",] == 0 & x["MUC5B",] != 0)) / ncol(x))
	bp_mat[,4]<-sapply(currData_hs,function(x)length(which(x["MUC5AC",] == 0 & x["MUC5B",] != 0)) / ncol(x))
	#Proportion of cells expressing MUC5AC AND MUC5B
	bp_mat[,5]<-sapply(currData_ns,function(x)length(which(x["MUC5AC",] != 0 & x["MUC5B",] != 0)) / ncol(x))
	bp_mat[,6]<-sapply(currData_hs,function(x)length(which(x["MUC5AC",] != 0 & x["MUC5B",] != 0)) / ncol(x))
	#Proportion of cells expressing neither mucus
	bp_mat[,7]<-sapply(currData_ns,function(x)length(which(x["MUC5AC",] == 0 & x["MUC5B",] == 0)) / ncol(x))
	bp_mat[,8]<-sapply(currData_hs,function(x)length(which(x["MUC5AC",] == 0 & x["MUC5B",] == 0)) / ncol(x))

	#Names
	names(bp_mat)<-c("MUC5AC_ns","MUC5AC_hs","MUC5B_ns","MUC5B_hs","Both_ns","Both_hs","Neither_ns","Neither_hs")
	
	#Get standard error
	sderror_vec<-apply(bp_mat,2,function(x)sd(x)/sqrt(length(x)))
	
	#Make bar plot
	if(do.boxplot == F){
		barplot2(apply(bp_mat,2,mean),las=1,col=c("black","red"),ylim=c(0,0.7),names.arg="",ylab="Proportion of cells",
			ci.l = apply(bp_mat,2,mean) - sderror_vec, ci.u = apply(bp_mat,2,mean) + sderror_vec, plot.ci=T, ci.color="grey45",)
		text(x=c(1.5,3.8,6.1,8.6),y=par()$usr[3]-0.08*(par()$usr[4]-par()$usr[3]),labels=c("MUC5AC","MUC5B","Both","Neither"),srt=45,adj=1,xpd=T)
	}
	
	#Make box plots
	if(do.boxplot == T){
		currData<-melt(bp_mat)
		par(bty="l")
		boxplot(currData$value~currData$variable,las=1,names=rep("",8),
			border=c("black","red"),medcol=c("black","red"),yaxt="n",
			ylim=c(min(currData$value),max(currData$value)*1.1))
		points(x=c(rep(1:8,each=6)),y=currData$value,col=colorVec,pch=16,cex=1)
		text(x=seq(1.5,8.5,by=2),y=par()$usr[3]-0.12*(par()$usr[4]-par()$usr[3]),labels=c("MUC5AC","MUC5B","Both","Neither"),
			srt=45,adj=1,xpd=T,cex=1)
		axis(2,at=round(seq(from=min(currData$value),to=max(currData$value),length.out=5),digits=6),las=1,
			labels=formatC(seq(from=min(currData$value),to=max(currData$value),length.out=5),format="f",digits=2))
		title(ylab="Proportion of cells", line=3.3, cex.lab=1)
	}
	
	#Now test the significance in differences between groups (p-values are "0.03004" "0.18882" "0.81118" "0.59091")
	count<-1
	for(j in 1:(ncol(bp_mat) / 2)){
		currData<-bp_mat[,count:(count+1)]
		if(j==1){
			pvalList[[i]]<-formatC(wilcox.test(currData[,1],currData[,2],alternative="less")$p.value,digits=4)
		}else{
			pvalList[[i]]<-append(pvalList[[i]],formatC(wilcox.test(currData[,1],currData[,2],alternative="greater")$p.value,digits=5))
		}
		count<-count + 2
	}
}	















############################Make pie chart of MUC panels
proportionsList<-list()
for(i in 1:1){
	currData<-as.matrix(GetAssayData(T15_smgSec,assay="SCT"))[c("MUC5AC","MUC5B"),]
	proportions<-c()
	proportions<-append(proportions,length(which(currData["MUC5AC",] != 0 & currData["MUC5B",] == 0)) / ncol(currData))
	proportions<-append(proportions,length(which(currData["MUC5AC",] == 0 & currData["MUC5B",] != 0)) / ncol(currData))
	proportions<-append(proportions,length(which(currData["MUC5AC",] != 0 & currData["MUC5B",] != 0)) / ncol(currData))
	proportions<-append(proportions,length(which(currData["MUC5AC",] == 0 & currData["MUC5B",] == 0)) / ncol(currData))
	names(proportions)<-c("MUC5AC","MUC5B","Both","Neither")
	proportionsList[[length(proportionsList) + 1]]<-proportions
}

#Plot pie
#pdf("T15_smgSec_MUCs_piechart.pdf")
for(i in 1:length(proportionsList)){
	pct <- round(proportionsList[[i]],3)*100
	lbls <- paste(pct,"%",sep="") # ad % to labels 
	colors<-c("cornflowerblue","forestgreen","darkorange","black")
	pie(proportionsList[[i]],labels = names(proportionsList[[i]]), col=colors, main = names(proportionsList)[[i]])
}

















############################## Make MUC5AC and MUC5B box plots for surface and SMG secretory cell groups
#First, divide cells into surface, smg, and other
treat.col<-data.frame("surfaceSMG"=T15_int@meta.data$clusters_16a,row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col$surfaceSMG[which(T15_int@meta.data$clusters_16a %in% c("c5M"))] <-"Surface secretory"
treat.col$surfaceSMG[which(T15_int@meta.data$clusters_16a %in% c("c7"))] <-"SMG secretory"
treat.col$surfaceSMG[which(T15_int@meta.data$clusters_16a %in% grep("c7|c5",unique(T15_int@meta.data$clusters_16a),value=T,invert=T))] <-"Other"
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="surfaceSMG")
T15_int@meta.data$surfaceSMG<-factor(T15_int@meta.data$surfaceSMG,levels=c("Surface secretory","SMG secretory","Other"))

#Define matrix
T15_sct<-as.matrix(GetAssayData(T15_int,assay="RNA"))

#Now make box plots
#pdf("SurfaceVsSMG_boxplot_MUC5AC.pdf")
dev.new(width=2,height=3)
par(bty="l")
boxplot(T15_sct["MUC5AC",]~T15_int@meta.data$surfaceSMG,pch=16,las=1,col=c("brown","saddlebrown","grey"),medcol="white",cex=0.4)

#pdf("SurfaceVsSMG_boxplot_MUC5B.pdf")
dev.new(width=2,height=3)
par(bty="l")
boxplot(T15_sct["MUC5B",]~T15_int@meta.data$surfaceSMG,pch=16,las=1,col=c("brown","saddlebrown","grey"),medcol="white",cex=0.4)
	
	
	



######## Get fold changes between the MUCs in surface versus SMG secretory
#Now do DE for each mature secretory group (c5M and c7) against all non-secretory cells (c5 and c7)
Idents(T15_int) <- T15_int@meta.data$clusters_16a
markers.list <- list()
markers.list[[length(markers.list) + 1]] <- FindMarkers(T15_int, ident.1="c5M", ident.2="c7",
	min.pct=0.1, logfc.threshold=0.25, only.pos=F, assay="SCT", test.use = "LR", latent.vars = "smoke")

#Combine tables
comparisonVec<-c("c5M_v_c7")
#Run OrganizeDF_list
SurfSMG_DEGs<-OrganizeDF_list(datasetList=markers.list,comparisonVec=comparisonVec)
	
SurfSMG_DEGs[which(SurfSMG_DEGs$gene == "MUC5AC" | 	SurfSMG_DEGs$gene == "MUC5B"),]
exp(3.2367) #25.4496
exp(0.3612) #-1.4351






























############################### Make enrichment box plots for supplement for shared surface and SMG markers

sct_exp<-as.matrix(GetAssayData(T15_int,assay="SCT"))

#Define gene sets
#First, Transmembrane transport of small molecules_Homo sapiens_R-HSA-382551 (based on shared genes up in both c7 and c5M compared to non-sec cells
transmembrane<-c("SLC12A2","SLC34A2","SLC26A2","SLC44A4","SLC35C1","ATP2A3","ATP10B","SLC5A1","ATP12A","CP","SLC4A4","SLC5A8","AZGP1","LCN2","SLC39A7","SLC26A4","SLC16A3","TRPM4")
#Next get "Mucosal defense" genes (randomly selected by me)
defense<-c("BPIFB1","BPIFA1","C3","CP","CD55","CYP2F1","GSTA1","LCN2","LYN","LYZ","PIGR","SLPI","WFDC2")

#Now get mean expression and make box plots
mean_exp<-CalculateMeanExpression(sct_exp,list(transmembrane,defense),method="geometric")
T15_int<-AddMetaData(T15_int,metadata=mean_exp,col.name=c("transmembrane","defense"))

#Divide into three groups
treat.col<-data.frame("SMGrelatedGroups"=T15_int@meta.data$clusters_16a,row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col$SMGrelatedGroups[which(T15_int@meta.data$clusters_16a %in% c("c5M"))] <-"Surface"
treat.col$SMGrelatedGroups[which(T15_int@meta.data$clusters_16a %in% c("c7"))] <-"SMG"
treat.col$SMGrelatedGroups[which(T15_int@meta.data$clusters_16a %in% unique(T15_int@meta.data$clusters_16a)[-c(5,12)])] <-"Other"

T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="SMGrelatedGroups")
T15_int@meta.data$SMGrelatedGroups<-factor(T15_int@meta.data$SMGrelatedGroups,levels=c("Surface","SMG","Other"))

#Now make box plots
#pdf('T15_surfaceVSsmg_transmembraneTerm_boxplot.pdf',width=2, height=3)
dev.new(height=3,width=1.8)
par(bty="n")
boxplot(T15_int@meta.data$transmembrane~T15_int@meta.data$SMGrelatedGroups,las=1,ylab="Mean normalized expression",
	main="Transmembrane",cex.main=0.3,col=c("yellow2","saddlebrown","gray"),medcol="white",cex=0.4,names=c("Surface","SMG","Other"),outlwd=0.5,outcex=0.3)

#pdf('T15_surfaceVSsmg_defenseGenes_boxplot.pdf',width=2, height=3)
dev.new(height=3,width=1.8)
par(bty="n")
boxplot(T15_int@meta.data$defense~T15_int@meta.data$SMGrelatedGroups,las=1,ylab="Mean normalized expression",
	main="defense",cex.main=0.3,col=c("yellow2","saddlebrown","gray"),medcol="white",cex=0.4,names=c("Surface","SMG","Other"),outlwd=0.5,outcex=0.3)
	
	mtext(paste("p = ",formatC(wilcox.test(T15_int@meta.data$transmembrane[which(T15_int@meta.data$SMGrelatedGroups=="Surface")],
	T15_int@meta.data$transmembrane[which(T15_int@meta.data$SMGrelatedGroups=="Other")],alternative="greater")$p.value,
	format="e",digits=2),sep=""),side=3,line=0.2,at=2.5,adj=1,cex=1)
dev.off()


































###################################### Make dot plots of unique smoking response genes in SMG secretory cells

#Genes to plot
up_unique<-c("HLA-F","HLA-DQA2","IFI44L","IL6","B2M","PIP","SCGB1A1","TM7SF2","STATH","WBP2","DNAJB9")
down_unique<-c("MUC5B","FCGBP","BPIFB1","S100A8","S100A9","HP","SLPI","GSTA2")

#Subset dataset to remove light and excluded donors
T15_smgSmoke<-subset(T15_int,cells=rownames(T15_int@meta.data)[grep("heavy|never",T15_int@meta.data$clusters16a_smoke)])

#Define new metadata that divides SMG-nonsmoke, SMG-smoke, other-nonsmoke, and other-smoke
treat.col<-data.frame("SMG_smoke"=T15_smgSmoke@meta.data$clusters16a_smoke,row.names=rownames(T15_smgSmoke@meta.data),stringsAsFactors=F)
treat.col$SMG_smoke[which(T15_smgSmoke@meta.data$clusters16a_smoke %in% c("c7_never"))] <-"SMG nonsmoker"
treat.col$SMG_smoke[which(T15_smgSmoke@meta.data$clusters16a_smoke %in% c("c7_heavy"))] <-"SMG smoker"
treat.col$SMG_smoke[which(T15_smgSmoke@meta.data$clusters16a_smoke %in% grep("never",unique(T15_smgSmoke@meta.data$clusters16a_smoke),value=T)[-14])] <-"Other nonsmoker"
treat.col$SMG_smoke[which(T15_smgSmoke@meta.data$clusters16a_smoke %in% grep("heavy",unique(T15_smgSmoke@meta.data$clusters16a_smoke),value=T)[-12])] <-"Other smoker"
T15_smgSmoke<-AddMetaData(T15_smgSmoke,metadata=treat.col,col.name="SMG_smoke")
T15_smgSmoke@meta.data$SMG_smoke<-factor(T15_smgSmoke@meta.data$SMG_smoke,levels=c("SMG nonsmoker","SMG smoker","Other nonsmoker","Other smoker"))

#Set desired clusters to ident
Idents(T15_smgSmoke)<-T15_smgSmoke@meta.data$SMG_smoke

#pdf("T15_SMGsecSmokeGenes_dotPlot_up.pdf")
dev.new(width=4.9,height=4.9)
dotPlots<-DotPlot(T15_smgSmoke,features=up_unique,plot.legend=T,col.min=-1,col.max=1,dot.scale=9.5,
	cols = c("lightgrey", "blue"),do.return=T,assay="SCT")
dotPlots + theme(axis.text.x=element_text(angle = -45, hjust = 0, size = 12)) + coord_flip()

#pdf("T15_SMGsecSmokeGenes_dotPlot_down.pdf")
dev.new(width=4.9,height=3.8)
dotPlots<-DotPlot(T15_smgSmoke,features=down_unique,plot.legend=T,col.min=-1,col.max=1,dot.scale=9.5,
	cols = c("lightgrey", "blue"),do.return=T,assay="SCT")
dotPlots + theme(axis.text.x=element_text(angle = -45, hjust = 0, size = 12)) + coord_flip()



















	