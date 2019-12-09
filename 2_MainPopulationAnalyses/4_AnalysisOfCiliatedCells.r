#########################################
#############Seurat analysis#############
#########################################

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




################################################Create box plots in Figure 5e
#Read in superfine markers
sf_markers<-read.table("Gene_lists/Superfine_biomarkers.txt",sep="\t",header=T,stringsAsFactors=F)

#Read in DEG list (but need to replace with with c5 markers when available)
#T0_markers<-read.table("DEG_Tables/T15_1againstAll_DEGs.txt",sep="\t",header=T,stringsAsFactors=F) #DEGs
T0_markers<-read.table("ConservedBiomarkers/Biomarkers_main/Conserved_DEGs.txt",sep="\t",header=T,stringsAsFactors=F) #Markers

#Pull out markers to plot
markerList<-list("O"=sf_markers$gene[which(sf_markers$comparison == "08_FOXN4 early ciliating")],
	"N"=sf_markers$gene[which(sf_markers$comparison == "09_Late ciliating")],
	"M"=sf_markers$gene[which(sf_markers$comparison == "10_Mature ciliated" )],
	"MucSec"=T0_markers$gene[which(T0_markers$comparison == "c5")])

#Now get mean expression of these markers across all the new T15 clusters and subclusters
mean_exp_df<-CalculateMeanExpression(dataset=as.matrix(GetAssayData(T15_int,assay="SCT")),genesList=markerList,method="geometric")

#Add these markers to metadata
T15_int<-AddMetaData(T15_int,mean_exp_df)

#How do these look in feature plots?
DefaultAssay(T15_int) <- "SCT"
pt1<-FeaturePlot(T15_int,c("O"),cols=c("gray", "blue"), pt.size=0.1, min.cutoff="q5", max.cutoff="q95") + NoAxes() + NoLegend()
pt2<-FeaturePlot(T15_int,c("N"),cols=c("gray", "blue"), pt.size=0.1, min.cutoff="q5", max.cutoff="q95") + NoAxes() + NoLegend()
pt3<-FeaturePlot(T15_int,c("M"),cols=c("gray", "blue"), pt.size=0.1, min.cutoff="q5", max.cutoff="q95") + NoAxes() + NoLegend()
pt4<-FeaturePlot(T15_int,c("MucSec"),cols=c("gray", "blue"), pt.size=0.1, min.cutoff="q5", max.cutoff="q95") + NoAxes() + NoLegend()
#dev.new(height = 4, width = 5)
png("FeaturePlots_ciliatedSec_markers.png",width=5,height=4,units="in",res=600)
plot_grid(pt1,pt2,pt3,pt4)
dev.off()


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





############## Box plots of all clusters
colorVec<-c("#1619BA","mediumpurple3","#3281FF","slategrey","aquamarine","turquoise4","tomato","darkgoldenrod",
	"darkorange2","yellowgreen","yellow3","wheat","saddlebrown","greenyellow","darkolivegreen","green3")

######O markers
#pdf("Boxplots_allCells_Omarkers.pdf")
dev.new(height=4,width=5)
par(bty="l")
boxplot(mean_exp_df$O~T15_int@meta.data$clusters_16a,las=1,col=colorVec,pch=16,cex=0.5,ylab="Mean normalized expression",
	names=rep("",16),medcol="black",outcex=0.3,medlwd=1,lwd=0.8,main="Early ciliating markers")
text(x=seq(1:16),y=par()$usr[3]-0.12*(par()$usr[4]-par()$usr[3]),labels=sort(unique(T15_int@meta.data$clusters_16a)),
	srt=45,adj=1,xpd=T,cex=0.8)
	
	
######N markers
#pdf("Boxplots_allCells_Nmarkers.pdf")
dev.new(height=4,width=5)
par(bty="l")
boxplot(mean_exp_df$N~T15_int@meta.data$clusters_16a,las=1,col=colorVec,pch=16,cex=0.5,ylab="Mean normalized expression",
	names=rep("",16),medcol="black",outcex=0.3,medlwd=1,lwd=0.8,main="Later ciliating markers")
text(x=seq(1:16),y=par()$usr[3]-0.12*(par()$usr[4]-par()$usr[3]),labels=sort(unique(T15_int@meta.data$clusters_16a)),
	srt=45,adj=1,xpd=T,cex=0.8)


######M markers
#pdf("Boxplots_allCells_Mmarkers.pdf")
dev.new(height=4,width=5)
par(bty="l")
boxplot(mean_exp_df$M~T15_int@meta.data$clusters_16a,las=1,col=colorVec,pch=16,cex=0.5,ylab="Mean normalized expression",
	names=rep("",16),medcol="black",outcex=0.3,medlwd=1,lwd=0.8,main="Mature ciliated markers")
text(x=seq(1:16),y=par()$usr[3]-0.12*(par()$usr[4]-par()$usr[3]),labels=sort(unique(T15_int@meta.data$clusters_16a)),
	srt=45,adj=1,xpd=T,cex=0.8)


######Secretory markers
#pdf("Boxplots_allCells_Secmarkers.pdf")
dev.new(height=4,width=5)
par(bty="l")
boxplot(mean_exp_df$MucSec~T15_int@meta.data$clusters_16a,las=1,col=colorVec,pch=16,cex=0.5,ylab="Mean normalized expression",
	names=rep("",16),medcol="black",outcex=0.3,medlwd=1,lwd=0.8,main="Secretory cell markers")
text(x=seq(1:16),y=par()$usr[3]-0.12*(par()$usr[4]-par()$usr[3]),labels=sort(unique(T15_int@meta.data$clusters_16a)),
	srt=45,adj=1,xpd=T,cex=0.8)






























############## Box plots of just the ciliated clusters (combine c8b and c8c)

#First need to make new metadata that combines all non-ciliated populations
treat.col<-data.frame("CiliatedGroups"=T15_int@meta.data$clusters_16a,row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col$CiliatedGroups[which(T15_int@meta.data$clusters_16a %in% c("c1","c1c","c2","c3a","c3b","c3c","c4","c7","c6a","c6b","c6c"))] <-"Non-ciliated"
treat.col$CiliatedGroups[which(T15_int@meta.data$clusters_16a %in% c("c5","c5M"))] <-"Mucus secretory"
treat.col$CiliatedGroups[which(T15_int@meta.data$clusters_16a %in% c("c8a"))] <-"Hybrid"
treat.col$CiliatedGroups[which(T15_int@meta.data$clusters_16a %in% c("c8b","c8c"))] <-"Mature"
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="CiliatedGroups")
T15_int@meta.data$CiliatedGroups<-factor(T15_int@meta.data$CiliatedGroups,levels=c("Hybrid","Mature","Non-ciliated","Mucus secretory"))

colorVec<-c("greenyellow","green3","grey","orange")

#pdf("Boxplots_ciliatedCells2_SFCilmarkers.pdf")
dev.new(height=6,width=3.7)
par(mfrow=c(2,2))
######O markers
#pdf("Boxplots_ciliatedCells2_Omarkers.pdf")
#dev.new(height=4,width=2.2)
par(bty="l")
boxplot(mean_exp_df$O~T15_int@meta.data$CiliatedGroups,las=1,col=colorVec,pch=16,cex=0.5,ylab="Mean normalized expression",
	names=rep("",4),medcol="black",outcex=0.3,medlwd=1.5,lwd=0.8,main="Early ciliating markers")
text(x=seq(1:4),y=par()$usr[3]-0.12*(par()$usr[4]-par()$usr[3]),labels=sort(unique(T15_int@meta.data$CiliatedGroups)),
	srt=45,adj=1,xpd=T,cex=0.8)
	
	
######N markers
#pdf("Boxplots_ciliatedCells2_Nmarkers.pdf")
#dev.new(height=4,width=2.2)
par(bty="l")
boxplot(mean_exp_df$N~T15_int@meta.data$CiliatedGroups,las=1,col=colorVec,pch=16,cex=0.5,ylab="Mean normalized expression",
	names=rep("",4),medcol="black",outcex=0.3,medlwd=1.5,lwd=0.8,main="Later ciliating markers")
text(x=seq(1:4),y=par()$usr[3]-0.12*(par()$usr[4]-par()$usr[3]),labels=sort(unique(T15_int@meta.data$CiliatedGroups)),
	srt=45,adj=1,xpd=T,cex=0.8)


######M markers
#pdf("Boxplots_ciliatedCells2_Mmarkers.pdf")
#dev.new(height=4,width=2.2)
par(bty="l")
boxplot(mean_exp_df$M~T15_int@meta.data$CiliatedGroups,las=1,col=colorVec,pch=16,cex=0.5,ylab="Mean normalized expression",
	names=rep("",4),medcol="black",outcex=0.3,medlwd=1.5,lwd=0.8,main="Mature ciliated markers")
text(x=seq(1:4),y=par()$usr[3]-0.12*(par()$usr[4]-par()$usr[3]),labels=sort(unique(T15_int@meta.data$CiliatedGroups)),
	srt=45,adj=1,xpd=T,cex=0.8)


######Secretory markers
#pdf("Boxplots_ciliatedCells2_Secmarkers.pdf")
#dev.new(height=4,width=2.2)
par(bty="l")
boxplot(mean_exp_df$MucSec~T15_int@meta.data$CiliatedGroups,las=1,col=colorVec,pch=16,cex=0.5,ylab="Mean normalized expression",
	names=rep("",4),medcol="black",outcex=0.3,medlwd=1.5,lwd=0.8,main="Secretory cell markers")
text(x=seq(1:4),y=par()$usr[3]-0.12*(par()$usr[4]-par()$usr[3]),labels=sort(unique(T15_int@meta.data$CiliatedGroups)),
	srt=45,adj=1,xpd=T,cex=0.8)


































########################### Finally, look at smoking DEGs in these ciliated populations	

#First, set up smoking status column for each of the 16a clusters
treat.col<-data.frame("clusters16a_smoke"=paste(T15_int@meta.data$clusters_16a,T15_int@meta.data$smoke,sep="_"),row.names=rownames(T15_int@meta.data))
treat.col[,1]<-as.character(treat.col[,1])
#Get rid of T89, who is too young
treat.col[which(T15_int@meta.data$donor == "T89"),]<-sapply(strsplit(treat.col[which(T15_int@meta.data$donor == "T89"),],"_"),function(x)paste(x[1],"exclude",sep="_"))
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="clusters16a_smoke")

#Now get up smoking DEGs in non-hybrid cells 
Idents(T15_int)<-T15_int@meta.data$clusters16a_smoke
c8_smoke_up<-FindMarkers(T15_int,ident.1 = c("c8b_heavy","c8c_heavy"),ident.2 = c("c8b_never","c8c_never"),min.pct = 0.1, 
	only.pos=T, logfc.threshold = 0.25, assay = "SCT", test.use = "wilcox")
	
#Now get up smoking DEGs in c8b non-hybrid cells 
Idents(T15_int)<-T15_int@meta.data$clusters16a_smoke
c8b_smoke_up<-FindMarkers(T15_int,ident.1 = "c8b_heavy",ident.2 = "c8b_never",min.pct = 0.1, 
	only.pos=T, logfc.threshold = 0.25, assay = "SCT", test.use = "wilcox")

#Now get up smoking DEGs in c8c non-hybrid cells 
Idents(T15_int)<-T15_int@meta.data$clusters16a_smoke
c8c_smoke_up<-FindMarkers(T15_int,ident.1 = "c8c_heavy",ident.2 = "c8c_never",min.pct = 0.1, 
	only.pos=T, logfc.threshold = 0.25, assay = "SCT", test.use = "wilcox")

#Now get down smoking DEGs in hybrid cells 
Idents(T15_int)<-T15_int@meta.data$clusters16a_smoke
c8a_hybrid_smoke_down<-FindMarkers(T15_int,ident.1 = "c8a_never",ident.2 = "c8a_heavy",min.pct = 0.1, 
	only.pos=T, logfc.threshold = 0.25, assay = "SCT", test.use = "wilcox")







####### Now bring in unique smoking DEGs stratified into c8a and c8bc
c8bc_up_unique<-read.table("DEG_Tables/Smoking_up_DEGs_coreUniqueDesignations.txt",header=T,stringsAsFactors=F)
c8a_down_unique<-read.table("DEG_Tables/Smoking_down_DEGs_coreUniqueDesignations.txt",header=T,stringsAsFactors=F)
c8bc_up_unique<-c8bc_up_unique[which(c8bc_up_unique$comparison=="c8bc" & c8bc_up_unique$main_cil == "Unique"),]
c8a_down_unique<-c8a_down_unique[which(c8a_down_unique$comparison=="c8a" & c8a_down_unique$main_cil == "Unique"),]




#Now get mean expression for each of these sets of unique genes
#Pull out markers to plot
markerList<-list("c8bc_up_unique"=c8bc_up_unique$gene,"c8a_down_unique"=c8a_down_unique$gene)

#Now get mean expression of these markers
mean_exp_df<-CalculateMeanExpression(dataset=as.matrix(GetAssayData(T15_int,assay="SCT")),genesList=markerList,method="geometric")

#Add these markers to metadata
T15_int<-AddMetaData(T15_int,mean_exp_df)








#Now get metadata to isolate hybrid, mature, and non-ciliated cells
#First need to make new metadata that combines all non-ciliated populations
treat.col<-data.frame("CiliatedGroups_mergeMature"=T15_int@meta.data$clusters_16a,row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col$CiliatedGroups_mergeMature[which(T15_int@meta.data$clusters_16a %in% c("c1","c1c","c2","c3a","c3b","c3c","c4","c5","c5M","c7","c6a","c6b","c6c"))] <-"Non-ciliated"
treat.col$CiliatedGroups_mergeMature[which(T15_int@meta.data$clusters_16a %in% c("c8a"))] <-"Hybrid"
treat.col$CiliatedGroups_mergeMature[which(T15_int@meta.data$clusters_16a %in% c("c8b","c8c"))] <-"Mature"
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="CiliatedGroups_mergeMature")
T15_int@meta.data$CiliatedGroups_mergeMature<-factor(T15_int@meta.data$CiliatedGroups_mergeMature,levels=c("Hybrid","Mature","Non-ciliated"))

#Then make combined CiliatedGroups + smoking metadata column
treat.col<-data.frame("CiliatedGroups_mergeMature_smoke"=paste(T15_int@meta.data$CiliatedGroups_mergeMature,T15_int@meta.data$smoke,sep="_"),row.names=rownames(T15_int@meta.data))
treat.col[,1]<-as.character(treat.col[,1])
#Get rid of T89, who is too young
treat.col[which(T15_int@meta.data$donor == "T89"),]<-sapply(strsplit(treat.col[which(T15_int@meta.data$donor == "T89"),],"_"),function(x)paste(x[1],"exclude",sep="_"))
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="CiliatedGroups_mergeMature_smoke")

#Then isolate dataset that exclude the light and T89 smoker
T15_int_temp<-subset(T15_int,cells=rownames(T15_int@meta.data)[-grep("light|exclude",T15_int@meta.data$CiliatedGroups_mergeMature_smoke)])
T15_int_temp@meta.data$CiliatedGroups_mergeMature_smoke<-factor(T15_int_temp@meta.data$CiliatedGroups_mergeMature_smoke,
	levels=c("Hybrid_never","Hybrid_heavy","Mature_never","Mature_heavy","Non-ciliated_never","Non-ciliated_heavy"))

colorVec<-rep(c("black","red"),3)


dev.new(height=3,width=5.5)
par(mfrow=c(1,2),bty="l")
#pdf("Boxplots_ciliatedCells_smokingEffects_mergeMature.pdf")
#dev.new(height=4,width=2.2)
boxplot(T15_int_temp@meta.data$c8bc_up_unique~T15_int_temp@meta.data$CiliatedGroups_mergeMature_smoke,las=1,col=colorVec,
	pch=16,cex=0.5,ylab="Mean normalized expression",names=rep("",6),medcol="white",outcex=0.3,medlwd=1.5,lwd=0.8,
	main="Unique smoke \nup mature",cex.main=0.8,at=c(1,2,4,5,7,8))
text(x=c(1,2,4,5,7,8),y=par()$usr[3]-0.12*(par()$usr[4]-par()$usr[3]),labels=sort(unique(T15_int_temp@meta.data$CiliatedGroups_mergeMature_smoke)),
	srt=45,adj=1,xpd=T,cex=0.8)
	
boxplot(T15_int_temp@meta.data$c8a_down_unique~T15_int_temp@meta.data$CiliatedGroups_mergeMature_smoke,las=1,col=colorVec,
	pch=16,cex=0.5,ylab="Mean normalized expression",names=rep("",6),medcol="white",outcex=0.3,medlwd=1.5,lwd=0.8,
	main="Unique smoke \ndown hybrid",cex.main=0.8,at=c(1,2,4,5,7,8))
text(x=c(1,2,4,5,7,8),y=par()$usr[3]-0.12*(par()$usr[4]-par()$usr[3]),labels=sort(unique(T15_int_temp@meta.data$CiliatedGroups_mergeMature_smoke)),
	srt=45,adj=1,xpd=T,cex=0.8)




#Do t tests of differences
#c8 smoke up (5.829435e-01 4.578781e-36 1.266509e-04)
pvals<-c()
count<-1
for(i in 1:(length(levels(T15_int_temp@meta.data$CiliatedGroups_mergeMature_smoke)) / 2)){
	pvals<-append(pvals,wilcox.test(T15_int_temp@meta.data$c8bc_up_unique[which(T15_int_temp@meta.data$CiliatedGroups_mergeMature_smoke == 
		levels(T15_int_temp@meta.data$CiliatedGroups_mergeMature_smoke)[count + 1])],
		T15_int_temp@meta.data$c8bc_up_unique[which(T15_int_temp@meta.data$CiliatedGroups_mergeMature_smoke == 
		levels(T15_int_temp@meta.data$CiliatedGroups_mergeMature_smoke)[count])],
		alternative="two.sided")$p.value)
	count<-count + 2
}


#c8 hybrid smoke down (5.505093e-09 8.863205e-01 1.021443e-02)
pvals<-c()
count<-1
for(i in 1:(length(levels(T15_int_temp@meta.data$CiliatedGroups_mergeMature_smoke)) / 2)){
	pvals<-append(pvals,wilcox.test(T15_int_temp@meta.data$c8a_down_unique[which(T15_int_temp@meta.data$CiliatedGroups_mergeMature_smoke == 
		levels(T15_int_temp@meta.data$CiliatedGroups_mergeMature_smoke)[count + 1])],
		T15_int_temp@meta.data$c8a_down_unique[which(T15_int_temp@meta.data$CiliatedGroups_mergeMature_smoke == 
		levels(T15_int_temp@meta.data$CiliatedGroups_mergeMature_smoke)[count])],
		alternative="two.sided")$p.value)
	count<-count + 2
}






