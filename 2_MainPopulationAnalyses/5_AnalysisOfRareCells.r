#########################################
#############Seurat analysis#############
#########################################

library(Seurat, lib.loc="/Users/jacksonna/R-dev")
library(plotrix)
library(scales)
library(ggplot2)
library(cowplot)
library(openxlsx)
library(vioplot)
library(beeswarm)
source("CommonFunctions.r")



#Download the processed Seurat R object, "Processed_invitro_seurat.Rdata" from GEO,
#obtainable from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134174
#Then load the object it contains ("T15_int") into the session.
#T15_int is a Seurat dataset that contains the processed, integrated dataset
load("Processed_invitro_seurat.Rdata")




###################### First, make beeswarms/violin plots for rare cell markers

#First need to make new metadata that allocates populations
treat.col<-data.frame("RareGroups"=T15_int@meta.data$clusters_16a,row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col$RareGroups[which(T15_int@meta.data$clusters_10 %in% c("c1","c1c","c2","c3","c4","c5","c7","c8"))] <-"Other"
treat.col$RareGroups[which(T15_int@meta.data$clusters_16a %in% c("c6a"))] <-"Ionocytes"
treat.col$RareGroups[which(T15_int@meta.data$clusters_16a %in% c("c6b"))] <-"Tuft-like"
treat.col$RareGroups[which(T15_int@meta.data$clusters_16a %in% c("c6c"))] <-"PNEC"
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="RareGroups")
T15_int@meta.data$RareGroups<-factor(T15_int@meta.data$RareGroups,levels=c("Ionocytes","Tuft-like","PNEC","Other"))

#Violin plot dataframe
T15_sct<-as.matrix(GetAssayData(T15_int,assay="SCT"))

#Make violin plots

###CLCNKB
#pdf("ViolinPlots_CLCNKB_rareCells.pdf")
dev.new(height=3.5,width=4)
vioplot(T15_sct["POU2F3",][which(T15_int@meta.data$RareGroups == "Ionocytes")],
	T15_sct["CLCNKB",][which(T15_int@meta.data$RareGroups == "Tuft-like")],
	T15_sct["CLCNKB",][which(T15_int@meta.data$RareGroups == "PNEC")],
	T15_sct["CLCNKB",][which(T15_int@meta.data$RareGroups == "Other")],
	drawRect=F,las=1,at=1:4,names=rep("",4),col=c("grey"))
jitterCoordinates<-jitter(as.numeric(T15_int@meta.data$RareGroups[which(T15_int@meta.data$RareGroups != "Other")]))
points(T15_sct["CLCNKB",which(T15_int@meta.data$RareGroups != "Other")]~jitterCoordinates,col="black",bg="black",pch=21,cex=0.8,lwd=0.5)
#beeswarm(T15_sct["CLCNKB",which(T15_int@meta.data$RareGroups != "Other")]~T15_int@meta.data$RareGroups[which(T15_int@meta.data$RareGroups != "Other")],
#	method="swarm", corral="random",add=T,pch=16,cex=0.8,col="black")

###POU2F3
#pdf("ViolinPlots_POU2F3_rareCells.pdf")
dev.new(height=3.5,width=3)
vioplot(T15_sct["POU2F3",][which(T15_int@meta.data$RareGroups == "Tuft-like")],
	T15_sct["POU2F3",][which(T15_int@meta.data$RareGroups == "Other")],
	drawRect=F,las=1,col=c("grey"))
jitterCoordinates<-jitter((as.numeric(T15_int@meta.data$RareGroups[which(T15_int@meta.data$RareGroups == "Tuft-like")]) - 1),amount=0.2)
points(T15_sct["POU2F3",which(T15_int@meta.data$RareGroups == "Tuft-like")]~jitterCoordinates,col="black",bg="black",pch=21,cex=0.8,lwd=0.5)
#beeswarm(T15_sct["POU2F3",which(T15_int@meta.data$RareGroups == "Tuft-like")]~T15_int@meta.data$RareGroups[which(T15_int@meta.data$RareGroups == "Tuft-like")],
#	method="swarm", corral="random",add=T,pch=16,cex=0.8,col="black")

###CALCA
#pdf("ViolinPlots_CALCA_rareCells.pdf")
dev.new(height=3.5,width=4)
vioplot(T15_sct["CALCA",][which(T15_int@meta.data$RareGroups == "Ionocytes")],
	T15_sct["CALCA",][which(T15_int@meta.data$RareGroups == "Tuft-like")],
	T15_sct["CALCA",][which(T15_int@meta.data$RareGroups == "PNEC")],
	T15_sct["CALCA",][which(T15_int@meta.data$RareGroups == "Other")],
	drawRect=F,las=1,at=1:4,names=rep("",4),col=c("grey"))
jitterCoordinates<-jitter(as.numeric(T15_int@meta.data$RareGroups[which(T15_int@meta.data$RareGroups != "Other")]))
points(T15_sct["CALCA",which(T15_int@meta.data$RareGroups != "Other")]~jitterCoordinates,col="black",bg="black",pch=21,cex=0.8,lwd=0.5)
#beeswarm(T15_sct["CALCA",which(T15_int@meta.data$RareGroups != "Other")]~T15_int@meta.data$RareGroups[which(T15_int@meta.data$RareGroups != "Other")],
#	method="swarm", corral="random",add=T,pch=16,cex=0.8,col="black")
	
	
	







	
	
	






############################### Create CFTR feature plot
#Based on Dyjack code:
library(DescTools)

theme_ele_tsne_2 <- theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_blank(),
                    axis.text = element_blank(), axis.ticks = element_blank(),
                    axis.title= element_blank(), plot.title = element_blank())

x <- "CFTR"
ggdat <- data.frame(x=Embeddings(T15_int, reduction = "umap")[,1],y=Embeddings(T15_int, reduction = "umap")[,2],
                expr=as.matrix(GetAssayData(T15_int,assay="RNA",slot="data"))[x,],exps=as.matrix(GetAssayData(T15_int,assay="SCT"))[x,])
p <- ggplot(ggdat,aes(x=x,y=y)) + theme_bw() + theme_ele_tsne_2
ggdat$exps <- Winsorize(ggdat$exps,probs = c(0.005, 0.995))
p <- p + geom_point(data=subset(ggdat,expr==0),size = 1,color="gray88") #to distinguish border and fill do shape = 21, colore = "black", fill = "white"
p <- p + geom_point(data=subset(ggdat,expr>0),aes(color=exps),size = 1)
p <- p + scale_color_gradient(x,low="#ffe6e6ff",high="#ff0000ff") + xlab("") + ylab("")
#pdf("T15_featureplot_CFTR.pdf",height=6,width=7)
dev.new(height=6,width=7)
print(p)
#dev.off()








############# Now calculate the proportion of cells in T15 populations expressing CFTR that are ionocytes (and compare with proportion of ionocytes in the group)
T15_counts<-as.matrix(GetAssayData(T15_int,assay="RNA",slot="counts"))
#Total number of reads
sum(T15_counts["CFTR",])

#Total number of reads expressed in c6a
sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_16a == "c6a")])
1104/9848 #11.2%

#Proportion of cells that are ionocytes
nrow(T15_int@meta.data[which(T15_int@meta.data$clusters_16a == "c6a"),]) / nrow(T15_int@meta.data) # 0.0028

#Number of CFTR reads per cell
#ionocytes
1104 / nrow(T15_int@meta.data[which(T15_int@meta.data$clusters_16a == "c6a"),]) # 10.93
#Non-ionocytes
8744 / nrow(T15_int@meta.data[-which(T15_int@meta.data$clusters_16a == "c6a"),]) # 0.24


############ Make table for Figure 5
cftr_table<-data.frame(matrix(nrow=11,ncol=4))
colnames(cftr_table)<-c("in vivo cell populations","Average CFTR UMI/cell","Cell equivalents needed per ionocyte","% of total CFTR UMI")
cftr_table[,1]<-c("Proliferating basal","Differentiating basal","KRT8Hight intermediate","Mucus secretory","Ciliated",
	"Ionocyte","Tuft-like","PNEC","SMG basal 1","SMG basal 2","SMG secretory")

#Get the average number of CFTR UMIs/cell
cftr_table[1,2]<-round(sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c1")]) / length(which(T15_int@meta.data$clusters_10 == "c1")),2)
cftr_table[2,2]<-round(sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c2")]) / length(which(T15_int@meta.data$clusters_10 == "c2")),2)
cftr_table[3,2]<-round(sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c4")]) / length(which(T15_int@meta.data$clusters_10 == "c4")),2)
cftr_table[4,2]<-round(sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c5")]) / length(which(T15_int@meta.data$clusters_10 == "c5")),2)
cftr_table[5,2]<-round(sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c8")]) / length(which(T15_int@meta.data$clusters_10 == "c8")),2)
cftr_table[6,2]<-round(sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_16a == "c6a")]) / length(which(T15_int@meta.data$clusters_16a == "c6a")),2)
cftr_table[7,2]<-round(sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_16a == "c6b")]) / length(which(T15_int@meta.data$clusters_16a == "c6b")),2)
cftr_table[8,2]<-round(sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_16a == "c6c")]) / length(which(T15_int@meta.data$clusters_16a == "c6c")),2)
cftr_table[9,2]<-round(sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c1c")]) / length(which(T15_int@meta.data$clusters_10 == "c1c")),2)
cftr_table[10,2]<-round(sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c3")]) / length(which(T15_int@meta.data$clusters_10 == "c3")),2)
cftr_table[11,2]<-round(sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c7")]) / length(which(T15_int@meta.data$clusters_10 == "c7")),2)

#Get number of cells for each cluster to make as much CFTR as a single ionocyte
cftr_table[1,3]<-round(10.93 / round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c1")]) / length(which(T15_int@meta.data$clusters_10 == "c1"))),2),1)
cftr_table[2,3]<-round(10.93 / round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c2")]) / length(which(T15_int@meta.data$clusters_10 == "c2"))),2),1)
cftr_table[3,3]<-round(10.93 / round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c4")]) / length(which(T15_int@meta.data$clusters_10 == "c4"))),2),1)
cftr_table[4,3]<-round(10.93 / round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c5")]) / length(which(T15_int@meta.data$clusters_10 == "c5"))),2),1)
cftr_table[5,3]<-round(10.93 / round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c8")]) / length(which(T15_int@meta.data$clusters_10 == "c8"))),2),1)
cftr_table[6,3]<-round(10.93 / round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_16a == "c6a")]) / length(which(T15_int@meta.data$clusters_16a == "c6a"))),2),1)
cftr_table[7,3]<-round(10.93 / round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_16a == "c6b")]) / length(which(T15_int@meta.data$clusters_16a == "c6b"))),2),1)
cftr_table[8,3]<-round(10.93 / round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_16a == "c6c")]) / length(which(T15_int@meta.data$clusters_16a == "c6c"))),2),1)
cftr_table[9,3]<-round(10.93 / round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c1c")]) / length(which(T15_int@meta.data$clusters_10 == "c1c"))),2),1)
cftr_table[10,3]<-round(10.93 / round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c3")]) / length(which(T15_int@meta.data$clusters_10 == "c3"))),2),1)
cftr_table[11,3]<-round(10.93 / round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c7")]) / length(which(T15_int@meta.data$clusters_10 == "c7"))),2),1)

#Get % of total CFTR UMIs per cluster
cftr_table[1,4]<-round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c1")])/9848) * 100,1)
cftr_table[2,4]<-round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c2")])/9848) * 100,1)
cftr_table[3,4]<-round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c4")])/9848) * 100,1)
cftr_table[4,4]<-round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c5")])/9848) * 100,1)
cftr_table[5,4]<-round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c8")])/9848) * 100,1)
cftr_table[6,4]<-round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_16a == "c6a")])/9848) * 100,1)
cftr_table[7,4]<-round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_16a == "c6b")])/9848) * 100,1)
cftr_table[8,4]<-round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_16a == "c6c")])/9848) * 100,1)
cftr_table[9,4]<-round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c1c")])/9848) * 100,1)
cftr_table[10,4]<-round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c3")])/9848) * 100,1)
cftr_table[11,4]<-round((sum(T15_counts["CFTR",which(T15_int@meta.data$clusters_10 == "c7")])/9848) * 100,1)

#Now export table
write.table(cftr_table,quote=F,row.names=F,sep="\t",file="CFTR_distribution_table.txt")



























########################## Make heat map of uniquely upregulated genes in each of the three rare subgroups 
########################## There are options for filtering for the integrated dataset if desired
#Define integrated dataset
T15_integrated<-as.matrix(GetAssayData(T15_int,assay="integrated"))

#First, read in the rare3 DEGs
rare3_DEGs<-read.table("Subclustering_results/rare3/DEG_Tables/Rare3_DEGs.txt",header=T,stringsAsFactors=F)
#Isolate rare and nonrare comparison datasets
rare3_DEGs_nonrare<-rare3_DEGs[-grep("vs_rare",rare3_DEGs$comparison),]
rare3_DEGs_rare<-rare3_DEGs[-grep("vs_nonrare",rare3_DEGs$comparison),]
#Toss FDR > 0.05
rare3_DEGs_nonrare<-rare3_DEGs_nonrare[which(rare3_DEGs_nonrare$p_adj_FDR < 0.05),]
rare3_DEGs_rare<-rare3_DEGs_rare[which(rare3_DEGs_rare$p_adj_FDR < 0.05),]


##Now get lists of genes uniquely up in each rare cluster compared to non-rare cells
#ion_genes<-rare3_DEGs_nonrare$gene[which(rare3_DEGs_nonrare$comparison == "Ion_vs_nonrare" & 
#	!(rare3_DEGs_nonrare$gene %in% rare3_DEGs_nonrare$gene[which(duplicated(rare3_DEGs_nonrare$gene))]))]
#tuft_genes<-rare3_DEGs_nonrare$gene[which(rare3_DEGs_nonrare$comparison == "Tuft_vs_nonrare" & 
#	!(rare3_DEGs_nonrare$gene %in% rare3_DEGs_nonrare$gene[which(duplicated(rare3_DEGs_nonrare$gene))]))]
#PNEC_genes<-rare3_DEGs_nonrare$gene[which(rare3_DEGs_nonrare$comparison == "PNEC_vs_nonrare" & 
#	!(rare3_DEGs_nonrare$gene %in% rare3_DEGs_nonrare$gene[which(duplicated(rare3_DEGs_nonrare$gene))]))]
#
##Even more stringently, get global DEGs that are also a local DEGs
#ion_genes<-rare3_DEGs_nonrare$gene[which(rare3_DEGs_nonrare$comparison == "Ion_vs_nonrare" & rare3_DEGs_nonrare$gene %in% 
#	rare3_DEGs_rare$gene[which(rare3_DEGs_rare$comparison == "Ion_vs_rare")])]
#tuft_genes<-rare3_DEGs_nonrare$gene[which(rare3_DEGs_nonrare$comparison == "Tuft_vs_nonrare" & rare3_DEGs_nonrare$gene %in% 
#	rare3_DEGs_rare$gene[which(rare3_DEGs_rare$comparison == "Tuft_vs_rare")])]
#PNEC_genes<-rare3_DEGs_nonrare$gene[which(rare3_DEGs_nonrare$comparison == "PNEC_vs_nonrare" & rare3_DEGs_nonrare$gene %in% 
#	rare3_DEGs_rare$gene[which(rare3_DEGs_rare$comparison == "PNEC_vs_rare")])]

#Most stringently, get global DEGs that are both unique to a single rare cell type AND are a local DEG
ion_genes<-rare3_DEGs_nonrare$gene[which(rare3_DEGs_nonrare$comparison == "Ion_vs_nonrare" & rare3_DEGs_nonrare$gene %in% 
	rare3_DEGs_rare$gene[which(rare3_DEGs_rare$comparison == "Ion_vs_rare")] & 
	!(rare3_DEGs_nonrare$gene %in% rare3_DEGs_nonrare$gene[which(duplicated(rare3_DEGs_nonrare$gene))]))]
tuft_genes<-rare3_DEGs_nonrare$gene[which(rare3_DEGs_nonrare$comparison == "Tuft_vs_nonrare" & rare3_DEGs_nonrare$gene %in% 
	rare3_DEGs_rare$gene[which(rare3_DEGs_rare$comparison == "Tuft_vs_rare")] & 
	!(rare3_DEGs_nonrare$gene %in% rare3_DEGs_nonrare$gene[which(duplicated(rare3_DEGs_nonrare$gene))]))]
PNEC_genes<-rare3_DEGs_nonrare$gene[which(rare3_DEGs_nonrare$comparison == "PNEC_vs_nonrare" & rare3_DEGs_nonrare$gene %in% 
	rare3_DEGs_rare$gene[which(rare3_DEGs_rare$comparison == "PNEC_vs_rare")] & 
	!(rare3_DEGs_nonrare$gene %in% rare3_DEGs_nonrare$gene[which(duplicated(rare3_DEGs_nonrare$gene))]))]

##Filter genes for those in the integrated dataset
#ion_genes<-ion_genes[which(ion_genes %in% rownames(T15_integrated))]
#tuft_genes<-tuft_genes[which(tuft_genes %in% rownames(T15_integrated))]
#PNEC_genes<-PNEC_genes[which(PNEC_genes %in% rownames(T15_integrated))]

##Now make heat map
#genesIcareAbout<-c("ASCL2","POU2F3","HCK","SOX4","LRMP","GRP","CALCA","CHGA","CHGB GP2","ASCL1","SYT1","SPOCK1",
#	"CFTR","ATP6V1B1","STAP1","RARRES2","DEFB1","SCNN1B","CEL","CLCNKB")
#genesIcareAbout %in% c(ion_genes,tuft_genes,PNEC_genes)

#Randomize order of genes
ion_genes<-sample(ion_genes,length(ion_genes),replace=FALSE)
tuft_genes<-sample(tuft_genes,length(tuft_genes),replace=FALSE)
PNEC_genes<-sample(PNEC_genes,length(PNEC_genes),replace=FALSE)

cols_clusters<-c("yellowgreen","yellow3","wheat")
cols_donor<-c("midnightblue","black","red","pink","tan","greenyellow","orangered4","yellow","lightblue",
	"grey","darkgreen","blue","purple4","palegreen","slateblue")

#Subset data
T15_rare<-subset(T15_int,cells=rownames(T15_int@meta.data)[grep("c6",T15_int@meta.data$clusters_16a)])

#First, need to isolate those genes we want and cells we want (in a randomized order to minimize donor differences)
###############
cellOrdering<-c(
	sample(which(T15_rare@meta.data$clusters_16a == "c6a"),
		length(which(T15_rare@meta.data$clusters_16a == "c6a")),replace=FALSE),
	sample(which(T15_rare@meta.data$clusters_16a == "c6b"),
		length(which(T15_rare@meta.data$clusters_16a == "c6b")),replace=FALSE),
	sample(which(T15_rare@meta.data$clusters_16a == "c6c"),
		length(which(T15_rare@meta.data$clusters_16a == "c6c")),replace=FALSE))		

#Get gene ordering based on DEGs
geneOrdering<-rev(c(ion_genes,tuft_genes,PNEC_genes))

#Get module ordering (cluster for targeted DEGs)
moduleOrdering<-c(rep("yellowgreen",length(ion_genes)),rep("yellow3",length(tuft_genes)),rep("wheat",length(PNEC_genes)))
moduleOrdering<-rev(moduleOrdering)

#Now make a heatmap using heatmap3 of expression across genes and cells
mat<-as.matrix(GetAssayData(T15_rare,assay="SCT",slot="data"))[geneOrdering,cellOrdering]
#mat<-as.matrix(GetAssayData(T15_rare,assay="integrated",slot="data"))[geneOrdering,cellOrdering]

#Get side colors for population
#Note that the colors need to be in numerical order, not in the order of cellOrdering
clusters<-color.scale(as.numeric(as.factor(T15_rare@meta.data[cellOrdering,]$clusters_16a)),extremes=cols_clusters,color.spec="rgb")
#donors<-color.scale(as.numeric(as.factor(T15_rare@meta.data[cellOrdering,]$donor)),extremes=cols_donor,color.spec="rgb")
ColSideColors<-data.frame(clusters=clusters)

#Get row colors for modules
modules<-color.scale(as.numeric(as.factor(moduleOrdering)),extremes=sort(unique(moduleOrdering)), color.spec="rgb")
RowSideColors<-data.frame(modules=modules)

mycols = colorRampPalette(c("blue","white","red"))(1000)
breakscale = c(-8,seq(-1.8,1.8, length.out=length(mycols)-1),8)

pdf('T15_rare_markersSCT_heatmap.pdf',height=24,width=18)
heatmap3(mat, scale="row",Rowv=NA,Colv=NA, 
	col=mycols,ColSideColors=as.matrix(ColSideColors),RowSideColors=as.matrix(RowSideColors),
	cexRow=0.2,cexCol=1,breaks=breakscale,labCol="",useRaster=T)
dev.off()






#Write heat map DEGs to files
writeDEGsToExcel(rare3_DEGs_nonrare[which(rare3_DEGs_nonrare$gene %in% c(ion_genes,tuft_genes,PNEC_genes)),],"Rare_heatmap_DEGs_SCT.xlsx")
writeDEGsToExcel(rare3_DEGs_nonrare[which(rare3_DEGs_nonrare$gene %in% c(ion_genes,tuft_genes,PNEC_genes) &
	rare3_DEGs_nonrare$gene %in% rownames(T15_integrated)),],"Rare_heatmap_DEGs_integrated.xlsx")

#Look at transcription regulators among these
IPA<-read.table("Gene_lists/database_genes_IPA.txt",sep="\t",header=T,stringsAsFactors=F)
blah<-rare3_DEGs_nonrare[which(rare3_DEGs_nonrare$gene %in% c(ion_genes,tuft_genes,PNEC_genes) &
	rare3_DEGs_nonrare$gene %in% rownames(T15_integrated)),]
blah<-blah[which(blah$gene %in% IPA$Gene[which(IPA$Type_IPA == "transcription regulator")]),]

#Do enrichments on these gene sets
#Write all to single txt file
write.table(rare3_DEGs_nonrare[which(rare3_DEGs_nonrare$gene %in% c(ion_genes,tuft_genes,PNEC_genes)),],sep="\t",quote=FALSE,col.names=TRUE,
	row.names=FALSE,file="Enrichr/Rare_heatmap_DEGs_SCT_forEnrichr.txt")
write.table(rare3_DEGs_nonrare[which(rare3_DEGs_nonrare$gene %in% c(ion_genes,tuft_genes,PNEC_genes) &
	rare3_DEGs_nonrare$gene %in% rownames(T15_integrated)),],sep="\t",quote=FALSE,col.names=TRUE,
	row.names=FALSE,file="Enrichr/Rare_heatmap_DEGs_integrated_forEnrichr.txt")

#Do enrichments
DEG_table<-read.table("Enrichr/Rare_heatmap_DEGs_SCT_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"Rare_heatmap_DEGs_SCT"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)

DEG_table<-read.table("Enrichr/Rare_heatmap_DEGs_integrated_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"Rare_heatmap_DEGs_integrated"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)






























################################ Make violin plot of FOXI1 expression
#Violin plot dataframe
T15_sct<-as.matrix(GetAssayData(T15_int,assay="SCT"))
#T15_rna<-as.matrix(GetAssayData(T15_int,assay="RNA")) #Use RNA to determine whether is gene is expressed or not
#First need to make new metadata that allocates clusters of interest
treat.col<-data.frame("Rare_groups"=T15_int@meta.data$clusters_16a,row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col$Rare_groups[which(T15_int@meta.data$clusters_16a %in% c("c6a"))] <-"Ionocytes"
treat.col$Rare_groups[which(T15_int@meta.data$clusters_16a %in% c("c6b"))] <-"Tuft-like"
treat.col$Rare_groups[which(T15_int@meta.data$clusters_16a %in% c("c6c"))] <-"PNEC"
treat.col$Rare_groups[which(T15_int@meta.data$clusters_16a %in% grep("c6",unique(T15_int@meta.data$clusters_16a),value=T,invert=T))] <-"Other"
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="Rare_groups")
T15_int@meta.data$Rare_groups<-factor(T15_int@meta.data$Rare_groups,levels=c("Ionocytes","Tuft-like","PNEC","Other"))

dev.new(height=3.5,width=4)
vioplot(T15_sct["FOXI1",][which(T15_int@meta.data$Rare_groups == "Ionocytes")],
	T15_sct["FOXI1",][which(T15_int@meta.data$Rare_groups == "Tuft-like")],
	T15_sct["FOXI1",][which(T15_int@meta.data$Rare_groups == "PNEC")],
	T15_sct["FOXI1",][which(T15_int@meta.data$Rare_groups == "Other")],
	drawRect=F,las=1,at=1:4,names=rep("",4),col=c("grey"))

#Get color vecs for CFTR, POU2F3, and FOXI1 expression profiles
cftrPos<-colnames(T15_sct[,which(T15_sct["CFTR",] > 0 & T15_sct["POU2F3",] == 0)])
pou2f3Pos<-colnames(T15_sct[,which(T15_sct["POU2F3",] > 0 & T15_sct["CFTR",] == 0)])
ctfrPou2f3Pos<-colnames(T15_sct[,which(T15_sct["POU2F3",] > 0 & T15_sct["CFTR",] > 0)])
foxi1Pos<-colnames(T15_sct[,which(T15_sct["FOXI1",] > 0 & T15_sct["CFTR",] == 0 & T15_sct["POU2F3",] == 0)])
#Exclude "Other" since too many to plot
ColorVec<-rep("grey66",ncol(T15_sct[,which(T15_int@meta.data$Rare_groups != "Other")]))
ColorVec[which(colnames(T15_sct[,which(T15_int@meta.data$Rare_groups != "Other")]) %in% cftrPos)]<-"red"
ColorVec[which(colnames(T15_sct[,which(T15_int@meta.data$Rare_groups != "Other")]) %in% pou2f3Pos)]<-"green"
ColorVec[which(colnames(T15_sct[,which(T15_int@meta.data$Rare_groups != "Other")]) %in% ctfrPou2f3Pos)]<-"black"
ColorVec[which(colnames(T15_sct[,which(T15_int@meta.data$Rare_groups != "Other")]) %in% foxi1Pos)]<-"white"

#Now add in beeswarm
#Get jitter coordinates
#pdf(Rare_FOXI1_violinPlots.pdf")
jitterCoordinates<-jitter(as.numeric(T15_int@meta.data$Rare_groups[which(T15_int@meta.data$Rare_groups != "Other")]))
points(T15_sct["FOXI1",which(T15_int@meta.data$Rare_groups != "Other")]~jitterCoordinates,col="black",bg=ColorVec,pch=21,cex=1,lwd=0.5)


















#Make pie chart for proportion of FOXI1- ionocytes with CFTR and without CFTR
#and the proportion of FOXI1- Tufts with Pou2F3 and without Pou2F3

#First, get all the FOXI1- cells for ionos and tufts
#I tried RNA, and the proportions are very similar
T15_sct<-as.matrix(GetAssayData(T15_int,assay="SCT"))
T15_iono_noFOXI1<-subset(T15_int,cells=rownames(T15_int@meta.data)[which(T15_int@meta.data$clusters_16a == "c6a" & T15_sct["FOXI1",] == 0)])
T15_tuft_noFOXI1<-subset(T15_int,cells=rownames(T15_int@meta.data)[which(T15_int@meta.data$clusters_16a == "c6b" & T15_sct["FOXI1",] == 0)])




#Get proportions for ionocytes
T15_sct_iono_noFOXI1<-as.matrix(GetAssayData(T15_iono_noFOXI1,assay="SCT"))
proportionsList<-list()
for(i in 1:1){
	currData<-T15_sct_iono_noFOXI1[c("CFTR","POU2F3"),]
	proportions<-c()
	proportions<-append(proportions,length(which(currData["CFTR",] != 0 & currData["POU2F3",] == 0)) / ncol(currData))
	proportions<-append(proportions,length(which(currData["POU2F3",] != 0 & currData["CFTR",] == 0)) / ncol(currData))
	proportions<-append(proportions,length(which(currData["CFTR",] != 0 & currData["POU2F3",] != 0)) / ncol(currData))
	proportions<-append(proportions,length(which(currData["CFTR",] == 0 & currData["POU2F3",] == 0)) / ncol(currData))
	names(proportions)<-c("CFTR only","POU2F3 only","Both","Neither")
	proportionsList[[length(proportionsList) + 1]]<-proportions
}

#Plot pie
#pdf("T15_ionocytes_CFTR_piechart.pdf")
for(i in 1:length(proportionsList)){
	pct <- round(proportionsList[[i]],3)*100
	lbls <- paste(pct,"%",sep="") # ad % to labels 
	colors<-c("red","green","black","grey66")
	pie(proportionsList[[i]],labels = names(proportionsList[[i]]), col=colors, main = names(proportionsList)[[i]])
}




#Get proportions for Tufts
T15_sct_tuft_noFOXI1<-as.matrix(GetAssayData(T15_tuft_noFOXI1,assay="SCT"))
proportionsList<-list()
for(i in 1:1){
	currData<-T15_sct_tuft_noFOXI1[c("CFTR","POU2F3"),]
	proportions<-c()
	proportions<-append(proportions,length(which(currData["CFTR",] != 0 & currData["POU2F3",] == 0)) / ncol(currData))
	proportions<-append(proportions,length(which(currData["POU2F3",] != 0 & currData["CFTR",] == 0)) / ncol(currData))
	proportions<-append(proportions,length(which(currData["CFTR",] != 0 & currData["POU2F3",] != 0)) / ncol(currData))
	proportions<-append(proportions,length(which(currData["CFTR",] == 0 & currData["POU2F3",] == 0)) / ncol(currData))
	names(proportions)<-c("CFTR only","POU2F3 only","Both","Neither")
	proportionsList[[length(proportionsList) + 1]]<-proportions
}

#Plot pie
#pdf("T15_tuft_POU2F3_piechart.pdf")
for(i in 1:length(proportionsList)){
	pct <- round(proportionsList[[i]],3)*100
	lbls <- paste(pct,"%",sep="") # ad % to labels 
	colors<-c("red","green","black","grey66")
	pie(proportionsList[[i]],labels = names(proportionsList[[i]]), col=colors, main = names(proportionsList)[[i]])
}


















#Repete the pie chart, this time for proportion of FOXI1+ ionocytes with CFTR and without CFTR
#and the proportion of FOXI1- Tufts with Pou2F3 and without Pou2F3

#First, get all the FOXI1+ cells for ionos and tufts
#I tried RNA, and the proportions are very similar
T15_sct<-as.matrix(GetAssayData(T15_int,assay="SCT"))
T15_iono_yesFOXI1<-subset(T15_int,cells=rownames(T15_int@meta.data)[which(T15_int@meta.data$clusters_16a == "c6a" & T15_sct["FOXI1",] != 0)])
T15_tuft_yesFOXI1<-subset(T15_int,cells=rownames(T15_int@meta.data)[which(T15_int@meta.data$clusters_16a == "c6b" & T15_sct["FOXI1",] != 0)])




#Get proportions for ionocytes
T15_sct_iono_yesFOXI1<-as.matrix(GetAssayData(T15_iono_yesFOXI1,assay="SCT"))
proportionsList<-list()
for(i in 1:1){
	currData<-T15_sct_iono_yesFOXI1[c("CFTR","POU2F3"),]
	proportions<-c()
	proportions<-append(proportions,length(which(currData["CFTR",] != 0 & currData["POU2F3",] == 0)) / ncol(currData))
	proportions<-append(proportions,length(which(currData["POU2F3",] != 0 & currData["CFTR",] == 0)) / ncol(currData))
	proportions<-append(proportions,length(which(currData["CFTR",] != 0 & currData["POU2F3",] != 0)) / ncol(currData))
	proportions<-append(proportions,length(which(currData["CFTR",] == 0 & currData["POU2F3",] == 0)) / ncol(currData))
	names(proportions)<-c("CFTR only","POU2F3 only","Both","Neither")
	proportionsList[[length(proportionsList) + 1]]<-proportions
}

#Plot pie
#pdf("T15_ionocytes_CFTR_piechart_FOXI1+cells.pdf")
for(i in 1:length(proportionsList)){
	pct <- round(proportionsList[[i]],3)*100
	lbls <- paste(pct,"%",sep="") # ad % to labels 
	colors<-c("red","green","black","white")
	pie(proportionsList[[i]],labels = names(proportionsList[[i]]), col=colors, main = names(proportionsList)[[i]])
}




#Get proportions for Tufts
T15_sct_tuft_yesFOXI1<-as.matrix(GetAssayData(T15_tuft_yesFOXI1,assay="SCT"))
proportionsList<-list()
for(i in 1:1){
	currData<-T15_sct_tuft_yesFOXI1[c("CFTR","POU2F3"),]
	proportions<-c()
	proportions<-append(proportions,length(which(currData["CFTR",] != 0 & currData["POU2F3",] == 0)) / ncol(currData))
	proportions<-append(proportions,length(which(currData["POU2F3",] != 0 & currData["CFTR",] == 0)) / ncol(currData))
	proportions<-append(proportions,length(which(currData["CFTR",] != 0 & currData["POU2F3",] != 0)) / ncol(currData))
	proportions<-append(proportions,length(which(currData["CFTR",] == 0 & currData["POU2F3",] == 0)) / ncol(currData))
	names(proportions)<-c("CFTR only","POU2F3 only","Both","Neither")
	proportionsList[[length(proportionsList) + 1]]<-proportions
}

#Plot pie
#pdf("T15_tuft_POU2F3_piechart_FOXI1+cells.pdf")
for(i in 1:length(proportionsList)){
	pct <- round(proportionsList[[i]],3)*100
	lbls <- paste(pct,"%",sep="") # ad % to labels 
	colors<-c("red","green","black","white")
	pie(proportionsList[[i]],labels = names(proportionsList[[i]]), col=colors, main = names(proportionsList)[[i]])
}

































######################### Make box plots in which I look at Tuft, Iono, PNEC, and Basal cell DEGs in each of the populations

#Repeat the above, except using basal markers instead of tuft markers
treat.col<-data.frame("rare_vs_other_basal"=rep("all",nrow(T15_int@meta.data)),row.names=rownames(T15_int@meta.data))
treat.col[,1]<-as.character(treat.col[,1])
treat.col$rare_vs_other_basal[which(T15_int@meta.data$clusters_16a == "c6a")] <- "ion"
treat.col$rare_vs_other_basal[which(T15_int@meta.data$clusters_16a == "c6b")] <- "tft"
treat.col$rare_vs_other_basal[which(T15_int@meta.data$clusters_16a == "c6c")] <- "pnc"
treat.col$rare_vs_other_basal[which(T15_int@meta.data$clusters_16a == "c1" | T15_int@meta.data$clusters_16a == "c2")] <- "basal"
treat.col$rare_vs_other_basal<-factor(treat.col$rare_vs_other,levels=c("tft","ion","pnc","basal","all"))
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="rare_vs_other_basal")

#Read in DEGs for basal cells
conservedDEGs<-read.table("DEG_Tables/T15_1againstAll_DEGs.txt",header=T,stringsAsFactors=F)
c1DEGs<-conservedDEGs[which(conservedDEGs$comparison == "c1"),][1:100,]
c2DEGs<-conservedDEGs[which(conservedDEGs$comparison == "c2"),][1:100,]
conservedDEGs<-rbind(c1DEGs,c2DEGs)

#Bring in DEGs for the rare cells
rareDEGs<-read.table("Subclustering_results/rare3/DEG_Tables/Rare3_DEGs.txt",header=T,stringsAsFactors=F)
iono<-rareDEGs[which(rareDEGs$p_adj_FDR < 0.05 & rareDEGs$comparison == "Ion_vs_rare"),]
tuft<-rareDEGs[which(rareDEGs$p_adj_FDR < 0.05 & rareDEGs$comparison == "Tuft_vs_rare"),]
pnec<-rareDEGs[which(rareDEGs$p_adj_FDR < 0.05 & rareDEGs$comparison == "PNEC_vs_rare"),]

#Create list of genes
genesList<-list("ionoDEGmean"=iono$gene,"pnecDEGmean"=pnec$gene,"tuftDEGmean"=tuft$gene,"basalDEGmean"=conservedDEGs$gene)

#Get mean expression table
#T15_int<-NormalizeData(T15_int,assay="RNA")
mean_exp_df<-CalculateMeanExpression(dataset=as.matrix(GetAssayData(T15_int,assay="SCT")),genesList=genesList)

#Add table to meta.data
T15_int<-AddMetaData(T15_int,metadata=mean_exp_df,col.name=colnames(mean_exp_df))

#Also, box plots
rareMeansDF<-data.frame("exp"=c(T15_int@meta.data$pnecDEGmean,T15_int@meta.data$ionoDEGmean,T15_int@meta.data$tuftDEGmean,T15_int@meta.data$basalDEGmean),
	"cellStatus"=c(rep(T15_int@meta.data$rare_vs_other_basal,4)),
	"category"=c(rep("pnec",nrow(T15_int@meta.data)),rep("iono",nrow(T15_int@meta.data)),rep("tuft",nrow(T15_int@meta.data)),rep("basal",nrow(T15_int@meta.data))),
	stringsAsFactors=F)
catVec<-c("tuft","iono","pnec","basal")
nameVec<-c("Tuft DEGs","Ionocyte DEGs","PNEC DEGs","basal DEGs")
#pdf("Boxplots_rareCellMeanDEGs_andBasalMarkers_acrossGroups.pdf")
dev.new(height=6,width=5)
par(mfrow=c(2,2),bty="n")
for(i in 1:4){
	currData<-rareMeansDF[which(rareMeansDF$category == catVec[i]),]
	boxplot(currData$exp)~as.factor(currData$cellStatus),las=1,col=c("yellow3","yellow4","wheat","cornflowerblue","grey"),
		ylab="Normalized epxression",xaxt="n",outcex=0.7,main=nameVec[i])
	text(x=c(0.5,1.7,2.85,4.0,5.2),y=par()$usr[3]-0.08*(par()$usr[4]-par()$usr[3]),labels=c("tuft","iono","pnec","basal","other"),
		srt=45,adj=0.8,xpd=T)
}























############################ Now make violin plots for DEGs up in each of the rare cells for Supplementary Table 3

#First bring in the rare cell DEGs
T15_sct<-as.matrix(GetAssayData(T15_int,assay="SCT"))
rareDEGs<-read.table("Subclustering_results/rare3/DEG_Tables/Rare3_DEGs.txt",header=T,stringsAsFactors=F)
rareDEGs<-rareDEGs[which(rareDEGs$p_adj_FDR < 0.05),]

#Get groups
iono_rare<-rareDEGs$gene[which(rareDEGs$comparison == "Ion_vs_rare")]
tuft_rare<-rareDEGs$gene[which(rareDEGs$comparison == "Tuft_vs_rare")]
pnec_rare<-rareDEGs$gene[which(rareDEGs$comparison == "PNEC_vs_rare")]
iono_nonrare<-rareDEGs$gene[which(rareDEGs$comparison == "Ion_vs_nonrare")]
tuft_nonrare<-rareDEGs$gene[which(rareDEGs$comparison == "Tuft_vs_nonrare")]
pnec_nonrare<-rareDEGs$gene[which(rareDEGs$comparison == "PNEC_vs_nonrare")]

#Isolated DEGs up in all 3 versus nonrare (66 DEGs)
rare_core<-intersect(pnec_nonrare,intersect(iono_nonrare,tuft_nonrare))

#Isolate DEGs shared between c6a and c6b (to the exclusion of c6c) either versus rare and/or nonrare
iono_tuft<-unique(c(setdiff(intersect(iono_nonrare,tuft_nonrare),pnec_nonrare),setdiff(intersect(iono_rare,tuft_rare),pnec_rare)))

#Isolate DEGs shared between c6a and c6c (to the exclusion of c6b) either versus rare and/or nonrare
iono_pnec<-unique(c(setdiff(intersect(iono_nonrare,pnec_nonrare),tuft_nonrare),setdiff(intersect(iono_rare,pnec_rare),tuft_rare)))

#Isolate DEGs shared between c6b and c6c (to the exclusion of c6a) either versus rare and/or nonrare
tuft_pnec<-unique(c(setdiff(intersect(tuft_nonrare,pnec_nonrare),iono_nonrare),setdiff(intersect(tuft_rare,pnec_rare),iono_rare)))

#Now get means of these genes
mean_exp_df<-CalculateMeanExpression(dataset=T15_sct,genesList=list(rare_core,iono_tuft,iono_pnec,tuft_pnec),method="geometric")
colnames(mean_exp_df)<-c("rare_core","iono_tuft","iono_pnec","tuft_pnec")

#Add table to meta.data
T15_int<-AddMetaData(T15_int,metadata=mean_exp_df,col.name=colnames(mean_exp_df))


#Now allocate clusters of interest
treat.col<-data.frame("Rare_groups"=T15_int@meta.data$clusters_16a,row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col$Rare_groups[which(T15_int@meta.data$clusters_16a %in% c("c6a"))] <-"Ionocytes"
treat.col$Rare_groups[which(T15_int@meta.data$clusters_16a %in% c("c6b"))] <-"Tuft-like"
treat.col$Rare_groups[which(T15_int@meta.data$clusters_16a %in% c("c6c"))] <-"PNEC"
treat.col$Rare_groups[which(T15_int@meta.data$clusters_16a %in% grep("c6",unique(T15_int@meta.data$clusters_16a),value=T,invert=T))] <-"Other"
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="Rare_groups")
T15_int@meta.data$Rare_groups<-factor(T15_int@meta.data$Rare_groups,levels=c("Ionocytes","Tuft-like","PNEC","Other"))




#Make violin plots
#pdf(Rare_core_violinPlots.pdf")
dev.new(height=3.5,width=4)
vioplot(T15_int@meta.data$rare_core[which(T15_int@meta.data$Rare_groups == "Ionocytes")],
	T15_int@meta.data$rare_core[which(T15_int@meta.data$Rare_groups == "Tuft-like")],
	T15_int@meta.data$rare_core[which(T15_int@meta.data$Rare_groups == "PNEC")],
	T15_int@meta.data$rare_core[which(T15_int@meta.data$Rare_groups == "Other")],
	drawRect=F,las=1,at=1:4,names=rep("",4),col=c("grey"))
#Now add in beeswarm
#Get jitter coordinates
jitterCoordinates<-jitter(as.numeric(T15_int@meta.data$Rare_groups[which(T15_int@meta.data$Rare_groups != "Other")]),amount=0.3)
points(T15_int@meta.data$rare_core[which(T15_int@meta.data$Rare_groups != "Other")]~jitterCoordinates,col="black",bg="black",pch=21,cex=0.6,lwd=0.5)



#pdf(Iono_tuft_violinPlots.pdf")
dev.new(height=3.5,width=4)
vioplot(T15_int@meta.data$iono_tuft[which(T15_int@meta.data$Rare_groups == "Ionocytes")],
	T15_int@meta.data$iono_tuft[which(T15_int@meta.data$Rare_groups == "Tuft-like")],
	T15_int@meta.data$iono_tuft[which(T15_int@meta.data$Rare_groups == "PNEC")],
	T15_int@meta.data$iono_tuft[which(T15_int@meta.data$Rare_groups == "Other")],
	drawRect=F,las=1,at=1:4,names=rep("",4),col=c("grey"))
#Now add in beeswarm
#Get jitter coordinates
jitterCoordinates<-jitter(as.numeric(T15_int@meta.data$Rare_groups[which(T15_int@meta.data$Rare_groups != "Other")]),amount=0.3)
points(T15_int@meta.data$iono_tuft[which(T15_int@meta.data$Rare_groups != "Other")]~jitterCoordinates,col="black",bg="black",pch=21,cex=0.6,lwd=0.5)



#pdf(Iono_pnec_violinPlots.pdf")
dev.new(height=3.5,width=4)
vioplot(T15_int@meta.data$iono_pnec[which(T15_int@meta.data$Rare_groups == "Ionocytes")],
	T15_int@meta.data$iono_pnec[which(T15_int@meta.data$Rare_groups == "Tuft-like")],
	T15_int@meta.data$iono_pnec[which(T15_int@meta.data$Rare_groups == "PNEC")],
	T15_int@meta.data$iono_pnec[which(T15_int@meta.data$Rare_groups == "Other")],
	drawRect=F,las=1,at=1:4,names=rep("",4),col=c("grey"))
#Now add in beeswarm
#Get jitter coordinates
jitterCoordinates<-jitter(as.numeric(T15_int@meta.data$Rare_groups[which(T15_int@meta.data$Rare_groups != "Other")]),amount=0.3)
points(T15_int@meta.data$iono_pnec[which(T15_int@meta.data$Rare_groups != "Other")]~jitterCoordinates,col="black",bg="black",pch=21,cex=0.6,lwd=0.5)



#pdf(Tuft_pnec_violinPlots.pdf")
dev.new(height=3.5,width=4)
vioplot(T15_int@meta.data$tuft_pnec[which(T15_int@meta.data$Rare_groups == "Ionocytes")],
	T15_int@meta.data$tuft_pnec[which(T15_int@meta.data$Rare_groups == "Tuft-like")],
	T15_int@meta.data$tuft_pnec[which(T15_int@meta.data$Rare_groups == "PNEC")],
	T15_int@meta.data$tuft_pnec[which(T15_int@meta.data$Rare_groups == "Other")],
	drawRect=F,las=1,at=1:4,names=rep("",4),col=c("grey"))
#Now add in beeswarm
#Get jitter coordinates
jitterCoordinates<-jitter(as.numeric(T15_int@meta.data$Rare_groups[which(T15_int@meta.data$Rare_groups != "Other")]),amount=0.3)
points(T15_int@meta.data$tuft_pnec[which(T15_int@meta.data$Rare_groups != "Other")]~jitterCoordinates,col="black",bg="black",pch=21,cex=0.6,lwd=0.5)





########### Do enrichments for the four lists
rare_suppViolin_forEnrichr<-data.frame("gene"=c(rare_core,iono_tuft,iono_pnec,tuft_pnec),"comparison"=c(rep("rare_core",length(rare_core)),
	rep("iono_tuft",length(iono_tuft)),rep("iono_pnec",length(iono_pnec)),rep("tuft_pnec",length(tuft_pnec))))
write.table(rare_suppViolin_forEnrichr,quote=F,row.names=F,file="Enrichr/rare_suppViolin_forEnrichr.txt")

#Do enrichments
DEG_table<-read.table("Enrichr/rare_suppViolin_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"rare_suppViolin"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)

























############################ Make violin plots for Figure S6a with Montoro rare cell markers and TFs
#Read in montoro markers
rare_mont<-readExcelTabs("Gene_lists/Monotoro_rareCell_markersAndTFs.xlsx")

#Convert mouse genes to human ones
rare_mont<-lapply(rare_mont,function(x)convertMouseGeneList(x$gene))

#Make new metadata that allocates populations to rare and non-rare cells
treat.col<-data.frame("RareGroups"=T15_int@meta.data$clusters_16a,row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col$RareGroups[which(T15_int@meta.data$clusters_10 %in% c("c1","c1c","c2","c3","c4","c5","c7","c8"))] <-"Other"
treat.col$RareGroups[which(T15_int@meta.data$clusters_16a %in% c("c6a"))] <-"Ionocytes"
treat.col$RareGroups[which(T15_int@meta.data$clusters_16a %in% c("c6b"))] <-"Tuft-like"
treat.col$RareGroups[which(T15_int@meta.data$clusters_16a %in% c("c6c"))] <-"PNEC"
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="RareGroups")
T15_int@meta.data$RareGroups<-factor(T15_int@meta.data$RareGroups,levels=c("Ionocytes","Tuft-like","PNEC","Other"))

#Violin plot dataframe
T15_sct<-as.matrix(GetAssayData(T15_int,assay="SCT"))

#Fill in list with rare expression tables
rare_mont_means_list<-list()
for(i in 1:length(rare_mont)){
	rare_mont_means<-data.frame()
	for(j in 1:length(levels(T15_int@meta.data$RareGroups))){
		if(j==1){
			rare_mont_means<-cbind(CalculateMeanExpression(T15_sct[,which(T15_int@meta.data$RareGroups == 
				levels(T15_int@meta.data$RareGroups)[j])],list(rare_mont[[i]]),method = "geometric"),
				"comparison"=rep(levels(T15_int@meta.data$RareGroups)[j],length(which(T15_int@meta.data$RareGroups == 
				levels(T15_int@meta.data$RareGroups)[j]))))
		}else{
			rare_mont_means<-rbind(rare_mont_means,cbind(CalculateMeanExpression(T15_sct[,which(T15_int@meta.data$RareGroups == 
				levels(T15_int@meta.data$RareGroups)[j])],list(rare_mont[[i]]),method = "geometric"),
				"comparison"=rep(levels(T15_int@meta.data$RareGroups)[j],length(which(T15_int@meta.data$RareGroups == 
				levels(T15_int@meta.data$RareGroups)[j])))))
		}
	}
	names(rare_mont_means)[1]<-"mean_exp"
	rare_mont_means_list[[length(rare_mont_means_list) + 1]]<-rare_mont_means
	names(rare_mont_means_list)[length(rare_mont_means_list)]<-names(rare_mont)[i]
}
		
#Now make violin plots for each of the six gene lists
library(ggplot2)
p1 <- ggplot(rare_mont_means_list[[1]], aes(x=comparison, y=mean_exp, fill=comparison)) + geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) + 
	scale_fill_manual(values=c("gold4","yellow3","wheat","grey")) + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + 
	ylab("Normalized expression") + ggtitle("Ionocyte markers")
p2 <- ggplot(rare_mont_means_list[[2]], aes(x=comparison, y=mean_exp, fill=comparison)) + geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) + 
	scale_fill_manual(values=c("gold4","yellow3","wheat","grey")) + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + 
	ylab("Normalized expression") + ggtitle("Tuft-like markers")
p3 <- ggplot(rare_mont_means_list[[3]], aes(x=comparison, y=mean_exp, fill=comparison)) + geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) + 
	scale_fill_manual(values=c("gold4","yellow3","wheat","grey")) + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + 
	ylab("Normalized expression") + ggtitle("PNEC markers")
p4 <- ggplot(rare_mont_means_list[[4]], aes(x=comparison, y=mean_exp, fill=comparison)) + geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) + 
	scale_fill_manual(values=c("gold4","yellow3","wheat","grey")) + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + 
	ylab("Normalized expression") + ggtitle("Ionocyte TFs")
p5 <- ggplot(rare_mont_means_list[[5]], aes(x=comparison, y=mean_exp, fill=comparison)) + geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) + 
	scale_fill_manual(values=c("gold4","yellow3","wheat","grey")) + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + 
	ylab("Normalized expression") + ggtitle("Tuft-like TFs")
p6 <- ggplot(rare_mont_means_list[[6]], aes(x=comparison, y=mean_exp, fill=comparison)) + geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) + 
	scale_fill_manual(values=c("gold4","yellow3","wheat","grey")) + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + 
	ylab("Normalized expression") + ggtitle("PNEC TFs")

#pdf("T15_violinPlots_rareCells_MonotroMarkers.pdf")	
dev.new(width=10,height=5)
par(mfrow=c(2,3))
plot_grid(p1,p2,p3,p4,p5,p6)
	





########################## Also, make a table of the Montoro genes that describes 
#1) is the gene in our dataset, 
#2) is it a DEG (compared to the rest of the epithelium) in the right population (expected pattern)
#3) is it a DEG in the right population plus one other (expected + 1 other)
#4) is it a DEG in ALL the rare populations (shared in all)
#5) is it a DEG (compared to the rest of the epithelum) in the right population, but a rare cell DEG (compared to the target cell) in the wrong rare cell as well?

#Bring in the rare cell DEGs
rare_DEGs<-read.table("Subclustering_results/rare3/DEG_Tables/Rare3_DEGs.txt",header=T,stringsAsFactors=F)
rare_DEGs<-rare_DEGs[which(rare_DEGs$p_adj_FDR < 0.05),]

#Get vector of cells 
currVec<-list(c("Ion_vs_nonrare","Tuft_vs_nonrare","PNEC_vs_nonrare"),c("Tuft_vs_nonrare","Ion_vs_nonrare","PNEC_vs_nonrare"),
	c("PNEC_vs_nonrare","Ion_vs_nonrare","Tuft_vs_nonrare"),c("Ion_vs_nonrare","Tuft_vs_nonrare","PNEC_vs_nonrare"),
	c("Tuft_vs_nonrare","Ion_vs_nonrare","PNEC_vs_nonrare"),c("PNEC_vs_nonrare","Ion_vs_nonrare","Tuft_vs_nonrare"))
currVec_rare<-list(c("Ion_vs_rare","Tuft_vs_rare","PNEC_vs_rare"),c("Tuft_vs_rare","Ion_vs_rare","PNEC_vs_rare"),
	c("PNEC_vs_rare","Ion_vs_rare","Tuft_vs_rare"),c("Ion_vs_rare","Tuft_vs_rare","PNEC_vs_rare"),
	c("Tuft_vs_rare","Ion_vs_rare","PNEC_vs_rare"),c("PNEC_vs_rare","Ion_vs_rare","Tuft_vs_rare"))

#Create datasets
rare_mont_list<-list()
for(i in 1:length(rare_mont)){
	rare_mont_list[[length(rare_mont_list) + 1]]<-data.frame("Gene"=rare_mont[[i]],"Expressed"=FALSE,"Expected_pattern"=FALSE,
		"Expected_pattern_plus1"=FALSE,"Shared_in_all"=FALSE,"Different_order"=FALSE)
	
	#Get those expressed
	rare_mont_list[[length(rare_mont_list)]]$Expressed<-rare_mont_list[[length(rare_mont_list)]]$Gene %in% rownames(T15_sct)
	#Get those that are only up in the target
	rare_mont_list[[length(rare_mont_list)]]$Expected_pattern<-rare_mont_list[[length(rare_mont_list)]]$Gene %in% 
		setdiff(rare_DEGs$gene[grep(currVec[[i]][1],rare_DEGs$comparison)],c(rare_DEGs$gene[grep(currVec[[i]][2],rare_DEGs$comparison)],
		rare_DEGs$gene[grep(currVec[[i]][3],rare_DEGs$comparison)]))
	
	#Get those shared between target and other 1
	shared1<-c(rare_DEGs$gene[grep(currVec[[i]][1],rare_DEGs$comparison)],rare_DEGs$gene[grep(currVec[[i]][2],rare_DEGs$comparison)])[which(
			duplicated(c(rare_DEGs$gene[grep(currVec[[i]][1],rare_DEGs$comparison)],rare_DEGs$gene[grep(currVec[[i]][2],rare_DEGs$comparison)])))]
	#Get those shared between target and other 2
	shared2<-c(rare_DEGs$gene[grep(currVec[[i]][1],rare_DEGs$comparison)],rare_DEGs$gene[grep(currVec[[i]][3],rare_DEGs$comparison)])[which(
			duplicated(c(rare_DEGs$gene[grep(currVec[[i]][1],rare_DEGs$comparison)],rare_DEGs$gene[grep(currVec[[i]][3],rare_DEGs$comparison)])))]	
	#Get genes only in 1 of these
	rare_mont_list[[length(rare_mont_list)]]$Expected_pattern_plus1<-rare_mont_list[[length(rare_mont_list)]]$Gene %in% 
		c(setdiff(shared1,shared2),setdiff(shared2,shared1))
	#Get genes shared by all three
	rare_mont_list[[length(rare_mont_list)]]$Shared_in_all<-rare_mont_list[[length(rare_mont_list)]]$Gene %in% 
		intersect(shared1,shared2)
	#Get genes up in the target compared to non-rare cells, but a DEG (when compared to the rare cells) in a non-target, but not the target
	rare_mont_list[[length(rare_mont_list)]]$Different_order<-rare_mont_list[[length(rare_mont_list)]]$Gene %in% 
		setdiff(intersect(rare_DEGs$gene[grep(currVec[[i]][1],rare_DEGs$comparison)],
		unique(c(rare_DEGs$gene[grep(currVec_rare[[i]][2],rare_DEGs$comparison)],rare_DEGs$gene[grep(currVec_rare[[i]][3],rare_DEGs$comparison)]))),
		rare_DEGs$gene[grep(currVec_rare[[i]][1],rare_DEGs$comparison)])
}
names(rare_mont_list)<-names(rare_mont)

#Write to file
write.xlsx(rare_mont_list,"T15_overlap_rareCells_MontoroMarkers.xlsx",firstRow=T,colWidths="auto")




