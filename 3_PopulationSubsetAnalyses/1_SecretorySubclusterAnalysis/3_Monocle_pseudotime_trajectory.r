#########################################
#############Monocle2 analysis#############
#########################################

library(data.table)
library(heatmap3)
library(plotrix)
library(devtools)
library(monocle)
library(Seurat, lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library2")
library(dplyr)
library(Matrix)
library(cowplot)
library(openxlsx)
source("../../CommonFunctions.r")


#######################Get data

#Bring in the subsetted Seurat object
attach("../../Subclustering_reintegrating_rda/T15_lumSec_3000.rda")

#Bring in both the SCT and integrated datasets
#SCT dataset 
T15_sct<-as.data.frame(as(as.matrix(GetAssayData(T15_lumSec,assay="SCT",slot="data")),'sparseMatrix'))
#Integrated dataset (use scale.data, but then don't re-scale when ordering cells)
T15_intd<-as.data.frame(as(as.matrix(GetAssayData(T15_lumSec,assay="integrated",slot="scale.data")),'sparseMatrix'))


#######################Read in attributes table for the cells (based on Seurat metadata) and read dataset into Monocle
#Set up metadata
currMetadata<-T15_lumSec@meta.data

#Make sample sheet
sample_sheet = currMetadata

#Make SCT object
gene_names = data.frame("gene_short_name"=rownames(T15_sct),stringsAsFactors=F,row.names=rownames(T15_sct))
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_names) 
HSMM_sct<-newCellDataSet(as.matrix(T15_sct),phenoData = pd,featureData = fd, expressionFamily=negbinomial())

#Make integrated object
gene_names = data.frame("gene_short_name"=rownames(T15_intd),stringsAsFactors=F,row.names=rownames(T15_intd))
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_names) 
HSMM_int<-newCellDataSet(as.matrix(T15_intd),phenoData = pd,featureData = fd, expressionFamily=negbinomial())

###Means and dispersions	
HSMM_sct <- estimateSizeFactors(HSMM_sct)
HSMM_sct <- estimateDispersions(HSMM_sct)
HSMM_int <- estimateSizeFactors(HSMM_int)
HSMM_int <- estimateDispersions(HSMM_int)








#######################Get genes for ordering
#For INT, just use all the genes, since there are only 3000
unsup_clustering_genes_int<-VariableFeatures(T15_lumSec,assay="integrated")
HSMM_int <- setOrderingFilter(HSMM_int, unsup_clustering_genes_int)
length(unsup_clustering_genes_int)













#######################Pseudotime plotting (INT)
#For ordering based on the most differentially expressed genes among the
#ordering_genes_int <-row.names(clustering_DEG_genes_int)[order(clustering_DEG_genes_int$qval)][1:2500]
#ordering_genes_int <- row.names (subset(clustering_DEG_genes_int, qval < 0.01))
ordering_genes_int <-unsup_clustering_genes_int
HSMM_int <- setOrderingFilter(HSMM_int, ordering_genes_int)
#~30 min on f00
HSMM_int <- reduceDimension(HSMM_int, max_components=2,reduction_method="DDRTree",norm_method="none", pseudo_expr = 0,
	relative_expr=F)
#~15 min on f00
HSMM_int <- orderCells(HSMM_int,reverse=F)





#Plot trajectories
pdf('cellTrajectory_int_topDEGs_2500_monocle.pdf')
p1<-plot_cell_trajectory(HSMM_int, color_by = "subclusters_10",show_branch_points = F,cell_size=0.1)+ 
	scale_colour_manual(values=c("lightcoral","deeppink","deeppink4","rosybrown","gold","pink","black","darkgoldenrod","darkorange2","brown"))
	
p2<-plot_cell_trajectory(HSMM_int, color_by="Pseudotime",show_branch_points = F,cell_size=0.1) + 
	scale_colour_gradient(low="snow2",high="springgreen4")

p3<-plot_cell_trajectory(HSMM_int, color_by="donor",show_branch_points = FALSE,cell_size=0.1) + 
	scale_colour_manual(values=c("midnightblue","yellow","orange","pink","tan","darkgreen",
	"violetred4","slateblue4","blue","black","orangered4","yellowgreen","greenyellow","wheat","darkslategray"))

p4<-plot_cell_trajectory(HSMM_int, color_by="smoke",show_branch_points = FALSE,cell_size=0.1) + 
	scale_colour_manual(values=c("red","orange","black"))

p5<-plot_cell_trajectory(HSMM_int, color_by="State",show_branch_points = FALSE,cell_size=0.1) + 
	scale_colour_manual(values=c("red3","cornflowerblue","orange","turquoise","midnightblue","violet","green3","yellow","darkgreen","purple","black"))

plot_grid(p1,p2,p3,p5)
dev.off()











#Now, re-set the root to 5
HSMM_int <- orderCells(HSMM_int,root_state=5)












######Inserting the integrated trajectory into the sct dataset
HSMM_sct@reducedDimS<-HSMM_int@reducedDimS
HSMM_sct@reducedDimW<-HSMM_int@reducedDimW
HSMM_sct@reducedDimK<-HSMM_int@reducedDimK
HSMM_sct@minSpanningTree<-HSMM_int@minSpanningTree
HSMM_sct@cellPairwiseDistances<-HSMM_int@cellPairwiseDistances
HSMM_sct@phenoData@data$Pseudotime<-HSMM_int@phenoData@data$Pseudotime
HSMM_sct@phenoData@data$State<-HSMM_int@phenoData@data$State
HSMM_sct@phenoData@varMetadata<-HSMM_int@phenoData@varMetadata
HSMM_sct@dim_reduce_type<-HSMM_int@dim_reduce_type























#Re-plot trajectories using the SCT dataset, just to make sure the embedding transferred ok
pdf('cellTrajectory_lumSec_monocle.pdf')
p1<-plot_cell_trajectory(HSMM_sct, color_by = "subclusters_10",show_branch_points = F,cell_size=0.1)+ 
	scale_colour_manual(values=c("lightcoral","deeppink","deeppink4","rosybrown","gold","pink","black","darkgoldenrod","darkorange2","brown"))
	
p2<-plot_cell_trajectory(HSMM_sct, color_by="Pseudotime",show_branch_points = F,cell_size=0.1) + 
	scale_colour_gradient(low="snow2",high="springgreen4")

p3<-plot_cell_trajectory(HSMM_sct, color_by="donor",show_branch_points = FALSE,cell_size=0.1) + 
	scale_colour_manual(values=c("midnightblue","yellow","orange","pink","tan","darkgreen",
	"violetred4","slateblue4","blue","black","orangered4","yellowgreen","greenyellow","wheat","darkslategray"))

p4<-plot_cell_trajectory(HSMM_sct, color_by="smoke",show_branch_points = FALSE,cell_size=0.1) + 
	scale_colour_manual(values=c("red","orange","black"))

p5<-plot_cell_trajectory(HSMM_sct, color_by="State",show_branch_points = FALSE,cell_size=0.1) + 
	scale_colour_manual(values=c("red3","cornflowerblue","orange","turquoise","midnightblue","violet","green3","yellow","darkgreen","purple","black"))

plot_grid(p1,p2,p3,p5)
dev.off()





















#################################### Do differential expression

#Tossing 4,717 genes not expressed in at least 1% of cells (leaving 15,124 genes to test)
#Of all the pseudotime DEGs in the previous analysis, only two were among these, so no loss
#sort(rownames(as.matrix(exprs(HSMM_sct)[which(rowSums(as.matrix(exprs(HSMM_sct)) != 0) / ncol(as.matrix(exprs(HSMM_sct))) < 0.01),])))

#Select the genes to test
to_be_tested<-rownames(as.matrix(exprs(HSMM_sct)[which(rowSums(as.matrix(exprs(HSMM_sct)) != 0) / ncol(as.matrix(exprs(HSMM_sct))) > 0.01),]))
#to_be_tested <- row.names(subset(fData(HSMM_sct),gene_short_name %in% rownames(clustering_DEG_genes_sct[which(clustering_DEG_genes_sct$qval < 0.05),])))
cds_subset <- HSMM_sct[to_be_tested,]

#Now do the DE test, with smoke as a covariate to subtract off its effects  
diff_test_sct_regressSmoke <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime, df=3) + smoke",
	reducedModelFormulaStr="~smoke", cores = detectCores())
diff_test_sct_regressSmoke<-diff_test_sct_regressSmoke[order(diff_test_sct_regressSmoke$qval),]

write.table(diff_test_sct_regressSmoke,row.names=F,sep="\t",quote=F,file="PseudotimeGenes_sct_regressSmoke.txt")

#Repeat, without regressing out smoke
diff_test_sct <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)" cores = detectCores())
diff_test_sct<-diff_test_sct[order(diff_test_sct$qval),]

write.table(diff_test_sct,row.names=F,sep="\t",quote=F,file="PseudotimeGenes_sct.txt")


















########## Plot TFs across pseudotime

######### Now make a smooth curve pseudotime plot with genes from each heat map cluster overlaid (scale between 0 and 1)
library(plyr)
cds_subsets<-list(HSMM_sct[c("TP63")],HSMM_sct[c("NOTCH1")],HSMM_sct[c("NOTCH3")],HSMM_sct[c("HES1")],
	HSMM_sct[c("KLF3")],HSMM_sct[c("NKX3-1")],HSMM_sct[c("SPDEF")],HSMM_sct[c("CREB3L1")],
	HSMM_sct[c("XBP1")],HSMM_sct[c("MESP1")],HSMM_sct[c("FOXA3")],HSMM_sct[c("SCGB1A1")],
	HSMM_sct[c("MUC5B")],HSMM_sct[c("MUC5AC")])
pdf("TF_curves_pseudotime_lumSec_scaling_overlaid.pdf",height = 2.5, width = 4.5)
#dev.new(height = 2.5, width = 4.5)
prepPseudotimePlotDatasets(cds_subsets=cds_subsets,color_by = "clusters",min_expr=0.01,
	color_vec = c("slategrey","turquoise","saddlebrown"),
	nrow = 1,ncol = 1,curve_cols=c("red","darkorange","brown","yellow3","yellowgreen","forestgreen","darkturquoise","cornflowerblue",
	"blue","purple","magenta","saddlebrown","black","grey60"),fix_y = F,
	overlay="genes",print_points=F,cell_size = 1.5,relative_expr=T,zeroToOneScale=T)
dev.off()











#Plot TFs across pseudotime one at a time using the Monocle function
genes<-c("TP63","NOTCH1","NOTCH3","HES1","KLF3","NKX3-1","SPDEF","CREB3L1","XBP1","MESP1","FOXA3","SCGB1A1","MUC5B","MUC5AC")
pdf("PsedotimeTFs_oneATATime.pdf")
cds_subset <- HSMM_sct[genes,]
print(plot_genes_in_pseudotime(cds_subset, color_by = "clusters_10",nrow=4,ncol=4,panel_order=genes,min_expr=0.01) + 
	scale_color_manual(values=c("red","darkorange","brown","yellow3","yellowgreen","forestgreen","darkturquoise",
	"cornflowerblue","blue","purple","magenta","saddlebrown","black","grey60")))
dev.off()










#Plot TFs across pseudotime one at a time using tailored function
cds_subsets<-list(HSMM_sct[c("TP63")],HSMM_sct[c("NOTCH1")],HSMM_sct[c("NOTCH3")],HSMM_sct[c("HES1")],
	HSMM_sct[c("KLF3")],HSMM_sct[c("NKX3-1")],HSMM_sct[c("SPDEF")],HSMM_sct[c("CREB3L1")],
	HSMM_sct[c("XBP1")],HSMM_sct[c("MESP1")],HSMM_sct[c("FOXA3")],HSMM_sct[c("SCGB1A1")],
	HSMM_sct[c("MUC5B")],HSMM_sct[c("MUC5AC")])
colorVec<-pData(HSMM_sct)$clusters_10
colorVec[which(colorVec=="c4")]<-"tomato"
colorVec[which(colorVec=="c5")]<-"orange"

pdf("TF_curves_pseudotime_lumSec_individualWithPoints_scaled.pdf")
par(mfrow=c(4,4))
for(i in 1:length(cds_subsets)){
	plotIndividualCurves(cds_subset=cds_subsets[[i]],colorVec=colorVec,zeroToOneScale=T,plot_points=T,plot_interval=F)
}
dev.off()

png("TF_curves_pseudotime_lumSec_individualWithPoints_sqrt.png",width=6.5,height=6.5,units="in",res=600)
#dev.new(width=6.5,height=6.5)
par(mfrow=c(3,5))
for(i in 1:length(cds_subsets)){
	par(mar = c(3,3.5,3,1))
	plotIndividualCurves(cds_subset=cds_subsets[[i]],colorVec=colorVec,zeroToOneScale=F,plot_points=T,plot_interval=F)
}
dev.off()



























######################################################### HEAT MAPS #####################################################

#For plotting heat map of genes, go to separate scripts
####For all 11,000+ pseudotime dependent genes (FDR < 0.05), see Heatmap_PseudotimeGenes.r
####For 5,000+ pseudotime dependent genes (FDR < 1e-10), see Heatmap_topPseudotimeGenes.r
####For for top 3,000 pseudotime dependent genes, see Heatmap_top3000PseudotimeGenes.r

#### The ones with more genes give better enrichments, but yield noisier heat maps. So I'm focusing on the one with 5,000+ DEGs (topDEGs),
#### which compromises between these two factors. These DEGs contain all the TFs in the paper and yield descent enrichments

##########################################################################################################################























############################################ Find cutoff for defining "Mature secretory cells" using pseudotime

####### So we don't need to rely on subclustering to pick out the mature secretory cells, want to isolate those cells
####### most defined by clusters N and O, which are at the end of pseudotime.

###### Looking at the heat map, it looks like average expression of cells will be bimodal (with everthing to the left, yellow or blue and
###### everything to the right, orange or red). However, when I plot out mean expression, it becomes clear that there is
###### an awful lot of noise before and after the point of where we want the threshold to be. Setting a mean expression threshold will either need to
###### be so high as to exclude most cells in the mature group or will need to be so low as to include lots of cells before the mature group.
###### What does appear to be bimodal however is the density of cells along pseudotime, which is low along the branches, and then high at the "end states". 
###### So, what I can do is fit a mixed gaussian model to the data, assuming there is a mixture of two distributions (branches and end states). I can get the
###### mean and standard deviation of the end state distribution and then select all cells within a standard deviation of the end point group mean.

library(mixtools)

#First, bring in the DEGs (topDEGs version)
NO_genes<-read.table("DEG_Tables/PseudotimeGenes_sct_regressSmoke_ordered_16_topDEGs.txt",header=T,stringsAsFactors=F)
NO_genes<-NO_genes$gene[which(NO_genes$cluster == "N" | NO_genes$cluster == "O")]

#Second, bring in the Seurat dataset so I have access to the normalized expression for these genes
attach("../../Subclustering_reintegrating_rda/T15_lumSec_3000.rda")

#Third, bring in the monocle dataset, we we have access to the pseudotime
attach("SAVED_T15_lumSec_monocle.rda")

#Now add the pseudotime as metadata to the T15 Seurat object
all(rownames(pData(HSMM_sct)) == rownames(T15_lumSec@meta.data))
T15_lumSec<-AddMetaData(T15_lumSec,metadata = data.frame(pData(HSMM_sct)$Pseudotime,row.names=rownames(pData(HSMM_sct))),col.name = "Pseudotime")

#Also add Monocle State to the metadata, so I can subset only the final differentiation state (state 1)
T15_lumSec<-AddMetaData(T15_lumSec,metadata = data.frame(pData(HSMM_sct)$State,row.names=rownames(pData(HSMM_sct))),col.name = "Monocle_States")

#Also, add mean expression of these NO genes as metadata
NO_means<-CalculateMeanExpression(dataset=as.matrix(GetAssayData(T15_lumSec,assay="SCT",slot="data")),genesList=list("NO_genes"=NO_genes),method="geometric")
T15_lumSec<-AddMetaData(T15_lumSec,metadata = NO_means,col.name = "NO_means")

#Now plot the mean expression of these genes in order of pseudotime
#Get subcluster_10 colors set up
colorValues<-c("lightcoral","deeppink","deeppink4","rosybrown","gold","pink","black","darkgoldenrod","darkorange2","brown")
colorVec<-T15_lumSec@meta.data$subclusters_10
for(i in 1:length(sort(unique(T15_lumSec@meta.data$subclusters_10)))){
	colorVec[which(colorVec == sort(unique(T15_lumSec@meta.data$subclusters_10))[i])] <- colorValues[i]
}

#Now subselect the cells that are in State 1
pseudotime_subset<-T15_lumSec@meta.data$Pseudotime[which(T15_lumSec@meta.data$Monocle_States == "1")]
NO_means_subset<-T15_lumSec@meta.data$NO_means[which(T15_lumSec@meta.data$Monocle_States == "1")]
colorVec_subset<-colorVec[which(T15_lumSec@meta.data$Monocle_States == "1")]

#Now do scatter plot with loess line (showing high noise in mean expression across pseudotime, but a bimodal distribution of point density between branch and endpoint)
#pdf("ScatterPlot_meanN0_againstPseudotime_state1.pdf")
scatter.smooth(x=pseudotime_subset,y=NO_means_subset,lpars =list(col = "red", lwd = 3, lty = 3),col=colorVec_subset,las=1,
	ylab="Mean expression of NO DEGs",xlab="Pseudotime",bty="l")

#So now plot histogram of pseudotime, where you can see the bimodality
hist(pseudotime_subset, breaks=10)

#In contrast, the distribution of mean expression in unimodal
hist(NO_means_subset, breaks=10)

#So, fit a mixed gaussian model and plot the estimated densities for the two distributions
normMix<-normalmixEM(pseudotime_subset)

#Get means and standard deviations (STATE 1: second mean is 17.87407; second sigma is 0.5247)
#Get means and standard deviations (STATE 5: first mean is 0.990796; first sigma is 0.189020)
summary(normMix)

#plot density
par(mfrow=c(1,2))
plot(normMix,density=T)

plot(density(normMix$posterior[,2] * normMix$x))




######## Plot the histogram and curves myself so I can plot the line where the sigma cutoff will be
#pdf("MixedGaussianCurves_state1.pdf")
par(bty="l")
a <- hist(normMix$x, plot = FALSE)
maxy <- max(max(a$density), 0.3989 * normMix$lambda/normMix$sigma)
hist(normMix$x, prob = TRUE, main = "", xlab = "Pseudotime", xlim= c(min(normMix$x),max(normMix$x)), 
     ylim = c(0, maxy), border = "azure3",las=1)
for (i in 1:ncol(normMix$posterior)) {
  curve(normMix$lambda[i] * dnorm(x, mean = normMix$mu[i], sd = normMix$sigma[i]), 
        col = 1 + i, lwd = 1.5, add = TRUE)
}
#lines(density(normMix$x), lty=2, lwd=0.8)
abline(v=normMix$mu[2] - (normMix$sigma[2] * 2),col="orange",lty=2,lwd=2)
#abline(v=normMix$mu[2] + (normMix$sigma[2] * 3),col="orange",lty=2,lwd=2)
box(lty=1)





########## Finally, export the list of samples in this "mature secretory cell" group
########## defined as those within at least 2 standard deviations of the second distribution
pseudotimeCutoff<-normMix$mu[2] - (normMix$sigma[2] * 2)
matureSecCells<-rownames(pData(HSMM_sct))[which(pData(HSMM_sct)$Pseudotime > pseudotimeCutoff)]
cat(matureSecCells,sep="\n",file="Mature_secretory_cells_monocle.txt")














###### Do differential expression among each of the states in the trajectory

#First, bring in the T15_lumSec Seurat object
attach("../../Subclustering_reintegrating_rda/T15_lumSec_3000.rda")

#Insert pseudotime states into metadata
all(rownames(pData(HSMM_sct)) == rownames(T15_lumSec@meta.data))
T15_lumSec<-AddMetaData(T15_lumSec,metadata = data.frame(pData(HSMM_sct)$State,row.names=rownames(pData(HSMM_sct))),col.name = "Monocle_States")

#Do DE
Idents(T15_lumSec) <- T15_lumSec@meta.data$Monocle_States
clusterList<-sort(as.character(unique(Idents(T15_lumSec))))
markers.list <- list()
for(i in 1:length(clusterList)) {
    markers.list[[i]] <- FindMarkers(T15_lumSec, ident.1=clusterList[i], min.pct=0.1, logfc.threshold=0.25, only.pos=T, 
    	assay="SCT", test.use = "LR", latent.vars = "smoke")
}

#Combine tables
Idents(T15_lumSec) <- T15_lumSec@meta.data$Monocle_States
clusterList<-sort(as.character(unique(Idents(T15_lumSec))))
comparisonVec<-clusterList

#Get vector of clusters for calculating pct values from RNA
pctIndexList<-list()
for(i in 1:length(clusterList)){
	currClusterIndex<-which(clusterList == clusterList[i])
	pctIndexList[[length(pctIndexList) + 1]]<-list(clusterList[i],clusterList[-currClusterIndex])
}

#Get RNA counts and cluster assignments
counts<-as.matrix(GetAssayData(T15_lumSec,assay="RNA",slot="data"))
clusterAssignments<-T15_lumSec@meta.data$Monocle_States

#Run OrganizeDF_list
T15_lumSec_monocleStates_DEGs<-OrganizeDF_list(datasetList=markers.list,comparisonVec=comparisonVec,pctIndexList=pctIndexList,
	counts=counts,clusterAssignments=clusterAssignments)

#Write all to single txt file
write.table(T15_lumSec_monocleStates_DEGs,sep="\t",quote=FALSE,row.names=FALSE,file="DEG_Tables/T15_lumSec_monocleStates_DEGs.txt")

#Write to excel
writeDEGsToExcel(DEG_table = T15_lumSec_monocleStates_DEGs, savedFile = "DEG_Tables/T15_lumSec_monocleStates_DEGs.xlsx")

#Export table for Enrichr
T15_lumSec_monocleStates_DEGs_forEnrichr<-T15_lumSec_monocleStates_DEGs[which(T15_lumSec_monocleStates_DEGs$p_adj_FDR < 0.05),] #Only good DEGs

#Write all to single txt file
write.table(T15_lumSec_monocleStates_DEGs_forEnrichr,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,
	file="Enrichr/T15_lumSec_monocleStates_DEGs_forEnrichr.txt")

#Do enrichments
DEG_table<-read.table("Enrichr/T15_lumSec_monocleStates_DEGs_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"T15_lumSec_monocleStates_DEGs"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)







