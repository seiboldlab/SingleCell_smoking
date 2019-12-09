#########################################
#############Monocle2 analysis#############
#########################################

#################

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

#Bring in subsetted Seurat object
load("../../Subclustering_reintegrated_rda/T15_smg_noC1c_allGenes.rda")

#Subset data such that myoepithelial cells are removed
T15_smg<-subset(T15_smg,cells=rownames(T15_smg@meta.data)[-which(T15_smg@meta.data$subclusters_4 == "c3c")])
T15_smg<-ScaleData(T15_smg, assay="integrated")
T15_smg<-NormalizeData(T15_smg,assay="SCT")
T15_smg<-ScaleData(T15_smg,assay="SCT")

##Bring in SCT and integrated datasets
#SCT dataset
T15_sct<-as.data.frame(as(as.matrix(GetAssayData(T15_smg,assay="SCT",slot="data")),'sparseMatrix'))
#Integrated dataset (use scale.data, but then don't re-scale when ordering cells)
T15_intd<-as.data.frame(as(as.matrix(GetAssayData(T15_smg,assay="integrated",slot="scale.data")),'sparseMatrix'))


#######################Read in attributes table for the cells (based on Seurat metadata) and read dataset into Monocle
#Set up metadata
currMetadata<-T15_smg@meta.data

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



###Toss genes not expressed in at least 1% of cells
#HSMM_cil <- detectGenes(HSMM_cil, min_expr = 0.1)
#expressed_genes<-row.names(subset(fData(HSMM_cil), num_cells_expressed >= 0.01 * ncol(SPDEF_cil_raw)))
#HSMM_cil<-HSMM_cil[expressed_genes,]






#######################Get genes for ordering
#For INT, take a subselection of all the top variable genes (10,354 in this case)
T15_smg<-FindVariableFeatures(T15_smg,assay="integrated",nfeatures=nrow(as.matrix(GetAssayData(T15_smg,assay="integrated"))))
unsup_clustering_genes_int<-VariableFeatures(T15_smg,assay="integrated")
HSMM_int <- setOrderingFilter(HSMM_int, unsup_clustering_genes_int[1:5000])
length(unsup_clustering_genes_int)



















#######################Pseudotime plotting (INT)
ordering_genes_int <-unsup_clustering_genes_int[1:1000]
HSMM_int <- setOrderingFilter(HSMM_int, ordering_genes_int)
#~30 min on f00
HSMM_int <- reduceDimension(HSMM_int, max_components=2,reduction_method="DDRTree",norm_method="none", pseudo_expr = 0,
	relative_expr=F,scaling=F)
#~15 min on f00
HSMM_int <- orderCells(HSMM_int,reverse=T)







#Plot trajectories
pdf('cellTrajectory_int_topGenes_1000_monocle.pdf')
p1<-plot_cell_trajectory(HSMM_int, color_by = "subclusters_4",show_branch_points = F,cell_size=0.1)+ 
	scale_colour_manual(values=c("slategrey","aquamarine","saddlebrown"))
	
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
pdf('cellTrajectory_smg_monocle.pdf')
p1<-plot_cell_trajectory(HSMM_sct, color_by = "subclusters_4",show_branch_points = F,cell_size=0.1)+ 
	scale_colour_manual(values=c("slategrey","aquamarine","saddlebrown"))
	
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

#Tossing 4,124 genes not expressed in at least 1% of cells (leaving 14,771 genes to test)
#Select the genes to test
to_be_tested<-rownames(as.matrix(exprs(HSMM_sct)[which(rowSums(as.matrix(exprs(HSMM_sct)) != 0) / ncol(as.matrix(exprs(HSMM_sct))) > 0.01),]))
#to_be_tested <- row.names(subset(fData(HSMM_sct),gene_short_name %in% rownames(clustering_DEG_genes_sct[which(clustering_DEG_genes_sct$qval < 0.05),])))
cds_subset <- HSMM_sct[to_be_tested,]

#Now do the DE test, with smoke as a covariate to subtract off its effects  
diff_test_sct_regressSmoke <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime, df=3) + smoke",
	reducedModelFormulaStr="~smoke", cores = detectCores())
diff_test_sct_regressSmoke<-diff_test_sct_regressSmoke[order(diff_test_sct_regressSmoke$qval),]

write.table(diff_test_sct_regressSmoke,row.names=F,sep="\t",quote=F,file="PseudotimeGenes_sct_regressSmoke.txt")



















########## Plot TFs across pseudotime

######### Now make a smooth curve pseudotime plot with genes from each heat map cluster overlaid (scale between 0 and 1)
library(plyr)
cds_subsets<-list(HSMM_sct[c("SOX9")],HSMM_sct[c("FOXC1")],HSMM_sct[c("KRT14")],HSMM_sct[c("KRT5")],HSMM_sct[c("IL33")],HSMM_sct[c("FOSL1")],
	HSMM_sct[c("PTTG1")],HSMM_sct[c("NOTCH1")],HSMM_sct[c("HEY1")],HSMM_sct[c("NOTCH3")],HSMM_sct[c("GRHL1")],
	HSMM_sct[c("ELF3")],HSMM_sct[c("SPDEF")],HSMM_sct[c("BARX2")],HSMM_sct[c("DMBT1")],HSMM_sct[c("SCGB3A1")],HSMM_sct[c("MUC5B")])
pdf("TF_curves_pseudotime_smg_noC1c_scaling_overlaid_SUBSET1.pdf",height = 2.5, width = 4.5)
#dev.new(height = 2.5, width = 4.5)
prepPseudotimePlotDatasets(cds_subsets=cds_subsets,color_by = "subclusters_4",min_expr=0.01,
	color_vec = c("slategrey","turquoise","saddlebrown"),
	nrow = 1,ncol = 1,curve_cols=c("red","darkorange","brown","yellow3","yellowgreen","green",
	"forestgreen","darkturquoise","cornflowerblue","lightblue","blue","purple","magenta","saddlebrown","lightgoldenrod4","black","grey60"),
	fix_y = F,overlay="genes",print_points=F,cell_size = 1.5,relative_expr=T,zeroToOneScale=T)
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
cds_subsets<-list(HSMM_sct[c("SOX9")],HSMM_sct[c("KRT14")],HSMM_sct[c("KRT5")],HSMM_sct[c("IL33")],HSMM_sct[c("TP63")],HSMM_sct[c("FOSL1")],
	HSMM_sct[c("SOX7")],HSMM_sct[c("SOX15")],HSMM_sct[c("PTTG1")],HSMM_sct[c("HMGA1")],HSMM_sct[c("NOTCH1")],HSMM_sct[c("NRARP")],
	HSMM_sct[c("HEY1")],HSMM_sct[c("FOXC1")],HSMM_sct[c("PYCARD")],HSMM_sct[c("LBH")],HSMM_sct[c("NOTCH3")],HSMM_sct[c("ZNF750")],HSMM_sct[c("BARX2")],
	HSMM_sct[c("ELF3")],HSMM_sct[c("PRRX2")],HSMM_sct[c("GRHL1")],HSMM_sct[c("SPDEF")],HSMM_sct[c("CREB3L1")],HSMM_sct[c("ELF5")],
	HSMM_sct[c("SCGB3A1")],HSMM_sct[c("MUC5B")],HSMM_sct[c("DMBT1")])
colorVec<-pData(HSMM_sct)$subclusters_4
colorVec[which(colorVec=="c3a")]<-"slategrey"
colorVec[which(colorVec=="c3b")]<-"aquamarine"
colorVec[which(colorVec=="c7")]<-"saddlebrown"

pdf("TF_curves_pseudotime_smg_individualWithPoints_scaled.pdf")
par(mfrow=c(4,4))
for(i in 1:length(cds_subsets)){
	plotIndividualCurves(cds_subset=cds_subsets[[i]],colorVec=colorVec,zeroToOneScale=T,plot_points=T,plot_interval=F)
}
dev.off()

#png("TF_curves_pseudotime_smg_individualWithPoints_sqrt.png",width=6.5,height=6.5,units="in",res=600)
pdf("TF_curves_pseudotime_smg_individualWithPoints_sqrt.pdf",width=6.5,height=6.5)
#dev.new(width=6.5,height=6.5)
par(mfrow=c(3,5))
for(i in 1:length(cds_subsets)){
	par(mar = c(3,3.5,3,1))
	plotIndividualCurves(cds_subset=cds_subsets[[i]],colorVec=colorVec,zeroToOneScale=F,plot_points=T,plot_interval=F)
}
dev.off()




#Repete for subset
cds_subsets<-list(HSMM_sct[c("SOX9")],HSMM_sct[c("FOXC1")],HSMM_sct[c("KRT14")],HSMM_sct[c("KRT5")],HSMM_sct[c("IL33")],HSMM_sct[c("FOSL1")],
	HSMM_sct[c("PTTG1")],HSMM_sct[c("NOTCH1")],HSMM_sct[c("HEY1")],HSMM_sct[c("NOTCH3")],HSMM_sct[c("GRHL1")],
	HSMM_sct[c("ELF3")],HSMM_sct[c("SPDEF")],HSMM_sct[c("BARX2")],HSMM_sct[c("DMBT1")],HSMM_sct[c("SCGB3A1")],HSMM_sct[c("MUC5B")])
colorVec<-pData(HSMM_sct)$subclusters_4
colorVec[which(colorVec=="c3a")]<-"slategrey"
colorVec[which(colorVec=="c3b")]<-"aquamarine"
colorVec[which(colorVec=="c7")]<-"saddlebrown"

pdf("TF_curves_pseudotime_smg_individualWithPoints_scaled_SUBSET1.pdf")
par(mfrow=c(4,4))
for(i in 1:length(cds_subsets)){
	plotIndividualCurves(cds_subset=cds_subsets[[i]],colorVec=colorVec,zeroToOneScale=T,plot_points=T,plot_interval=F)
}
dev.off()

#png("TF_curves_pseudotime_smg_individualWithPoints_sqrt.png",width=6.5,height=6.5,units="in",res=600)
pdf("TF_curves_pseudotime_smg_individualWithPoints_sqrt_SUBSET1.pdf",width=6.5,height=6.5)
#dev.new(width=6.5,height=6.5)
par(mfrow=c(3,5))
for(i in 1:length(cds_subsets)){
	par(mar = c(3,3.5,3,1))
	plotIndividualCurves(cds_subset=cds_subsets[[i]],colorVec=colorVec,zeroToOneScale=F,plot_points=T,plot_interval=F)
}
dev.off()


