#bsub -e test.err -o test.out -R "rusage[mem=64000]" "module load R/3.5.1 && env R_MAX_VSIZE=700Gb R CMD BATCH subcluster_smg_noC1c_reintegrate_allGenes.r"
library(Seurat, lib.loc="/Seibold/home/jacksonna/R-dev")
library(plotrix)
library(scales)
library(ggplot2)
library(cowplot)
library(openxlsx)


#Download the processed Seurat R object, "Processed_invitro_seurat.Rdata" from GEO,
#obtainable from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134174
#Then load the object it contains ("T15_int") into the session.
#T15_int is a Seurat dataset that contains the processed, integrated dataset
load("Processed_invitro_seurat.Rdata")



#Get subset
Idents(T15_int) <- T15_int@meta.data$clusters_10
T15_3_7<-subset(T15_int,idents=c("c3","c7"))

#To construct a reference, we will identify ‘anchors’ between the individual datasets. First, we split the combined object into a list, with each dataset as an element.
T15.list<-SplitObject(T15_3_7,split.by="donor")

#Do Seurat's scTransform plug in in liu of the normal normalization and feature selection 
#NOTE that with SCTransform::vst, the default is to skip genes not expressed in at least 5 cells (min_cells = 5), as these rare genes violate the
#negative binomial. However, this seems to strict for when subsetting, as I only have a few cells per donor. To keep as many genes as possible
#(since I can only integrate across the genes that intersect across ALL donors), I'm moving min_cells to 2
for (i in 1:length(T15.list)) {
    T15.list[[i]] <- SCTransform(T15.list[[i]], verbose = FALSE, return.only.var.genes = FALSE, min_cells = 2)
}

#Next, select features for downstream integration, and run PrepSCTIntegration, which ensures that 
#all necessary Pearson residuals have been calculated
T15.features <- SelectIntegrationFeatures(object.list = T15.list, nfeatures = 3000)
T15.list <- PrepSCTIntegration(object.list = T15.list, anchor.features = T15.features, verbose = TRUE)


#Now find integration "anchors" (uses CCA) based on 30 dimensions
#Note that the default k.filter is 200, but this number can't be larger than the smallest number of cells per donor
#Check: table(T15_int@meta.data$donor[which(T15_int@meta.data$clusters_10=="c1c" | T15_int@meta.data$clusters_10=="c3" | T15_int@meta.data$clusters_10=="c7")])
T15.anchors <- FindIntegrationAnchors(T15.list, anchor.features = T15.features, dims = 1:30, normalization.method = "SCT", k.filter = 60)

save(list=c("T15.anchors"),file="/Seibold/proj/Single_Cell_Seq/10X_8Tracheal_donors/ANALYSIS_NJackson/10X_15Tracheal_donors_culled_FINAL/subcluster_smg_noC1c_reintegrate_allGenes.rda")

#Get the list of genes that intersect across all 15 donors
featuresToIntegrate <- Reduce(intersect,lapply(T15.list,function(x)rownames(GetAssayData(x,assay="SCT",slot="data"))))

#Now integrate the datasets
T15_smg <- IntegrateData(anchorset = T15.anchors, dims = 1:30, normalization.method = "SCT",features.to.integrate=featuresToIntegrate)

#Switch to integrated assay. The variable features of this assay are automatically set during IntegrateData
DefaultAssay(T15_smg) <- "integrated"

#Run pca
T15_smg<- RunPCA(T15_smg, npcs = 30, verbose = TRUE)

save(list=c("T15_smg"),file="/Seibold/proj/Single_Cell_Seq/10X_8Tracheal_donors/ANALYSIS_NJackson/10X_15Tracheal_donors_culled_FINAL/subcluster_smg_noC1c_reintegrate_allGenes.rda")

if(".RData" %in% system("ls -a",intern=T)){
	unlink(".RData")
}
