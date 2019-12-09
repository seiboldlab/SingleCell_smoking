#bsub -e test.err -o test.out -R "rusage[mem=64000]" "module load R/3.5.1 && env R_MAX_VSIZE=700Gb R CMD BATCH subcluster_c4_5_reintegrate_3000.r"
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
T15_c4_5<-subset(T15_int,idents=c("c4","c5"))

#To construct a reference, we will identify ‘anchors’ between the individual datasets. First, we split the combined object into a list, with each dataset as an element.
T15.list<-SplitObject(T15_c4_5,split.by="donor")

#Do Seurat's scTransform plug in in liu of the normal normalization and feature selection (~20 min on cluster)
for (i in 1:length(T15.list)) {
    T15.list[[i]] <- SCTransform(T15.list[[i]], verbose = FALSE)
}

#Next, select features for downstream integration, and run PrepSCTIntegration, which ensures that 
#all necessary Pearson residuals have been calculated
T15.features <- SelectIntegrationFeatures(object.list = T15.list, nfeatures = 3000)
T15.list <- PrepSCTIntegration(object.list = T15.list, anchor.features = T15.features, verbose = TRUE)


#Now find integration "anchors" (uses CCA) based on 30 dimensions
T15.anchors <- FindIntegrationAnchors(T15.list, anchor.features = T15.features, dims = 1:30, normalization.method = "SCT")

save(list=c("T15.anchors"),file="/Seibold/proj/Single_Cell_Seq/10X_8Tracheal_donors/ANALYSIS_NJackson/10X_15Tracheal_donors_culled_FINAL/subcluster_c4_5_reintegrate_3000.rda")

#Now integrate the datasets
T15_lumSec <- IntegrateData(anchorset = T15.anchors, dims = 1:30, normalization.method = "SCT")

#Switch to integrated assay. The variable features of this assay are automatically set during IntegrateData
DefaultAssay(T15_lumSec) <- "integrated"

#Run pca
T15_lumSec<- RunPCA(T15_lumSec, npcs = 30, verbose = TRUE)

save(list=c("T15.anchors","T15_lumSec"),file="/Seibold/proj/Single_Cell_Seq/10X_8Tracheal_donors/ANALYSIS_NJackson/10X_15Tracheal_donors_culled_FINAL/subcluster_c4_5_reintegrate_3000.rda")
