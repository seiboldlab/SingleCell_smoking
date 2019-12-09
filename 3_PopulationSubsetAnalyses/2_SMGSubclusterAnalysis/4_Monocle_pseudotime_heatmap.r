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

######################################################### MONCLE HEAT MAPS #####################################################

###################################### Now make heat map of pseudotime dependent genes

#bsub -e test.err -o test.out -R "rusage[mem=64000]" "module load R/3.5.1 && env R_MAX_VSIZE=700Gb R CMD BATCH 5_Heatmap_pseudotimeDEGs_lumSec_16.r"

#First, read in the genes
diff_test_sct_regressSmoke<-read.table("PseudotimeGenes_sct_regressSmoke.txt",header=T,stringsAsFactors=F)
genes_subset<-diff_test_sct_regressSmoke$gene_short_name[which(diff_test_sct_regressSmoke$qval < 1e-10)]

#Get annotations for different clusters
annotations1<-getBinnedValues(pData=pData(HSMM_sct),metadata_col_index=c(45))
annotations<-annotations1

#Plot a heat map, where rows are genes and where columns are cells, ordered
#by pseudotime (the add_annotation_col argument doesn't seem to work on monocle 2.4 on the cluster, so need
#to use R/3.5.1 which has monocle 2.8)
pdf("Heatmap_PseudotimeGenes_sct.pdf")
pheatmap_16<-plot_pseudotime_heatmap(HSMM_sct[genes_subset,],
	num_clusters = 16,
	cores = detectCores(),
	show_rownames = F,
	add_annotation_col=annotations,
	return_heatmap = T)
dev.off()


#Print tree
pdf("Heatmap_PseudotimeGenes_sct_16_topDEGs.pdf")
pheatmap_16
dev.off()








####################################### Now modify colors for the heat map
############### Here are where relevant colors are stored

#I     ####### Legend colors ########
##### Here's where the State legend colors are stored
#pheatmap_16[[4]][[1]][[7]]$children[[6]]$gp$fill

#### Here's where the subclusters_4 legend colors are stored
pheatmap_16[[4]][[1]][[7]]$children[[6]]$gp$fill

#### Here's where the gene clusters legend colors are stored
pheatmap_16[[4]][[1]][[7]]$children[[9]]$gp$fill



#II     ######### Annotation bar colors ##########
##### Here's where the column header colors are stored, in order
pheatmap_16[[4]][[1]][[3]]$gp$fill



#III     ########## Gene cluster colors ###########
#### Here's where the gene cluster colors are stored
pheatmap_16[[4]][[1]][[5]]$gp$fill






#########################################################Let's modify the colors first

########### So first, let's change the gene cluster colors
#New colors
colorVec<-c("blue","dodgerblue3","mediumpurple3","lightskyblue","purple","turquoise","slategrey","brown","yellow","deeppink4","pink",
	"yellow4","wheat","azure3","darkolivegreen","forestgreen","darkslategray","springgreen")
#Change the legend first
for(i in 1:length(pheatmap_16[[4]][[1]][[7]]$children[[9]]$gp$fill)){
	pheatmap_16[[4]][[1]][[7]]$children[[9]]$gp$fill[which(pheatmap_16[[4]][[1]][[7]]$children[[9]]$gp$fill == 
		unique(pheatmap_16[[4]][[1]][[5]]$gp$fill)[i])] <- colorVec[i]
}
#Now, let's change the tree colors
for(i in 1:length(unique(pheatmap_16[[4]][[1]][[5]]$gp$fill))){
	pheatmap_16[[4]][[1]][[5]]$gp$fill[which(pheatmap_16[[4]][[1]][[5]]$gp$fill == unique(pheatmap_16[[4]][[1]][[5]]$gp$fill)[i])]<-colorVec[i]
}




######sublusters_4
###First, need to change the annotation bar colors
#New colors
colorVec<-c("slategrey","turquoise","saddlebrown")
for(i in 1:length(pheatmap_16[[4]][[1]][[7]]$children[[6]]$gp$fill)){
	pheatmap_16[[4]][[1]][[3]]$gp$fill[,which(colnames(pheatmap_16[[4]][[1]][[3]]$gp$fill) == "subclusters_4")][
		which(pheatmap_16[[4]][[1]][[3]]$gp$fill[,which(colnames(pheatmap_16[[4]][[1]][[3]]$gp$fill) == "subclusters_4")] ==
		pheatmap_16[[4]][[1]][[7]]$children[[6]]$gp$fill[i])] <- colorVec[i]
}

###Then, change the legend
for(i in 1:length(pheatmap_16[[4]][[1]][[7]]$children[[6]]$gp$fill)){
	pheatmap_16[[4]][[1]][[7]]$children[[6]]$gp$fill[i] <- colorVec[i]
}



#Print tree
pdf("Heatmap_PseudotimeGenes_sct_16_topDEGs.pdf")
pheatmap_16
dev.off()













#########################################################Let's finally modify the order of the genes so progress from left to right
#To do this, we just need to identify the order we want 
cluster_order<-c(16,12,9,5,4,10,14,11,8,15,13,7,2,3,1,6)

genes_subset_clusters<-c()
for(i in 1:length(pheatmap_16[[4]][[1]][[7]]$children[[9]]$gp$fill)){
	#If I want to reverse the order of the genes, then do it
	if(i == 12 | i == 10 | i == 13 | i == 16){
		genes_subset_clusters<-append(genes_subset_clusters,
			rev(rownames(pheatmap_16[[4]][[1]][[5]]$gp$fill)[which(pheatmap_16[[4]][[1]][[5]]$gp$fill == 
			pheatmap_16[[4]][[1]][[7]]$children[[9]]$gp$fill[cluster_order[i]])]))
	}else{ #else don't
		genes_subset_clusters<-append(genes_subset_clusters,
			rownames(pheatmap_16[[4]][[1]][[5]]$gp$fill)[which(pheatmap_16[[4]][[1]][[5]]$gp$fill == 
			pheatmap_16[[4]][[1]][[7]]$children[[9]]$gp$fill[cluster_order[i]])])
	}
}









#################################Now, also create a new pseudotime DEG table that identifies the clusters to which they belong,
#################################keeping with the new ordering
#####Input to get genes for each cluster
cluster2colors<-pheatmap_16[[4]][[1]][[7]]$children[[9]]$gp$fill
gene2colors<-pheatmap_16[[4]][[1]][[5]]$gp$fill

#Get cluster orders
cluster_order<-c(16,12,9,5,4,10,14,11,8,15,13,7,2,3,1,6)
letterVec<-c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P")
mergedVec<-c(rep("A-C",3),rep("D-E",2),rep("F-G",2),rep("H-K",4),rep("L-M",2),rep("N-P",3))

#Prepare dataset
diff_test_sct_regressSmoke<-read.table("DEG_Tables/PseudotimeGenes_sct_regressSmoke.txt",header=T,stringsAsFactors=F)
diff_test_sct_regressSmoke_sig<-diff_test_sct_regressSmoke[which(diff_test_sct_regressSmoke$qval < 1e-10),c(5,3:4)]
colnames(diff_test_sct_regressSmoke_sig)[1]<-"gene"
diff_test_sct_regressSmoke_sig<-cbind(diff_test_sct_regressSmoke_sig,"hclust_num"=NA,"cluster"=NA,"clusters_merged"=NA)

#Assign the clusters
for(i in 1:length(cluster_order)){
	diff_test_sct_regressSmoke_sig$hclust_num[which(diff_test_sct_regressSmoke_sig$gene %in% 
		rownames(gene2colors)[which(gene2colors == cluster2colors[cluster_order[i]])])] <- cluster_order[i]
	diff_test_sct_regressSmoke_sig$cluster[which(diff_test_sct_regressSmoke_sig$gene %in% 
		rownames(gene2colors)[which(gene2colors == cluster2colors[cluster_order[i]])])] <- letterVec[i]
	diff_test_sct_regressSmoke_sig$clusters_merged[which(diff_test_sct_regressSmoke_sig$gene %in% 
		rownames(gene2colors)[which(gene2colors == cluster2colors[cluster_order[i]])])] <- mergedVec[i]
}

#Order
diff_test_sct_regressSmoke_sig<-diff_test_sct_regressSmoke_sig[
 	with(diff_test_sct_regressSmoke_sig,order(clusters_merged,qval)),
]

##Export table 
write.table(diff_test_sct_regressSmoke_sig,row.names=F,sep="\t",quote=F,file="DEG_Tables/PseudotimeGenes_sct_regressSmoke_ordered_16_topDEGs.txt")

#Make a version with just IPA transcription factors
TFs<-read.table("../../Gene_lists/database_genes_IPA.txt",sep="\t",header=T,stringsAsFactors=F)
write.table(diff_test_sct_regressSmoke_sig[which(diff_test_sct_regressSmoke_sig$gene %in% TFs$Gene[which(TFs$Type == "transcription regulator")]),],
	row.names=F,sep="\t",quote=F,file="DEG_Tables/PseudotimeGenes_sct_regressSmoke_ordered_16_topDEGs_TFsOnly.txt")

#Do enrichments on the DEGs for each cluster
#Export table for Enrichr
diff_test_sct_regressSmoke_sig_forEnrichr<-diff_test_sct_regressSmoke_sig
colnames(diff_test_sct_regressSmoke_sig_forEnrichr)[6]<-"comparison"

#Write all to single txt file
write.table(diff_test_sct_regressSmoke_sig_forEnrichr,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,
	file="Enrichr/diff_test_sct_regressSmoke_ordered_16_topDEGs_forEnrichr.txt")

#Do enrichments
DEG_table<-read.table("Enrichr/diff_test_sct_regressSmoke_ordered_16_topDEGs_forEnrichr.txt",header=T,stringsAsFactors=F)
dataset<-"diff_test_sct_regressSmoke_ordered_16_topDEGs"
EnrichrAPI_location<-"/usr/local/bin/enrichrAPI.py"
doEnrichOneAtATime(DEG_table=DEG_table,dataset=dataset,EnrichrAPI_location=EnrichrAPI_location)










##########################################################Now re-make heatmap in the order specified, without clustering
#bsub -e test.err -o test.out -R "rusage[mem=64000]" "module load R/3.5.1 && env R_MAX_VSIZE=700Gb R CMD BATCH 6_Heatmap_pseudotimeDEGs_lumSec_16_reorder.r"

#Get annotations, adding in State
annotations1<-getBinnedValues(pData=pData(HSMM_sct),metadata_col_index=c(45))
annotations<-annotations1

pheatmap_16_ordered<-plot_pseudotime_heatmap(HSMM_sct[genes_subset_clusters,],
	cluster_rows=F,
	cores = detectCores(),
	show_rownames = F,
	add_annotation_col=annotations,
	return_heatmap = T)

#Print tree
pdf("Heatmap_PseudotimeGenes_sct_16_ordered.pdf")
pheatmap_16_ordered
dev.off()













########################################### Now change all the colors again for the final plot

#I     ####### Legend colors ########
### Here's where the states legend subclusters_4 colors are stored
pheatmap_16_ordered[[4]][[1]][[4]]$children[[6]]$gp$fill


#II     ######### Annotation bar colors ##########
##### Here's where the column header colors are stored, in order
pheatmap_16_ordered[[4]][[1]][[2]]$gp$fill





#########################################################Let's modify the colors first


########## Change the pseudotime cell annotation colors
######Subclusters_4
###First, need to change the annotation bar colors
#New colors
colorVec<-c("slategrey","turquoise","saddlebrown")
for(i in 1:length(pheatmap_16_ordered[[4]][[1]][[4]]$children[[6]]$gp$fill)){
	pheatmap_16_ordered[[4]][[1]][[2]]$gp$fill[,which(colnames(pheatmap_16_ordered[[4]][[1]][[2]]$gp$fill) == "subclusters_4")][
		which(pheatmap_16_ordered[[4]][[1]][[2]]$gp$fill[,which(colnames(pheatmap_16_ordered[[4]][[1]][[2]]$gp$fill) == "subclusters_4")] ==
		pheatmap_16_ordered[[4]][[1]][[4]]$children[[6]]$gp$fill[i])] <- colorVec[i]
}
###Then, change the legend
for(i in 1:length(pheatmap_16_ordered[[4]][[1]][[4]]$children[[6]]$gp$fill)){
	pheatmap_16_ordered[[4]][[1]][[4]]$children[[6]]$gp$fill[i] <- colorVec[i]
}




########### Now print the final heat map
#Print tree
pdf("Heatmap_PseudotimeGenes_sct_16_topDEGs_ordered.pdf")
pheatmap_16_ordered
dev.off()

#Print tree
png("Heatmap_PseudotimeGenes_sct_16_topDEGs_ordered.png",height=5,width=5,units="in",res=300)
pheatmap_16_ordered
dev.off()












############# Finally, get the positions of gene clusters in the heat map for use in affinity designer

DEG_table<-read.table("DEG_Tables/PseudotimeGenes_sct_regressSmoke_ordered_16_topDEGs.txt",header=T,stringsAsFactors=F)
bottomLocation = 424

getCoordinateVec(DEG_table=DEG_table,bottomLocation=424)
