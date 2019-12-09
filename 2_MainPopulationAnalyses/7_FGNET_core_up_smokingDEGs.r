####FGNET for smoking response genes
library(dplyr)
library(FGNet)
library(igraph)
library(openxlsx)
library(Seurat, lib.loc="/Users/jacksonna/R-dev")
source("../../CommonFunctions.r")


################# First bring in all the smoking DEGs
##### Bring in relevant smoking DEG tables scattered about
smokeDEGs_up<-read.table("../../DEG_Tables/Smoking_up_DEGs_coreUniqueDesignations.txt",header=T,stringsAsFactors=F)
smokeDEGs_down<-read.table("../../DEG_Tables/Smoking_down_DEGs_coreUniqueDesignations.txt",header=T,stringsAsFactors=F)

smokeDEGs_up<-smokeDEGs_up[which(smokeDEGs_up$main_cil == "Core" & !duplicated(smokeDEGs_up$gene)),]
smokeDEGs_down<-smokeDEGs_down[which(smokeDEGs_down$main_cil == "Core" & !duplicated(smokeDEGs_down$gene)),]










######################################################FGNET

#############Isolate all core DEGs
gtl_smoke_up_core<-unique(smokeDEGs_up$gene)
gtl_smoke_down_core<-unique(smokeDEGs_down$gene)
cat(gtl_smoke_up_core,sep="\n")
cat(gtl_smoke_down_core,sep="\n")


###########Gene-Term Linker
###### Note if I want to pull down later from downloaded files, just run fea_gtLinker_getResults with alreadyDownloaded=T

#7025154 (minSupport = 3)
jobID_up <- fea_gtLinker(geneList=gtl_smoke_up_core,annotations = c("GO_Biological_Process", 		
	"GO_Molecular_Function","GO_Cellular_Component","KEGG_Pathways"),minSupport=5)
jobID_up<-1346682
feaResults_gtLinker_up<-fea_gtLinker_getResults(jobID=jobID_up)
















#Get html reports
FGNet_report(feaResults_gtLinker_up)
FGNet_report(feaResults_gtLinker_down)

















######feaResults_gtLinker_up
#Get incidence matrix with gene assignments to different metagroups

############ Bring in the job id I want from above
feaResults_gtLinker_up<-fea_gtLinker_getResults(jobID="7025154",alreadyDownloaded=T)

feaResults <- feaResults_gtLinker_up
incidMat <- fea2incidMat(feaResults)

#Let's filter the dataset the same way we did the report above
selectedGroups <- rep(FALSE, nrow(feaResults$metagroups))
selectedGroups[c(1,2,5,6,11,12,13)]<-TRUE

tmpFea <- feaResults
tmpFea$metagroups<-cbind(tmpFea$metagroups,select=selectedGroups)

#Name the incidence matrix
incidMatSelection <- fea2incidMat(tmpFea, 
  filterAttribute="select", filterOperator="!=",filterThreshold="TRUE")

#Make the network (exports an iGraph object)
tmpFea_igraph<-functionalNetwork(incidMatSelection,weighted=T)
FGNet_report(tmpFea,filterAttribute="select",filterOperator= "!=",
	filterThreshold="TRUE")

#Now export the iGraph objects to gml (for use in Cytoscape)
#Note that the gtLinker network is stored in commonGtSets
write_graph(tmpFea_igraph$iGraph$commonGtSets,"FGNet_network_up.gml",format="gml")




####################### Make node file
#Now create a node input file for cytoscape documenting the metagroup colors to use,
#DE results, etc.
#Note that to get betweenness and other network metrics, I can calculate these in cytoscape
#using Tools <- NetworkAnalyzer <- Network Analysis <- Analyze Network

#Make new data.frame to export to cytoscape
nodeColors<-data.frame("nodeName"=rownames(incidMatSelection$metagroupsMatr),"nodeColors"=NA,
	"groupNum"=NA,"singleton"=TRUE,"LFC"=NA,"FDR"=NA,"unique_status"=NA)

#First, all multi-cluster genes are white
nodeColors$nodeColors[which(rowSums(incidMatSelection$metagroupsMatr) > 1)] <- "white"
nodeColors$groupNum[which(rowSums(incidMatSelection$metagroupsMatr) > 1)] <- "multi"
nodeColors$singleton[which(rowSums(incidMatSelection$metagroupsMatr) > 1)] <- FALSE

#Now fill in the rest of the clusters
colVec<-c("blue","lightpink","yellow2","violet","darkolivegreen","tan","black","grey","magenta","purple","lightblue","darkgoldenrod")
clustVec<-colnames(incidMatSelection$metagroupsMatr)
for(i in 1:ncol(incidMatSelection$metagroupsMatr)){
	nodeColors$nodeColors[which(incidMatSelection$metagroupsMatr[,i] == 1 & nodeColors$singleton)] <- colVec[i]
	nodeColors$groupNum[which(incidMatSelection$metagroupsMatr[,i] == 1 & nodeColors$singleton)] <- clustVec[i]
}

#Now fill in all the rest of the info
currData<-smokeDEGs_up[-c(8:9,11)]
rownames(currData)<-currData$gene
nodeColors$LFC<-currData[as.character(nodeColors$nodeName),]$avg_logFC
nodeColors$FDR<--log(currData[as.character(nodeColors$nodeName),]$p_adj_FDR)
nodeColors$unique_status<-currData[as.character(nodeColors$nodeName),]$main_cil

#Write to cytoscape table
write.table(nodeColors,sep="\t",quote=F,row.names=F,file="FGNet_nodeColors_up.txt")





############# Now modify the edge file (to denote metagroup connections)
#Bring in the edge table from cytoscape
edgeTable<-read.csv("FGNet_core_network_up_edgeTable.csv",stringsAsFactors=F)

#Assign each node its value (based on alphabetical order, 0-(N-1) and add this to the node table
nodeColors<-data.frame("nodeNum"=seq(0,(nrow(nodeColors) - 1)),nodeColors)

#Now make a new "metagroupLink" column with the following characteristics
	#1. Edges connecting two nodes uniquely belonging to the same metagroup, assign that metagroup number ("groupNum")
	#2. Edges connecting two nodes, one unique to a metagroup and one multi, assign the metagroup number for the unique one
	#3. Edges connecting two multi nodes, assign to multi metagroup
	#Note there's never a case where connections involve unique metagroups that are not the same
#First, pull out the node connections into a list
connectionList<-lapply(strsplit(edgeTable$name," \\(interacts with\\) "),function(x)c(x[1],x[2]))
#Now, for each edge, specify a metagroupLink value based on the above rules
metagroupLink<-c()
for(i in 1:length(connectionList)){
	#Number 1 rule
	if(nodeColors$groupNum[which(nodeColors$nodeNum == connectionList[[i]][1])] != "multi" &
		nodeColors$groupNum[which(nodeColors$nodeNum == connectionList[[i]][2])] != "multi" &
		nodeColors$groupNum[which(nodeColors$nodeNum == connectionList[[i]][1])] == 
		nodeColors$groupNum[which(nodeColors$nodeNum == connectionList[[i]][2])]){
		metagroupLink<-append(metagroupLink,nodeColors$groupNum[which(nodeColors$nodeNum == connectionList[[i]][1])])
	}
	#Number 2 rule a
	if(nodeColors$groupNum[which(nodeColors$nodeNum == connectionList[[i]][1])] != "multi" &
		nodeColors$groupNum[which(nodeColors$nodeNum == connectionList[[i]][2])] == "multi"){
		metagroupLink<-append(metagroupLink,nodeColors$groupNum[which(nodeColors$nodeNum == connectionList[[i]][1])])
	}
	#Number 2 rule b
	if(nodeColors$groupNum[which(nodeColors$nodeNum == connectionList[[i]][1])] == "multi" &
		nodeColors$groupNum[which(nodeColors$nodeNum == connectionList[[i]][2])] != "multi"){
		metagroupLink<-append(metagroupLink,nodeColors$groupNum[which(nodeColors$nodeNum == connectionList[[i]][2])])
	}
	#Number 3 rule
	if(nodeColors$groupNum[which(nodeColors$nodeNum == connectionList[[i]][1])] == "multi" &
		nodeColors$groupNum[which(nodeColors$nodeNum == connectionList[[i]][2])] == "multi"){
		metagroupLink<-append(metagroupLink,"multi")
	}
}

#Now add the new column to the edge dataframe and export
#Then import the metagroupLink column to the existing cytoscape edge table, using "name" as the key
edgeTable<-data.frame(edgeTable[,1:4],"metagroupLink"=metagroupLink,edgeTable[,5:ncol(edgeTable)],stringsAsFactors=F)
write.table(edgeTable,row.names=F,quote=F,sep="\t",file="FGNet_core_network_up_edgeTable_withMetagroupLink.txt")







