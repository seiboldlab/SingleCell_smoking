library(Seurat, lib.loc="/Users/jacksonna/R-dev")
library(plotrix)
library(scales)
library(ggplot2)
library(cowplot)
library(openxlsx)
source("CommonFunctions.r")






##################################################### 1. Bring in 10X raw count matrices for each donor
######### These can be downloaded from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134174
T84<-read.table("T84_invivo_expression_matrix.txt",header=T,sep="\t")
T85<-read.table("T85_invivo_expression_matrix.txt",header=T,sep="\t")
T89<-read.table("T89_invivo_expression_matrix.txt",header=T,sep="\t")
T120<-read.table("T120_invivo_expression_matrix.txt",header=T,sep="\t")
T121<-read.table("T121_invivo_expression_matrix.txt",header=T,sep="\t")
T126<-read.table("T126_invivo_expression_matrix.txt",header=T,sep="\t")
T137<-read.table("T137_invivo_expression_matrix.txt",header=T,sep="\t")
T90<-read.table("T90_invivo_expression_matrix.txt",header=T,sep="\t")
T101<-read.table("T101_invivo_expression_matrix.txt",header=T,sep="\t")
T153<-read.table("T153_invivo_expression_matrix.txt",header=T,sep="\t")
T154<-read.table("T154_invivo_expression_matrix.txt",header=T,sep="\t")
T164<-read.table("T164_invivo_expression_matrix.txt",header=T,sep="\t")
T165<-read.table("T165_invivo_expression_matrix.txt",header=T,sep="\t")
T166<-read.table("T166_invivo_expression_matrix.txt",header=T,sep="\t")
T167<-read.table("T167_invivo_expression_matrix.txt",header=T,sep="\t")








##################################################### 2. Make summary stats table for each of the samples and decide what to keep
#### Make list of summary tables
sumTabList<-list()
namesVec<-c("T101","T120","T121","T126","T137","T153","T154","T164","T165","T166","T167","T84","T85","T89","T90")
for(i in 1:length(namesVec)){
	currData<-as.matrix(get(namesVec[i]))
	sumTabList[[length(sumTabList) + 1]]<-as.data.frame(matrix(NA,nrow=ncol(currData),ncol=4))
	names(sumTabList)[[length(sumTabList)]]<-namesVec[i]
	colnames(sumTabList[[length(sumTabList)]])<-c("ID","nUMI","nGene","propMT")
	sumTabList[[length(sumTabList)]]$ID<-colnames(currData)
	sumTabList[[length(sumTabList)]]$nUMI<-colSums(currData)
	sumTabList[[length(sumTabList)]]$nGene<-colSums(currData != 0)
	mtGenes<-grep("^MTAT|^MT-|^MTCO|^MTCY|^MTERF|^MTND|^MTRF|^MTRN|^MRPL|^MRPS",rownames(currData),value=T)
	sumTabList[[length(sumTabList)]]$propMT<-colSums(currData[mtGenes,]) / sumTabList[[length(sumTabList)]]$nUMI
	riboGenes<-grep("^RPL|^RPS",rownames(currData),value=T)
	sumTabList[[length(sumTabList)]]$propRibo<-colSums(currData[riboGenes,]) / sumTabList[[length(sumTabList)]]$nUMI
}

### Export these tables to xcel
write.xlsx(sumTabList, file = "10X_T15_summaryStats.xlsx")









########################### Histograms of summary stats (separate)

#nGenes
#ColorVector
colVec<-c("black","saddlebrown","red","yellow","blue","green","orange","purple","grey","tan","turquoise","lightblue","magenta","darkgreen","pink")
#Equal size bins
bins <- function(x){
	c(seq(min(x),max(x),by = 200),max(x))
}
pdf("10X_T15_nGene_hist_separate.pdf")
par(mfrow=c(3,3))
for(i in 1:length(namesVec)){	
	currData<-sumTabList[[i]]
	hist(currData$nGene,breaks=bins,las=1,main = names(sumTabList)[[i]],xlab="Number of genes",col=colVec[i],ylab="")
	abline(v=quantile(currData$nGene,0.99),col="red")
	abline(v=1500,col="red")
}	
dev.off()


#nUMI
colVec<-c("black","saddlebrown","red","yellow","blue","green","orange","purple","grey","tan","turquoise","lightblue","magenta","darkgreen","pink")
#Equal size bins
bins <- function(x){
	c(seq(min(x),max(x),by = 6000),max(x))
}
pdf("10X_T15_nUMI_hist_separate.pdf")
par(mfrow=c(3,3))
for(i in 1:length(namesVec)){	
	currData<-sumTabList[[i]]
	hist(currData$nUMI,breaks=bins,las=1,main = names(sumTabList)[[i]],xlab="Number of UMIs",col=colVec[i],ylab="")
	abline(v=quantile(currData$nUMI,0.99),col="red")
#	abline(v=quantile(currData$nUMI,0.01),col="red")
}	
dev.off()


#propMT
colVec<-c("black","saddlebrown","red","yellow","blue","green","orange","purple","grey","tan","turquoise","lightblue","magenta","darkgreen","pink")

pdf("10X_T15_propMT_hist_separate.pdf")
par(mfrow=c(3,3))
for(i in 1:length(namesVec)){	
	currData<-sumTabList[[i]]
	hist(currData$propMT,breaks=50,las=1,main = names(sumTabList)[[i]],xlab="Proportion MT genes",col=colVec[i],ylab="")
	abline(v=0.3,col="red")
}	
dev.off()










##################################################### 3. Filter out cells, merge datasets, and filter out genes

#Identify the samples that are outliers to be removed
removeList<-list()
for(i in 1:length(namesVec)){
removeList[[length(removeList) + 1]]<-sumTabList[[i]]$ID[which(sumTabList[[i]]$nGene > quantile(sumTabList[[i]]$nGene,0.99) | 
	sumTabList[[i]]$nGene < 1500 | sumTabList[[i]]$nUMI > quantile(sumTabList[[i]]$nUMI,0.99) | sumTabList[[i]]$propMT > 0.3)]
names(removeList)[[length(removeList)]]<-namesVec[i]
} 

#Original number of cells #68,370
NumCellsVec<-sapply(sumTabList,function(x)nrow(x))
#T101 T120 T121 T126 T137 T153 T154 T164 T165 T166 T167  T84  T85  T89  T90 
#5834 2871 2382 2702 2469 5565 3676 5650 4835 6167 5983 4224 7226 5282 3504 

#Number gone (29,661 total tossed)
NumGoneVec<-sapply(removeList,function(x)length(x)) 
#T101 T120 T121 T126 T137 T153 T154 T164 T165 T166 T167  T84  T85  T89  T90 
#2287  957  875  834 1042 2006 1831 2221 2444 3070 2887 2291 3345 2469 1133 

#Proportion gone
PropGoneVec<-NumGoneVec / NumCellsVec
#     T101      T120      T121      T126      T137      T153      T154      T164 
#0.3920123 0.3333333 0.3673384 0.3086603 0.4220332 0.3604672 0.4980958 0.3930973 
#     T165      T166      T167       T84       T85       T89       T90 
#0.5054809 0.4978109 0.4825338 0.5423769 0.4629117 0.4674366 0.3233447

#Number remaining (38,678 total left)
NumRemainVec<-NumCellsVec - NumGoneVec
#T101 T120 T121 T126 T137 T153 T154 T164 T165 T166 T167  T84  T85  T89  T90 
#3547 1914 1507 1868 1427 3559 1845 3429 2391 3097 3096 1933 3881 2813 2371  

#Proportion remaining
PropRemainVec<-NumRemainVec / NumCellsVec
#     T101      T120      T121      T126      T137      T153      T154      T164 
#0.6079877 0.6666667 0.6326616 0.6913397 0.5779668 0.6395328 0.5019042 0.6069027 
#     T165      T166      T167       T84       T85       T89       T90 
#0.4945191 0.5021891 0.5174662 0.4576231 0.5370883 0.5325634 0.6766553 


#Also, remove 2,461 cell identified as non-epithelial in the preliminary analysis
removedCells<-scan("Gene_lists/removedCells_scTransform_3000Features.txt",what="character")
for(i in 1:length(removeList)){
	currRemovedCells<-removedCells[grep(paste("_",i,"$",sep=""),removedCells)]
	currRemovedCells<-sapply(strsplit(currRemovedCells,"_"),function(x)x[1])
	removeList[[i]]<-unique(append(removeList[[i]],currRemovedCells))
} 

#Number gone (32,122 total tossed)
NumGoneVec<-sapply(removeList,function(x)length(x)) 
#T101 T120 T121 T126 T137 T153 T154 T164 T165 T166 T167  T84  T85  T89  T90 
#2539 1025  949  904 1118 2131 2026 2403 2717 3540 3116 2341 3451 2588 1274 

#Number remaining (36,248 total left)
NumRemainVec<-NumCellsVec - NumGoneVec
#T101 T120 T121 T126 T137 T153 T154 T164 T165 T166 T167  T84  T85  T89  T90 
#3295 1846 1433 1798 1351 3434 1650 3247 2118 2627 2867 1883 3775 2694 2230  


##### Filter out samples
T101f<-as.matrix(T101)[,-which(colnames(T101) %in% removeList$T101)]
T120f<-as.matrix(T120)[,-which(colnames(T120) %in% removeList$T120)]
T121f<-as.matrix(T121)[,-which(colnames(T121) %in% removeList$T121)]
T126f<-as.matrix(T126)[,-which(colnames(T126) %in% removeList$T126)]
T137f<-as.matrix(T137)[,-which(colnames(T137) %in% removeList$T137)]
T153f<-as.matrix(T153)[,-which(colnames(T153) %in% removeList$T153)]
T154f<-as.matrix(T154)[,-which(colnames(T154) %in% removeList$T154)]
T164f<-as.matrix(T164)[,-which(colnames(T164) %in% removeList$T164)]
T165f<-as.matrix(T165)[,-which(colnames(T165) %in% removeList$T165)]
T166f<-as.matrix(T166)[,-which(colnames(T166) %in% removeList$T166)]
T167f<-as.matrix(T167)[,-which(colnames(T167) %in% removeList$T167)]
T84f<-as.matrix(T84)[,-which(colnames(T84) %in% removeList$T84)]
T85f<-as.matrix(T85)[,-which(colnames(T85) %in% removeList$T85)]
T89f<-as.matrix(T89)[,-which(colnames(T89) %in% removeList$T89)]
T90f<-as.matrix(T90)[,-which(colnames(T90) %in% removeList$T90)]








############### Combine datasets and remove genes
#Merge the 8 samples into 1
T15f<-cbind(T101f,T120f,T121f,T126f,T137f,T153f,T154f,T164f,T165f,T166f,T167f,T84f,T85f,T89f,T90f)

#Add an identifier to barcodes, so there will be no duplicates for each donor
colnames(T15f)<-c(paste(colnames(T101f),"1",sep="_"),paste(colnames(T120f),"2",sep="_"),paste(colnames(T121f),"3",sep="_"),paste(colnames(T126f),"4",sep="_"),
	paste(colnames(T137f),"5",sep="_"),paste(colnames(T153f),"6",sep="_"),paste(colnames(T154f),"7",sep="_"),paste(colnames(T164f),"8",sep="_"),
	paste(colnames(T165f),"9",sep="_"),paste(colnames(T166f),"10",sep="_"),paste(colnames(T167f),"11",sep="_"),
	paste(colnames(T84f),"12",sep="_"),paste(colnames(T85f),"13",sep="_"),paste(colnames(T89f),"14",sep="_"),paste(colnames(T90f),"15",sep="_"))

#Filter out genes (toss 223 genes; 33315 genes remaining
T15f<-T15f[-c(grep("MTAT|MT-|MTCO|MTCY|MTERF|MTND|MTRF|MTRN|MRPL|MRPS|RPL|RPS",rownames(T15f))),]

#Initially, just toss genes with zero expression (toss 5,263; 28,052 remaining)
T15f<-T15f[-which(rowSums(T15f) == 0),]

#Make meta.data table
T15_metadata<-data.frame("donor"=c(rep("T101",ncol(T101f)),rep("T120",ncol(T120f)),rep("T121",ncol(T121f)),rep("T126",ncol(T126f)),
	rep("T137",ncol(T137f)),rep("T153",ncol(T153f)),rep("T154",ncol(T154f)),
	rep("T164",ncol(T164f)),rep("T165",ncol(T165f)),rep("T166",ncol(T166f)),rep("T167",ncol(T167f)),
	rep("T84",ncol(T84f)),rep("T85",ncol(T85f)),rep("T89",ncol(T89f)),rep("T90",ncol(T90f))),
	
	"smoke"=c(rep("heavy",ncol(T101f)),rep("heavy",ncol(T120f)),rep("light",ncol(T121f)),rep("never",ncol(T126f)),
	rep("never",ncol(T137f)),rep("never",ncol(T153f)),rep("heavy",ncol(T154f)),
	rep("never",ncol(T164f)),rep("never",ncol(T165f)),rep("never",ncol(T166f)),rep("heavy",ncol(T167f)),
	rep("light",ncol(T84f)),rep("heavy",ncol(T85f)),rep("never",ncol(T89f)),rep("heavy",ncol(T90f))),
	
	"pack_years"=c(rep("25",ncol(T101f)),rep("90",ncol(T120f)),rep("0.5",ncol(T121f)),rep("0",ncol(T126f)),
	rep("0",ncol(T137f)),rep("0",ncol(T153f)),rep("60",ncol(T154f)),
	rep("0",ncol(T164f)),rep("0",ncol(T165f)),rep("0",ncol(T166f)),rep("46",ncol(T167f)),
	rep("3",ncol(T84f)),rep("15",ncol(T85f)),rep("0",ncol(T89f)),rep("30",ncol(T90f))),
	
	"age"=c(rep("55",ncol(T101f)),rep("57",ncol(T120f)),rep("23",ncol(T121f)),rep("35",ncol(T126f)),
	rep("27",ncol(T137f)),rep("38",ncol(T153f)),rep("61",ncol(T154f)),
	rep("66",ncol(T164f)),rep("64",ncol(T165f)),rep("68",ncol(T166f)),rep("66",ncol(T167f)),
	rep("22",ncol(T84f)),rep("59",ncol(T85f)),rep("10",ncol(T89f)),rep("44",ncol(T90f))),
	
	"sex"=c(rep("M",ncol(T101f)),rep("F",ncol(T120f)),rep("F",ncol(T121f)),rep("F",ncol(T126f)),
	rep("M",ncol(T137f)),rep("M",ncol(T153f)),rep("M",ncol(T154f)),
	rep("M",ncol(T164f)),rep("M",ncol(T165f)),rep("F",ncol(T166f)),rep("F",ncol(T167f)),
	rep("M",ncol(T84f)),rep("M",ncol(T85f)),rep("F",ncol(T89f)),rep("M",ncol(T90f))),

	"nUMI"=c(sumTabList_f[[1]]$nUMI,sumTabList_f[[2]]$nUMI,sumTabList_f[[3]]$nUMI,sumTabList_f[[4]]$nUMI,
	sumTabList_f[[5]]$nUMI,sumTabList_f[[6]]$nUMI,sumTabList_f[[7]]$nUMI,sumTabList_f[[8]]$nUMI,
	sumTabList_f[[9]]$nUMI,sumTabList_f[[10]]$nUMI,sumTabList_f[[11]]$nUMI,
	sumTabList_f[[12]]$nUMI,sumTabList_f[[13]]$nUMI,sumTabList_f[[14]]$nUMI,sumTabList_f[[15]]$nUMI),
	
	"nGene"=c(sumTabList_f[[1]]$nGene,sumTabList_f[[2]]$nGene,sumTabList_f[[3]]$nGene,sumTabList_f[[4]]$nGene,
	sumTabList_f[[5]]$nGene,sumTabList_f[[6]]$nGene,sumTabList_f[[7]]$nGene,sumTabList_f[[8]]$nGene,
	sumTabList_f[[9]]$nGene,sumTabList_f[[10]]$nGene,sumTabList_f[[11]]$nGene,
	sumTabList_f[[12]]$nGene,sumTabList_f[[13]]$nGene,sumTabList_f[[14]]$nGene,sumTabList_f[[15]]$nGene),
	
	"propMT"=c(sumTabList_f[[1]]$propMT,sumTabList_f[[2]]$propMT,sumTabList_f[[3]]$propMT,sumTabList_f[[4]]$propMT,
	sumTabList_f[[5]]$propMT,sumTabList_f[[6]]$propMT,sumTabList_f[[7]]$propMT,sumTabList_f[[8]]$propMT,
	sumTabList_f[[9]]$propMT,sumTabList_f[[10]]$propMT,sumTabList_f[[11]]$propMT,
	sumTabList_f[[12]]$propMT,sumTabList_f[[13]]$propMT,sumTabList_f[[14]]$propMT,sumTabList_f[[15]]$propMT),
	
	"propRibo"=c(sumTabList_f[[1]]$propRibo,sumTabList_f[[2]]$propRibo,sumTabList_f[[3]]$propRibo,sumTabList_f[[4]]$propRibo,
	sumTabList_f[[5]]$propRibo,sumTabList_f[[6]]$propRibo,sumTabList_f[[7]]$propRibo,sumTabList_f[[8]]$propRibo,
	sumTabList_f[[9]]$propRibo,sumTabList_f[[10]]$propRibo,sumTabList_f[[11]]$propRibo,
	sumTabList_f[[12]]$propRibo,sumTabList_f[[13]]$propRibo,sumTabList_f[[14]]$propRibo,sumTabList_f[[15]]$propRibo),
	
	row.names=colnames(T15f))