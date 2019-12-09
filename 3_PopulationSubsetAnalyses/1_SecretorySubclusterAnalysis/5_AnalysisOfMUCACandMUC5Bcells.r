library(Seurat, lib.loc="/Users/jacksonna/R-dev")
library(gplots)
library(reshape)
library(ppcor)
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


#First, bring in the monocle-derived mature secretory cell group and add as metadata
matureSecCells<-scan("monocle/ANALYSIS_lumSec_reintegratedData/Mature_secretory_cells_monocle.txt",what="character")
treat.col<-data.frame("Mature_secretory_cells_monocle"=rep(FALSE,nrow(T15_int@meta.data)),row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col[which(rownames(treat.col) %in% matureSecCells),]<-TRUE
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="Mature_secretory_cells_monocle")

#Now subset data to include only those mature secretory cells
Idents(T15_int) <- T15_int@meta.data$Mature_secretory_cells_monocle
T15_sec<-subset(T15_int,idents=TRUE)






############################Make pie charts for mucin panels

proportionsList<-list()
for(i in 1:1){
	currData<-as.matrix(GetAssayData(T15_sec,assay="SCT"))[c("MUC5AC","MUC5B"),]
	proportions<-c()
	proportions<-append(proportions,length(which(currData["MUC5AC",] != 0 & currData["MUC5B",] == 0)) / ncol(currData))
	proportions<-append(proportions,length(which(currData["MUC5AC",] == 0 & currData["MUC5B",] != 0)) / ncol(currData))
	proportions<-append(proportions,length(which(currData["MUC5AC",] != 0 & currData["MUC5B",] != 0)) / ncol(currData))
	proportions<-append(proportions,length(which(currData["MUC5AC",] == 0 & currData["MUC5B",] == 0)) / ncol(currData))
	names(proportions)<-c("MUC5AC","MUC5B","Both","Neither")
	proportionsList[[length(proportionsList) + 1]]<-proportions
}

#Plot pie
#pdf("T15_sec_MUCs_piechart.pdf")
for(i in 1:length(proportionsList)){
	pct <- round(proportionsList[[i]],3)*100
	lbls <- paste(pct,"%",sep="") # ad % to labels 
	colors<-c("cornflowerblue","forestgreen","darkorange","black")
	pie(proportionsList[[i]],labels = names(proportionsList[[i]]), col=colors, main = names(proportionsList)[[i]])
}



















############################ Show MUC distribution panels on the UMAP
T15_sct<-as.matrix(GetAssayData(T15_sec,assay="SCT"))

#####First, encode mucin status into metadata
treat.col<-data.frame("MUCs"=rep("Both",nrow(T15_sec@meta.data),row.names=rownames(T15_sec@meta.data),stringsAsFactors=F))
treat.col[,1]<-as.character(treat.col[,1])
treat.col$MUCs[which(T15_sct["MUC5AC",] != 0 & T15_sct["MUC5B",] == 0)]<-"MUC5AC"
treat.col$MUCs[which(T15_sct["MUC5AC",] == 0 & T15_sct["MUC5B",] != 0)]<-"MUC5B"
treat.col$MUCs[which(T15_sct["MUC5AC",] == 0 & T15_sct["MUC5B",] == 0)]<-"Neither"
rownames(treat.col)<-rownames(T15_sec@meta.data)
treat.col$MUCs<-factor(treat.col$MUCs,levels=c("MUC5AC","MUC5B","Both","Neither"))
T15_sec<-AddMetaData(T15_sec,metadata=treat.col,col.name="MUCs")


#Now, read in the T15_int object and do the below
treat.col<-data.frame("MUCs"=rep("other",nrow(T15_int@meta.data)),row.names=rownames(T15_int@meta.data),stringsAsFactors=F)
treat.col[rownames(T15_sec@meta.data),]<-as.character(T15_sec@meta.data$MUCs)
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="MUCs")
T15_int@meta.data$MUCs<-as.character(T15_int@meta.data$MUCs)
T15_int@meta.data$MUCs[which(is.na(T15_int@meta.data$MUCs))]<-"Other"
T15_int@meta.data$MUCs<-factor(T15_int@meta.data$MUCs,levels=c("MUC5AC","MUC5B","Both","Neither","Other"))

#Plot
colorVec<-c("cornflowerblue","forestgreen","orange","black","gray70")
#pdf("T15_UMAP_MUCs_secretoryCellsOnly.pdf",height=5,width=7)
dev.new(height=5,width=6.5)
DimPlot(T15_int, reduction='umap', label=F, pt.size = 0.01, cols=colorVec,group.by="MUCs") + NoAxes()




############################ Repeat this for ALL cells
T15_sct<-as.matrix(GetAssayData(T15_int,assay="SCT"))

#####First, encode mucin status into metadata
treat.col<-data.frame("MUCs"=rep("Both",nrow(T15_int@meta.data),row.names=rownames(T15_int@meta.data),stringsAsFactors=F))
treat.col[,1]<-as.character(treat.col[,1])
treat.col$MUCs[which(T15_sct["MUC5AC",] != 0 & T15_sct["MUC5B",] == 0)]<-"MUC5AC"
treat.col$MUCs[which(T15_sct["MUC5AC",] == 0 & T15_sct["MUC5B",] != 0)]<-"MUC5B"
treat.col$MUCs[which(T15_sct["MUC5AC",] == 0 & T15_sct["MUC5B",] == 0)]<-"Neither"
rownames(treat.col)<-rownames(T15_int@meta.data)
treat.col$MUCs<-factor(treat.col$MUCs,levels=c("MUC5AC","MUC5B","Both","Neither"))
T15_int<-AddMetaData(T15_int,metadata=treat.col,col.name="MUCs")


#Plot
colorVec<-c("cornflowerblue","forestgreen","orange","black")
#pdf("T15_UMAP_MUCs_allCells.pdf",height=5,width=7)
dev.new(height=5,width=6.5)
DimPlot(T15_int, reduction='umap', label=F, pt.size = 0.01, cols=colorVec,group.by="MUCs") + NoAxes()

























#########################Also, make a bar graph, breaking out smoking status

#Define clusters to plot
T15_sct<-as.matrix(GetAssayData(T15_sec,assay="SCT"))

donorList_ns<-unique(T15_sec@meta.data$donor[which(T15_sec@meta.data$smoke_noT89 == "never")])
donorList_hs<-unique(T15_sec@meta.data$donor[which(T15_sec@meta.data$smoke_noT89 == "heavy")])
pvalList<-list()

#Donor colors
colorVec<-c("pink","tan","greenyellow","yellow","lightblue","grey","midnightblue","black","orangered4","darkgreen","purple4","slateblue")

#Do boxplot or do barplot
do.boxplot <- T
#Make plot
#pdf("T15_sec_MUC_smoking_barplots.pdf")
#pdf("T15_sec_MUC_smoking_boxplots.pdf")
dev.new(height=3,width=3)
for(i in 1:1){
	currData_ns<-list()
	currData_hs<-list()
	for(j in 1:length(donorList_ns)){
		currData_ns[[length(currData_ns) + 1]]<-T15_sct[c("MUC5AC","MUC5B"),which(T15_sec@meta.data$donor == donorList_ns[j])]
		currData_hs[[length(currData_hs) + 1]]<-T15_sct[c("MUC5AC","MUC5B"),which(T15_sec@meta.data$donor == donorList_hs[j])]
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
	
	#Now test the significance in differences between groups
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
	
	


















#################### Now get correlations of genes with MUC5AC and MUC5B in cells expressing both
#################### When doing so, partial out the smoking status (only looking at non- and heavy smoking cells)
T15_sct<-as.matrix(GetAssayData(T15_sec,assay="SCT"))

#####First, encode mucin status into metadata
treat.col<-data.frame("MUCs"=rep("Both",nrow(T15_sec@meta.data),row.names=rownames(T15_sec@meta.data),stringsAsFactors=F))
treat.col[,1]<-as.character(treat.col[,1])
treat.col$MUCs[which(T15_sct["MUC5AC",] != 0 & T15_sct["MUC5B",] == 0)]<-"MUC5AC"
treat.col$MUCs[which(T15_sct["MUC5AC",] == 0 & T15_sct["MUC5B",] != 0)]<-"MUC5B"
treat.col$MUCs[which(T15_sct["MUC5AC",] == 0 & T15_sct["MUC5B",] == 0)]<-"Neither"
rownames(treat.col)<-rownames(T15_sec@meta.data)
treat.col$MUCs<-factor(treat.col$MUCs,levels=c("MUC5AC","MUC5B","Both","Neither"))
T15_sec<-AddMetaData(T15_sec,metadata=treat.col,col.name="MUCs")

######### Now do correlations for MUC5AC and MUC5B, in turn
#MUC5AC in Both cells
#Isolate cells that 1) express both MUC5AC and MUC5B and 2) are either nonsmoking or heavy smoking
bothCells<-which(T15_sec@meta.data$MUCs == "Both" & (T15_sec@meta.data$smoke_noT89 == "never" | 
	T15_sec@meta.data$smoke_noT89 == "heavy"))
#Get expression matrix to use
tempSecMat<-T15_sct
#Get numeric vector of nonsmokers and heavy smokers identities
z<-as.numeric(factor(T15_sec@meta.data$smoke_noT89[bothCells],levels=c("never","heavy")))
#Get correlations of all genes with MUC5AC (except MUC5AC, because this throws a singularity error), controlling for smoke status
AC_cor_Both_cells_smoke<-apply(tempSecMat[which(rowSums(tempSecMat[,bothCells]) > 0),bothCells][
	-grep("MUC5AC",rownames(tempSecMat[which(rowSums(tempSecMat[,bothCells]) > 0),bothCells])),],
	1,function(x)pcor.test(tempSecMat["MUC5AC",bothCells],x,z,method="spearman"))
#Pull out the lists of correlations and pvalues and put into a dataframe
tempMat<-as.data.frame(matrix(ncol=2,nrow=0))
AC_cor_Both_cells_smoke<-t(sapply(AC_cor_Both_cells_smoke,function(x)rbind(tempMat,as.data.frame(matrix(c(x$estimate,x$p.value),nrow=1)))))
colnames(AC_cor_Both_cells_smoke)<-c("correlation","pvalue")
AC_cor_Both_cells_smoke<-as.data.frame(AC_cor_Both_cells_smoke)
AC_cor_Both_cells_smoke[,1]<-simplify2array(AC_cor_Both_cells_smoke[,1])
AC_cor_Both_cells_smoke[,2]<-simplify2array(AC_cor_Both_cells_smoke[,2])
#Add in MUC5AC as a gene, so we have them all
AC_cor_Both_cells_smoke<-rbind(AC_cor_Both_cells_smoke,data.frame("correlation"=1,"pvalue"=0,row.names="MUC5AC"))
#Order and rank
AC_cor_Both_cells_smoke<-AC_cor_Both_cells_smoke[order(AC_cor_Both_cells_smoke$correlation,decreasing=T),]
AC_cor_Both_cells_smoke<-data.frame("gene"=rownames(AC_cor_Both_cells_smoke),"cor"=AC_cor_Both_cells_smoke$correlation,
	"pvalue"=AC_cor_Both_cells_smoke$pvalue,"qvalue"=p.adjust(AC_cor_Both_cells_smoke$pvalue,method="fdr"),
	"rank"=seq(nrow(AC_cor_Both_cells_smoke)),row.names=seq(nrow(AC_cor_Both_cells_smoke)),stringsAsFactors=F)




#MUC5B in Both cells
#Isolate cells that 1) express both MUC5AC and MUC5B and 2) are either nonsmoking or heavy smoking
bothCells<-which(T15_sec@meta.data$MUCs == "Both" & (T15_sec@meta.data$smoke_noT89 == "never" | 
	T15_sec@meta.data$smoke_noT89 == "heavy"))
#Get expression matrix to use
tempSecMat<-T15_sct
#Get numeric vector of nonsmokers and heavy smokers identities
z<-as.numeric(factor(T15_sec@meta.data$smoke_noT89[bothCells],levels=c("never","heavy")))
#Get correlations of all genes with MUC5B (except MUC5B), controlling for smoke status
B_cor_Both_cells_smoke<-apply(tempSecMat[which(rowSums(tempSecMat[,bothCells]) > 0),bothCells][
	-grep("MUC5B",rownames(tempSecMat[which(rowSums(tempSecMat[,bothCells]) > 0),bothCells])),],
	1,function(x)pcor.test(tempSecMat["MUC5B",bothCells],x,z,method="spearman"))
#Pull out the lists of correlations and pvalues and put into a dataframe
tempMat<-as.data.frame(matrix(ncol=2,nrow=0))
B_cor_Both_cells_smoke<-t(sapply(B_cor_Both_cells_smoke,function(x)rbind(tempMat,as.data.frame(matrix(c(x$estimate,x$p.value),nrow=1)))))
colnames(B_cor_Both_cells_smoke)<-c("correlation","pvalue")
B_cor_Both_cells_smoke<-as.data.frame(B_cor_Both_cells_smoke)
B_cor_Both_cells_smoke[,1]<-simplify2array(B_cor_Both_cells_smoke[,1])
B_cor_Both_cells_smoke[,2]<-simplify2array(B_cor_Both_cells_smoke[,2])
#Add in MUC5B as a gene, so we have them all
B_cor_Both_cells_smoke<-rbind(B_cor_Both_cells_smoke,data.frame("correlation"=1,"pvalue"=0,row.names="MUC5B"))
#Order and rank
B_cor_Both_cells_smoke<-B_cor_Both_cells_smoke[order(B_cor_Both_cells_smoke$correlation,decreasing=T),]
B_cor_Both_cells_smoke<-data.frame("gene"=rownames(B_cor_Both_cells_smoke),"cor"=B_cor_Both_cells_smoke$correlation,
	"pvalue"=B_cor_Both_cells_smoke$pvalue,"qvalue"=p.adjust(B_cor_Both_cells_smoke$pvalue,method="fdr"),
	"rank"=seq(nrow(B_cor_Both_cells_smoke)),row.names=seq(nrow(B_cor_Both_cells_smoke)),stringsAsFactors=F)




#Now, assemble into lists of tables for export to excel file
#Write to excel file
datasetList<-list(AC_cor_Both_cells_smoke,B_cor_Both_cells_smoke)
write.xlsx(datasetList, file = "CorGenes_MUCs_rankings_T15_secretory_cells_smoke.xlsx")




##Now again, which top correlated genes with MUC5AC are also top correlated genes with MUC5B, and which are unique?
rownames(AC_cor_Both_cells_smoke)<-AC_cor_Both_cells_smoke$gene
rownames(B_cor_Both_cells_smoke)<-B_cor_Both_cells_smoke$gene

#Identify genes to plot 
AC_genes<-which(AC_cor_Both_cells_smoke$qval < 0.05 & AC_cor_Both_cells_smoke$cor > 0)
B_genes<-which(B_cor_Both_cells_smoke$qval < 0.05 & B_cor_Both_cells_smoke$cor > 0)

#Make table of shared genes
shared_AC<-AC_cor_Both_cells_smoke[AC_genes,][which(AC_cor_Both_cells_smoke[AC_genes,]$gene %in% B_cor_Both_cells_smoke[B_genes,]$gene),]
shared_B<-B_cor_Both_cells_smoke[B_genes,][which(B_cor_Both_cells_smoke[B_genes,]$gene %in% AC_cor_Both_cells_smoke[AC_genes,]$gene),]
shared_combined<-shared_AC
colnames(shared_combined)<-c("gene","cor_AC","p_AC","q_AC","rank_AC")
shared_combined<-data.frame(shared_combined,"cor_B"=NA,"p_B"=NA,"q_B"=NA,"rank_B"=NA)
shared_combined$cor_B<-shared_B[shared_combined$gene,]$cor
shared_combined$p_B<-shared_B[shared_combined$gene,]$pvalue
shared_combined$q_B<-shared_B[shared_combined$gene,]$qvalue
shared_combined$rank_B<-shared_B[shared_combined$gene,]$rank


#Make table of uniquely highly correlated genes with MUC5AC
unique_AC<-AC_cor_Both_cells_smoke[AC_genes,][-which(AC_cor_Both_cells_smoke[AC_genes,]$gene %in% B_cor_Both_cells_smoke[B_genes,]$gene),]
unique_AC<-unique_AC
colnames(unique_AC)<-c("gene","cor_AC","p_AC","q_AC","rank_AC")
unique_AC<-data.frame(unique_AC,"cor_B"=NA,"p_B"=NA,"q_B"=NA,"rank_B"=NA)
unique_AC$cor_B<-B_cor_Both_cells_smoke[unique_AC$gene,]$cor
unique_AC$p_B<-B_cor_Both_cells_smoke[unique_AC$gene,]$pvalue
unique_AC$q_B<-B_cor_Both_cells_smoke[unique_AC$gene,]$qvalue
unique_AC$rank_B<-B_cor_Both_cells_smoke[unique_AC$gene,]$rank


#Make table of uniquely highly correlated genes with MUC5B
unique_B<-B_cor_Both_cells_smoke[B_genes,][-which(B_cor_Both_cells_smoke[B_genes,]$gene %in% AC_cor_Both_cells_smoke[AC_genes,]$gene),]
unique_B<-unique_B
colnames(unique_B)<-c("gene","cor_B","p_B","q_B","rank_B")
unique_B<-data.frame(unique_B,"cor_AC"=NA,"p_AC"=NA,"q_AC"=NA,"rank_AC"=NA)
unique_B$cor_AC<-AC_cor_Both_cells_smoke[unique_B$gene,]$cor
unique_B$p_AC<-AC_cor_Both_cells_smoke[unique_B$gene,]$pvalue
unique_B$q_AC<-AC_cor_Both_cells_smoke[unique_B$gene,]$qvalue
unique_B$rank_AC<-AC_cor_Both_cells_smoke[unique_B$gene,]$rank


#Now export these to a excel file
datasetList<-list("shared_AC_and_B"=shared_combined,"unique_AC"=unique_AC,"unique_B"=unique_B)
write.xlsx(datasetList, file = "CorGenes_MUCs_BothCells_sharedAndUnique_smoke_all.xlsx")


#Also, export uniques again, but excluding genes with p values < 0.05 in the OTHER group
unique_AC_strong<-unique_AC[-which(unique_AC$p_B < 0.05),]
unique_B_strong<-unique_B[-which(unique_B$p_AC < 0.05),]
datasetList<-list("unique_AC"=unique_AC_strong,"unique_B"=unique_B_strong)
write.xlsx(datasetList, file = "CorGenes_MUCs_BothCells_stronglyUniqueOnly_smoke.xlsx")





#Make scatter plot for Figure 2
#Let's just plot the strongest genes
#For shared, plot only those genes with q < 0.05 and correlation > 0.15
#For unique, plot only those genes with q < 0.05 and correlation > 0.15 in one, but q > 0.05 in the other
#Exclude MUC5AC and MUC5B
cutoff<-0.15
genesToPlot<-shared_combined[which(shared_combined$cor_AC > cutoff & shared_combined$cor_B > cutoff),]
genesToPlot<-rbind(genesToPlot,unique_AC[which(unique_AC$cor_AC > cutoff),])
genesToPlot<-rbind(genesToPlot,unique_B[which(unique_B$cor_B > cutoff),])
genesToPlot<-cbind(genesToPlot,"colorVec"=c(rep("darkorange",length(which(shared_combined$cor_AC > cutoff & shared_combined$cor_B > cutoff))),
	rep("cornflowerblue",length(which(unique_AC$cor_AC > cutoff))),rep("forestgreen",length(which(unique_B$cor_B > cutoff)))))
genesToPlot<-genesToPlot[-grep("MUC5",rownames(genesToPlot)),]

#pdf("scatterplot_MUC_correlated_genes.pdf")
dev.new(height=5,width=5)
plot(genesToPlot$cor_B~genesToPlot$cor_A,las=1,ylab="MUC5B correlations",xlab="MUC5AC correlations",pch=16,
	col=as.character(genesToPlot$colorVec),bty="l")
#text(genesToPlot$cor_A,genesToPlot$cor_B+0.012,genesToPlot$gene,cex=0.2)
abline(h=0.2,lty=2)
abline(v=0.2,lty=2)


























################################# Now looking at whether the MUC5AC correlates go up with smoking and MUC5B go down
cutoff<-0.15
uniqueWithCutoff_AC<-c(unique_AC$gene[which(unique_AC$cor_AC > cutoff)])
uniqueWithCutoff_B<-c(unique_B$gene[which(unique_B$cor_B > cutoff)])
	
uniqueAndShared_AC<-c(unique_AC$gene[which(unique_AC$cor_AC > cutoff)],
	shared_combined$gene[which(shared_combined$cor_AC > cutoff & shared_combined$cor_B > cutoff)])
uniqueAndShared_B<-c(unique_B$gene[which(unique_B$cor_B > cutoff)],
	shared_combined$gene[which(shared_combined$cor_AC > cutoff & shared_combined$cor_B > cutoff)])

#Just take top 25
uniqueWithCutoff_AC<-uniqueWithCutoff_AC[1:25]
uniqueWithCutoff_B<-uniqueWithCutoff_B[1:25]

###############Get mean expression of genes
######### UNIQUE AND SHARED
#all
compareAC_B_smoke_all<-data.frame("exp"=c(exp(colMeans(log(as.matrix(GetAssayData(T15_sec,assay="SCT"))[uniqueWithCutoff_AC,which(T15_sec@meta.data$smoke_noT89 == "never")] + 1))),
	exp(colMeans(log(as.matrix(GetAssayData(T15_sec,assay="SCT"))[uniqueWithCutoff_AC,which(T15_sec@meta.data$smoke_noT89 == "heavy")] + 1))),
	exp(colMeans(log(as.matrix(GetAssayData(T15_sec,assay="SCT"))[uniqueWithCutoff_B,which(T15_sec@meta.data$smoke_noT89 == "never")] + 1))),
	exp(colMeans(log(as.matrix(GetAssayData(T15_sec,assay="SCT"))[uniqueWithCutoff_B,which(T15_sec@meta.data$smoke_noT89 == "heavy")] + 1)))),
	"category"=c(rep("AC_ns",length(which(T15_sec@meta.data$smoke_noT89 == "never"))),
	rep("AC_sh",length(which(T15_sec@meta.data$smoke_noT89 == "heavy"))),
	rep("B_ns",length(which(T15_sec@meta.data$smoke_noT89 == "never"))),
	rep("B_sh",length(which(T15_sec@meta.data$smoke_noT89 == "heavy")))),stringsAsFactors=F)
	


#Create boxplots
#pdf("Boxplots_MUC_correlates_all_smokingEffect.pdf")
dev.new(height=3.5,width=2)
par(bty="n")
boxplot(compareAC_B_smoke_all$exp~as.factor(compareAC_B_smoke_all$category),las=1,col=c("black","red"),
	ylab="Mean normalized epxression",xaxt="n",outcex=0.5,outpch=16,medcol="white")
mtext(paste("AC_p = ",formatC(wilcox.test(compareAC_B_smoke_all$exp[which(compareAC_B_smoke_all$category=="AC_sh")],
	compareAC_B_smoke_all$exp[which(compareAC_B_smoke_all$category=="AC_ns")],alternative="greater")$p.value,format="e",digits=2),sep=""),
	side=3,line=-10,at=3,adj=1,cex=0.7)
mtext(paste("B_p = ",formatC(wilcox.test(compareAC_B_smoke_all$exp[which(compareAC_B_smoke_all$category=="B_sh")],
	compareAC_B_smoke_all$exp[which(compareAC_B_smoke_all$category=="B_ns")],alternative="less")$p.value,format="e",digits=2),sep=""),
	side=3,line=-10.8,at=3,adj=1,cex=0.7)




#Make supplementary table without averaging
compareAC_B_smoke_table<-data.frame(rbind(t(as.matrix(GetAssayData(T15_sec,assay="SCT"))[uniqueWithCutoff_AC,which(T15_sec@meta.data$smoke_noT89 == "never")]),
	t(as.matrix(GetAssayData(T15_sec,assay="SCT"))[uniqueWithCutoff_AC,which(T15_sec@meta.data$smoke_noT89 == "heavy")]),
	t(as.matrix(GetAssayData(T15_sec,assay="SCT"))[uniqueWithCutoff_B,which(T15_sec@meta.data$smoke_noT89 == "never")]),
	t(as.matrix(GetAssayData(T15_sec,assay="SCT"))[uniqueWithCutoff_B,which(T15_sec@meta.data$smoke_noT89 == "heavy")])),
	"category"=c(rep("MUC5AC_nonsmokers",length(which(T15_sec@meta.data$smoke_noT89 == "never"))),
	rep("MUC5AC_smokers",length(which(T15_sec@meta.data$smoke_noT89 == "heavy"))),
	rep("MUC5B_nonsmokers",length(which(T15_sec@meta.data$smoke_noT89 == "never"))),
	rep("MUC5B_smokers",length(which(T15_sec@meta.data$smoke_noT89 == "heavy")))),stringsAsFactors=F)


























######################################## Subcluster 5M for the supplement

T15_sec<-ScaleData(T15_sec,assay="integrated")
T15_sec<-FindVariableFeatures(T15_sec, selection.method = "vst", nfeatures = 3000,assay="integrated")

#Run pca
T15_sec<- RunPCA(T15_sec, npcs = 30, verbose = TRUE)
ElbowPlot(T15_sec,ndims=30) #How many PCs

#Look at heat maps
pdf("DimHeatmaps.pdf")
DimHeatmap(T15_sec, reduction="pca", dims = 1:6, cells = 500, balanced = TRUE)
DimHeatmap(T15_sec, reduction="pca", dims = 7:12, cells = 500, balanced = TRUE)
DimHeatmap(T15_sec, reduction="pca", dims = 13:18, cells = 500, balanced = TRUE)
DimHeatmap(T15_sec, reduction="pca", dims = 19:24, cells = 500, balanced = TRUE)
DimHeatmap(T15_sec, reduction="pca", dims = 25:30, cells = 500, balanced = TRUE)
dev.off()




################################Clustering
#Input
nPCs<-30
colorVec<-c("yellow3","wheat","yellowgreen","tomato","cornflowerblue","green")
	
#Run umap across different values of hyperparameters
#n_neighbors - In general this parameter should often be in the range 5 to 50 (lower values optimize for local structure; default is 30)
#min_dist - Sensible values are in the range 0.001 to 0.5 (lower values optimize for local structure; default is 0.3).
T15_sec <- RunUMAP(T15_sec, reduction.use = "pca", dims=1:nPCs, n.neighbors = 10, min.dist = 0.5)
T15_sec <- FindNeighbors(object = T15_sec, dims = 1:nPCs)
T15_sec <- FindClusters(T15_sec, reduction.type="pca", resolution=0.1, algorithm=1)

#pdf("T15_sec_UMAP_subclusters.pdf",height=5,width=7)
dev.new(height=5,width=7)
DimPlot(T15_sec, reduction='umap', label=T, pt.size = 2, cols=colorVec)
#FeaturePlot(T151, c("nCount_RNA"), cols=c("gray", "blue"), pt.size=1, min.cutoff="q5", max.cutoff="q95")






############################ Overlay MUC distributions
T15_sct<-as.matrix(GetAssayData(T15_sec,assay="SCT"))

#####First, encode mucin status into metadata
treat.col<-data.frame("MUCs"=rep("Both",nrow(T15_sec@meta.data),row.names=rownames(T15_sec@meta.data),stringsAsFactors=F))
treat.col[,1]<-as.character(treat.col[,1])
treat.col$MUCs[which(T15_sct["MUC5AC",] != 0 & T15_sct["MUC5B",] == 0)]<-"MUC5AC"
treat.col$MUCs[which(T15_sct["MUC5AC",] == 0 & T15_sct["MUC5B",] != 0)]<-"MUC5B"
treat.col$MUCs[which(T15_sct["MUC5AC",] == 0 & T15_sct["MUC5B",] == 0)]<-"Neither"
rownames(treat.col)<-rownames(T15_sec@meta.data)
treat.col$MUCs<-factor(treat.col$MUCs,levels=c("MUC5AC","MUC5B","Both","Neither"))
T15_sec<-AddMetaData(T15_sec,metadata=treat.col,col.name="MUCs")


#Plot
colorVec<-c("cornflowerblue","forestgreen","orange","black")
#pdf("T15_UMAP_MUCs_c5M.pdf",height=5,width=7)
dev.new(height=3,width=4.8)
DimPlot(T15_sec, reduction='umap', label=F, pt.size = 1.5, cols=colorVec,group.by="MUCs") + NoAxes()


#Find markers
m0<-FindMarkers(T15_sec,ident.1 = c("0"), min.pct = 0.1, only.pos=T, logfc.threshold = 0.25, test.use = "wilcox",assay="SCT")
m1<-FindMarkers(T15_sec,ident.1 = c("1"), min.pct = 0.1, only.pos=T, logfc.threshold = 0.25, test.use = "wilcox",assay="SCT")



##########################Overlay Montoro goblet genes
#Do feature plots across
dev.new(height=3,width=4.5)
FeaturePlot(T15_sec, c("SCGB1A1"), cols=c("gray", "blue"), pt.size=1.5, min.cutoff="q5", max.cutoff="q95")

#Bring in Montoro goblet 1 and 2 subset markers, get mean expression, and plot feature plot
goblet1<-scan("Gene_lists/Montoro.etal_Goblet1.txt",what="character")
goblet2<-scan("Gene_lists/Montoro.etal_Goblet2.txt",what="character")

gobletMean<-CalculateMeanExpression(T15_sct,list(goblet1,goblet2),method="geometric")
T15_sec<-AddMetaData(T15_sec,metadata=gobletMean,col.name=c("goblet1","goblet2"))

#pdf("T15_UMAP_goblet1_c5M.pdf",height=5,width=7)
dev.new(height=3,width=4.5)
FeaturePlot(T15_sec, c("goblet1"), cols=c("gray", "blue"), pt.size=1.5, min.cutoff="q5", max.cutoff="q95")
#pdf("T15_UMAP_goblet2_c5M.pdf",height=5,width=7)
dev.new(height=3,width=4.5)
FeaturePlot(T15_sec, c("goblet2"), cols=c("gray", "blue"), pt.size=1.5, min.cutoff="q5", max.cutoff="q95")



############################## Now use the 85% quantile as a cutoff
treat.col<-data.frame("g1VSg2_HIGH_LOW"=rep(NA,nrow(T15_sec@meta.data)),row.names=rownames(T15_sec@meta.data))
treat.col[,1]<-as.character(treat.col[,1])
treat.col$g1VSg2_HIGH_LOW[which(T15_sec@meta.data$goblet1 > quantile(T15_sec@meta.data$goblet1,0.50) &
	T15_sec@meta.data$goblet2 <= quantile(T15_sec@meta.data$goblet2,0.50))]<-"High in Goblet-1 markers"
treat.col$g1VSg2_HIGH_LOW[which(T15_sec@meta.data$goblet2 > quantile(T15_sec@meta.data$goblet2,0.50) &
	T15_sec@meta.data$goblet1 <= quantile(T15_sec@meta.data$goblet1,0.50))]<-"High in Goblet-2 markers"
treat.col$g1VSg2_HIGH_LOW[which(T15_sec@meta.data$goblet1 > quantile(T15_sec@meta.data$goblet1,0.50) &
	T15_sec@meta.data$goblet2 > quantile(T15_sec@meta.data$goblet2,0.50))]<-"High in both"
treat.col$g1VSg2_HIGH_LOW[which(T15_sec@meta.data$goblet1 <= quantile(T15_sec@meta.data$goblet1,0.50) &
	T15_sec@meta.data$goblet2 <= quantile(T15_sec@meta.data$goblet2,0.50))]<-"Low in both"
T15_sec<-AddMetaData(T15_sec,metadata=treat.col,col.name="g1VSg2_HIGH_LOW")
T15_sec@meta.data$g1VSg2_HIGH_LOW<-factor(T15_sec@meta.data$g1VSg2_HIGH_LOW,levels=c("High in Goblet-1 markers","High in Goblet-2 markers","High in both","Low in both"))

#pdf('T15_UMAP_Goblet1_AND_Goblet2_c5M_featurePlot_overlay.pdf')
dev.new(height=3,width=6)
DimPlot(T15_sec,cols=c("blue","red","darkturquoise","grey"),group.by="g1VSg2_HIGH_LOW",pt.size=1.5) + NoAxes()




