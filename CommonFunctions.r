#Function for adding FDR and concatenating DEGs (list version)
OrganizeDF_list<-function(datasetList=NULL,comparisonVec=NULL,pctIndexList=NULL,counts=NULL,clusterAssignments=NULL){
	masterDataset<-data.frame()
	for(i in 1:length(datasetList)){
		geneNames<-rownames(get("datasetList")[[i]])
		markersTemp<-cbind("avg_logFC"=get("datasetList")[[i]]$avg_logFC,
			"p_val"=get("datasetList")[[i]]$p_val,"p_adj_FDR"=p.adjust(get("datasetList")[[i]]$p_val,method="fdr"),
			"pct.1"=get("datasetList")[[i]]$pct.1,"pct.2"=get("datasetList")[[i]]$pct.2)
		markersTemp<-as.data.frame(markersTemp)
		markersTemp<-cbind("gene"=geneNames,markersTemp)
		markersTemp<-cbind(markersTemp,"comparison"=comparisonVec[i])
		markersTemp$gene<-as.character(markersTemp$gene)
		#If you want to add in extra columns for raw count pct values (as the above now uses scTransform counts)
		if(!is.null(pctIndexList)){
			markersTemp<-cbind(markersTemp[1:6],"pct.1_counts"=NA,"pct.2_counts"=NA,"comparison"=markersTemp[7])
			markersTemp$pct.1_counts<-apply(counts[markersTemp$gene,which(clusterAssignments %in% unlist(pctIndexList[[i]][1]))],
				1,function(x)sum(x != 0)/length(x))
			markersTemp$pct.2_counts<-apply(counts[markersTemp$gene,which(clusterAssignments %in% unlist(pctIndexList[[i]][2]))],
				1,function(x)sum(x != 0)/length(x))
		}	
		masterDataset<-rbind(masterDataset,markersTemp)	
	}
	return(masterDataset)
}













##### Write DEG results to Excel
writeDEGsToExcel<-function(DEG_table,savedFile){
	wb<-createWorkbook()
	for(i in 1:length(unique(DEG_table$comparison))){
	 	addWorksheet(wb,as.character(unique(DEG_table$comparison)[i]))
	 	writeData(wb,i,DEG_table[which(DEG_table$comparison == unique(DEG_table$comparison)[i]),],keepNA=T)
	 	setColWidths(wb, sheet = i, cols = 1:ncol(DEG_table), widths = "auto")
	 	freezePane(wb,i,firstRow=T)
	 	saveWorkbook(wb,savedFile,overwrite = TRUE) 
	 } 	
}










##### Write list of DEG tables to Excel
writeListToExcel<-function(DEG_table_list,sheetNames,savedFile){
	wb<-createWorkbook()
	for(i in 1:length(DEG_table_list)){
	 	addWorksheet(wb,sheetNames[i])
	 	writeData(wb,i,DEG_table_list[[i]],keepNA=T)
	 	setColWidths(wb, sheet = i, cols = 1:ncol(DEG_table_list[[i]]), widths = "auto")
	 	freezePane(wb,i,firstRow=T)
	 	saveWorkbook(wb,savedFile,overwrite = TRUE) 
	 } 	
}
















#####Reads in excel tabs
#If returnSingleTable = T, adds a column with the tab info
#If returnSingleTable = F, tabs are in a list
readExcelTabs<-function(xlsxFile,returnSingleTable=F){
sheetNames<-openxlsx::getSheetNames(xlsxFile)
sheetList<-lapply(sheetNames,openxlsx::read.xlsx,xlsxFile=xlsxFile)
names(sheetList)<-sheetNames
  if(returnSingleTable){
  	sheetDF<-data.frame()
  	for(sn in sheetNames){
  		if(nrow(sheetList[[sn]]) > 0){
  			sheetDF<-rbind(sheetDF,cbind(sheetList[[sn]],"comparison"= sn))
  		}
  	}
  	return(sheetDF)
  }else{
  	return(sheetList)
  }
}




















##########ENRICHR FUNCTIONS
#First, export gene list where we take top (avg_logFC) x.genes number of genes for each comparison
	#Function for generating python input
runEnrichrAPIFromR<-function(ifile="Enrichr.geneList.txt",ofile="Enrichr.output.txt",
	libraries="Enrichr.libraries.txt",minOverlap=5,minAdjPval=0.05,sortParam="adjPval",
	EnrichrAPI_location="/Seibold/home/jacksonna/enrichrAPI.py",pythonVersion="python2.7"){
	
	command<-paste(pythonVersion,EnrichrAPI_location,"--ifile",ifile,"--ofile",ofile,"--libraries",libraries,
		"--minOverlap",minOverlap,"--minAdjPval",minAdjPval,"--summarize","--sort",sortParam,sep=" ")
	system(command,intern=T)
}


#This function takes a DEGs list from seurat and formats it for EnrichrAPI and then runs Enrichr
#If ngenes is a number, that many of genes are taken from each cluster/module. If ngenes is NULL, then all genes in the table are used	
#cluster: /Seibold/home/jacksonna/enrichrAPI.py (must have loaded python/2.7.12 - might need to be on the fat node to load python 2 and R/3.5.1 simultaneously)
#laptop: /usr/local/bin/enrichrAPI.py
DEG_enrich_now<-function(DEG_table,dataset,ngenes=NULL,Enrichr_dir="Enrichr",
	enrichr.libraries=c("GO_Cellular_Component_2018","GO_Biological_Process_2018","GO_Molecular_Function_2018",
	"KEGG_2019_Human","Reactome_2016"),minOverlap=3,minAdjPval=0.2,EnrichrAPI_location="/Seibold/home/jacksonna/enrichrAPI.py",
	pythonVersion="python2.7"){

	if(!dir.exists(Enrichr_dir)){system(paste("mkdir ",Enrichr_dir,sep=""))}
	DEG_table<-read.table(DEG_table,header=T,stringsAsFactors=F)
	if(!is.null(ngenes)){
		d <- data.table(DEG_table, key="comparison")
		DEG_table_top<-d[,head(.SD, ngenes), by=comparison]
		DEG_table_top<-as.data.frame(DEG_table_top)[c(2,1)]
		write.table(DEG_table_top,sep="\t",quote=FALSE,col.names=c("gene","module"),row.names=FALSE,
			file=paste(Enrichr_dir,"/",dataset,"_Enrichr.geneList.txt",sep=""))
	}else{
		write.table(DEG_table[,c(1,grep("comparison",colnames(DEG_table)))],sep="\t",quote=FALSE,col.names=c("gene","module"),row.names=FALSE,
			file=paste(Enrichr_dir,"/",dataset,"_Enrichr.geneList.txt",sep=""))
	}

	#Print enrichr libraries file
	cat(paste(enrichr.libraries,collapse="\n"),file=paste(Enrichr_dir,"/Enrichr.libraries.txt",sep=""))

	#Run pythonAPI
	ifile=paste(Enrichr_dir,"/",dataset,"_Enrichr.geneList.txt",sep="")
	ofile=paste(Enrichr_dir,"/",dataset,"_Enrichr.output.xlsx",sep="")
	libraries=paste(Enrichr_dir,"/","Enrichr.libraries.txt",sep="")
	runEnrichrAPIFromR(ifile=ifile,ofile=ofile,libraries=libraries,minOverlap=minOverlap,minAdjPval=minAdjPval,
		EnrichrAPI_location=EnrichrAPI_location)
}



#Run Enrichr
#Submit tables one at a time to avoid time out errors, then stich the results together later
doEnrichOneAtATime<-function(DEG_table,dataset,ngenes=NULL,Enrichr_dir="Enrichr",
	enrichr.libraries=c("GO_Cellular_Component_2018","GO_Biological_Process_2018","GO_Molecular_Function_2018",
	"KEGG_2019_Human","Reactome_2016"),minOverlap=3,minAdjPval=0.2,EnrichrAPI_location="/Seibold/home/jacksonna/enrichrAPI.py",
	pythonVersion="python2.7"){

	for(i in 1:length(unique(DEG_table$comparison))){
		currData<-DEG_table[which(DEG_table$comparison == unique(DEG_table$comparison)[i]),]
		DEG_table_curr<-paste("Enrichr/",dataset,"_forEnrichr_",i,".txt",sep="")
		dataset_curr<-paste(dataset,"_",i,sep="")
		write.table(currData,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,file=DEG_table_curr)
		DEG_enrich_now(DEG_table=DEG_table_curr,dataset=dataset_curr,ngenes=ngenes,Enrichr_dir=Enrichr_dir,EnrichrAPI_location=EnrichrAPI_location,
			enrichr.libraries=enrichr.libraries,minOverlap=minOverlap,minAdjPval=minAdjPval,pythonVersion=pythonVersion)
	}
	
	#Now read in all the outputs into a list and export them as a single file
	wb<-createWorkbook()
	for(i in 1:length(unique(DEG_table$comparison))){
		currPath<-paste("Enrichr/",dataset,"_",i,"_Enrichr.output.xlsx",sep="")
		a <- loadWorkbook(currPath)
		sheetNames <- names(a)
		#Bring in all sheets minus the summary
		for(j in 1:(length(sheetNames) - 1)){
	 		currData<-as.data.frame(readWorkbook(currPath,sheet = j))
	 	}
	 	#Now add the imported data to the workbook
	 	addWorksheet(wb,as.character(unique(DEG_table$comparison)[i]))
	 	writeData(wb, i, currData)
	 	widths<-c(25,80,8,8,8,12,15,100)
	 	setColWidths(wb, sheet = i, cols = 1:ncol(currData), widths = widths) #or use auto
	 	freezePane(wb,i,firstRow=T)
	 	saveWorkbook(wb, paste("Enrichr/",dataset,"_Enrichr.output.xlsx",sep=""), overwrite = TRUE)  	
	}
	
	#Remove all the single tables
	for(i in 1:length(unique(DEG_table$comparison))){
		file.remove(paste("Enrichr/",dataset,"_forEnrichr_",i,".txt",sep=""))
		file.remove(paste("Enrichr/",dataset,"_",i,"_Enrichr.output.xlsx",sep=""))
		file.remove(paste("Enrichr/",dataset,"_",i,"_Enrichr.geneList.txt",sep=""))
	}
}

















#Calculate mean expression across overlapping genes and place into meta.data
#Function takes the @data table and a named list of genes
#(method = arithmetic, geometric, or median)
CalculateMeanExpression<-function(dataset,genesList,method="arithmetic"){
	mean_exp_df<-data.frame(matrix(nrow=ncol(dataset),ncol=0))
	for(i in 1:length(genesList)){
		currData<-dataset[intersect(genesList[[i]],rownames(dataset)),]
		currData<-currData[which(rowSums(currData) > 0),]
		if(method == "geometric"){
			mean_exp_df<-cbind(mean_exp_df,data.frame(exp(colMeans(log(currData + 1)))))
			#mean_exp_df<-cbind(mean_exp_df,data.frame(geometric.mean(currData + 0.01))) #requires psych
		}
		if(method == "arithmetic"){
			mean_exp_df<-cbind(mean_exp_df,data.frame(colMeans(currData)))
		}
		if(method == "median"){
			mean_exp_df<-cbind(mean_exp_df,data.frame(robustbase::colMedians(currData)))
		}
	}
	colnames(mean_exp_df)<-names(genesList)
	return(mean_exp_df)
}













# Convert mouse to human gene names
convertMouseGeneList <- function(x){
	require("biomaRt")
	human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
 
	genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x, 
		mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
		
	humanx <- unique(genesV2[, 2])
	return(humanx)
}


















#Get conserved markers
##The approach here is, for each target cluster, compare to all the other clusters one at a time.
#From these, chose exemplar p-values, LFCs, and pct.2's based on the leaset differentiated comparisons.
#Output table
#Can set filters for genes to consider (based on min.pct and logfc.threshold)
#Only considers genes that show up for all comparisons
#Can either not downsample (NULL), downsample based on a number of cells (integer) or based on the median number of cells ("median")
GetConservedMarkers<-function(dataset,cluster.meta.data,min.pct=0.25,logfc.threshold=0.1,test.use="wilcox",
	latent.vars = NULL, downsample=NULL,print.excel=T){
	require("Seurat")
	if(print.excel){
		require("openxlsx")
	}
	DEGsCons<-data.frame()
	Idents(dataset)<-cluster.meta.data
	identVec<-unique(cluster.meta.data)
	
	#Set the maximum number of cells per identity class
	if(is.null(downsample)){max.cells.per.ident = Inf}
	if(is.numeric(downsample)){max.cells.per.ident = downsample}
	if(downsample == "median"){
		identSize<-c()
		for(k in 1:length(identVec)){
			identSize<-append(identSize,length(cluster.meta.data[which(cluster.meta.data == identVec[k])]))
		}
		max.cells.per.ident<-median(identSize)
	}
	
	#Run FindMarkers and save results
	for(i in 1:length(identVec)){
		masterMarker<-data.frame()
		for(j in 1:length(identVec[-i])){
			#Get current markers
			if(is.null(latent.vars)){
				markersCurr<-FindMarkers(dataset, ident.1 = identVec[i], ident.2 = identVec[-i][j], min.pct = min.pct, only.pos = T,
					logfc.threshold = logfc.threshold, test.use = "wilcox", max.cells.per.ident=max.cells.per.ident)
			}else{
				markersCurr<-FindMarkers(dataset, ident.1 = identVec[i], ident.2 = identVec[-i][j], min.pct = min.pct, only.pos = T,
					logfc.threshold = logfc.threshold, test.use = "wilcox", latent.vars = latent.vars, max.cells.per.ident=max.cells.per.ident)
			}
			#Calculate FDR
			geneNames<-rownames(markersCurr)
			markersCurr<-as.data.frame(cbind("avg_logFC"=markersCurr$avg_logFC,"p_val"=markersCurr$p_val,"p_adj_FDR"=
				p.adjust(markersCurr$p_val,method="fdr"),"pct.1"=markersCurr$pct.1,"pct.2"=markersCurr$pct.2))
			markersCurr<-cbind("gene"=geneNames,markersCurr,"comparison"=paste(identVec[i],identVec[-i][j],sep="_"),
				stringsAsFactors=F)
			#Add to master data.frame
			masterMarker<-rbind(masterMarker,markersCurr)
		}
		
		#Print master table
		if(!dir.exists("PairwiseDF")){dir.create("PairwiseDF")}
		if(!dir.exists("PairwiseDF/DEGs_txt")){dir.create("PairwiseDF/DEGs_txt")}
		write.table(masterMarker,sep="\t",quote=FALSE,row.names=FALSE,file=paste("PairwiseDF/DEGs_txt/DEGs_",identVec[i],".txt",sep=""))

		#Write to excel file
		if(print.excel){
			if(!dir.exists("PairwiseDF/DEGs_excel")){dir.create("PairwiseDF/DEGs_excel")}
			datasetList<-list()
			for(m in 1:length(unique(masterMarker$comparison))){
				datasetList[[length(datasetList) + 1]]<-masterMarker[which(masterMarker$comparison == unique(masterMarker$comparison)[m]),]
				names(datasetList)[length(datasetList)]<-as.character(unique(masterMarker$comparison)[m])
			}
			write.xlsx(datasetList,file=paste("PairwiseDF/DEGs_excel/DEGs_",identVec[i],".xlsx",sep=""))
		}
		
		#Now add to master table with all unique genes across comparisons for the current target cluster,using
		#the highest p-value, the lowest LFC, and largest pct.2 value
		#pct.2 mean gives the mean value of pct.2 rather than the largest value.
		#pct_ratio is pct.2 / pct.1, which gives the the difference in proportional expression between the two clusters (with lower values meaning more
		#discriminatory for the target cluster
		for(k in 1:length(unique(masterMarker$gene))){
			currGeneTab<-masterMarker[which(masterMarker$gene == unique(masterMarker$gene)[k]),]
			if(sum(!is.na(currGeneTab$p_val)) == length(unique(masterMarker$comparison))){
				DEGsCons<-rbind(DEGsCons,data.frame("gene"=currGeneTab$gene[1],"avg_logFC"=min(currGeneTab$avg_logFC),
					"p_val"=max(currGeneTab$p_val),"p_adj_FDR"=max(currGeneTab$p_adj_FDR),"pct.1"=mean(currGeneTab$pct.1),
					"pct.2"=max(currGeneTab$pct.2),"pct.2_mean"=mean(currGeneTab$pct.2),
					"pct_ratio"=max(currGeneTab$pct.2) / mean(currGeneTab$pct.1),"pct_ratio_mean"=mean(currGeneTab$pct.2) / mean(currGeneTab$pct.1),
					"comparison"=identVec[i],
					"comparison_logFC"=strsplit(currGeneTab$comparison[which(currGeneTab$avg_logFC == min(currGeneTab$avg_logFC))],"_")[[1]][2],
					"comparison_FDR"=strsplit(currGeneTab$comparison[which(currGeneTab$p_adj_FDR == max(currGeneTab$p_adj_FDR))],"_")[[1]][2],
					"comparison_pct.2"=strsplit(currGeneTab$comparison[which(currGeneTab$pct.2 == max(currGeneTab$pct.2))],"_")[[1]][2]))
			}
		}
	}
	
	#Re-calculate FDR based on the number of tests being considered
	DEGsCons$p_adj_FDR<-p.adjust(DEGsCons$p_val,method="fdr")
	
	#Sort table by comparison, the FDR	
	DEGsCons<-DEGsCons[
  		with(DEGsCons,order(comparison,p_adj_FDR)),
	]
	
	#Print conserved DEGs table
	write.table(DEGsCons,sep="\t",quote=FALSE,row.names=FALSE,file="PairwiseDF/Conserved_DEGs.txt")

	#Write to excel file
	if(print.excel){
		datasetList<-list()
		for(m in 1:length(unique(DEGsCons$comparison))){
			datasetList[[length(datasetList) + 1]]<-DEGsCons[which(DEGsCons$comparison == unique(DEGsCons$comparison)[m]),]
			names(datasetList)[length(datasetList)]<-as.character(unique(DEGsCons$comparison)[m])
		}
		write.xlsx(datasetList,file="PairwiseDF/Conserved_DEGs.xlsx")
	}
	return(DEGsCons)
}














#A sister function that takes a DEGsCons table and culls it for the most discriminating markers based on specified criteria 
#p_adj_FDR, avg_logFC, pct.1, and pct_ratio give hard cutoff values. These can be given as numeric values or as quantiles in the form of "q25", 
	#where the 25% quantile will be set as the cutoff
#Options for rank.by: "FDR", "pct_ratio", "FDR_pctRatio_combined"
#Note that single values for all clusters or vectors filled with cluster-specific values can be specified for minGene, maxGene, and rank.by
GetCulledConservedMarkersList<-function(DEGsCons,p_adj_FDR=1e-5,avg_logFC=0.25,pct.1=0.2,pct_ratio="q25",minGene=20,maxGene=20,rank.by="FDR_pctRatio_combined"){

	#Create list of genes
	genesList<-list()
	for(i in 1:length(unique(DEGsCons$comparison))){
		currDataOrig<-DEGsCons[which(DEGsCons$comparison == unique(DEGsCons$comparison)[i]),]
		
		#Get parameters (if quantiles)
		if(is.character(p_adj_FDR)){
			p_adj_FDR<-quantile(currDataOrig$p_adj_FDR,as.numeric(strsplit(p_adj_FDR,"q")[[1]][2]) * 0.01)
		}
		if(is.character(avg_logFC)){
			avg_logFC<-quantile(currDataOrig$avg_logFC,as.numeric(strsplit(avg_logFC,"q")[[1]][2]) * 0.01)
		}
		if(is.character(pct.1)){
			pct.1<-quantile(currDataOrig$pct.1,as.numeric(strsplit(pct.1,"q")[[1]][2]) * 0.01)
		}
		if(is.character(pct_ratio)){
			pct_ratio<-quantile(currDataOrig$pct_ratio,as.numeric(strsplit(pct_ratio,"q")[[1]][2]) * 0.01)
		}	
		
		#Determine whether minGene, maxGene, or rank.by are single values or vectors of values
		if(length(minGene)==length(unique(DEGsCons$comparison))){
			minGeneCurr<-minGene[i]
		}else{
			minGeneCurr<-minGene
		}
		
		if(length(maxGene)==length(unique(DEGsCons$comparison))){
			maxGeneCurr<-maxGene[i]
		}else{
			maxGeneCurr<-maxGene
		}
		
		if(length(rank.by)==length(unique(DEGsCons$comparison))){
			rank.byCurr<-rank.by[i]
		}else{
			rank.byCurr<-rank.by
		}
	
		#Sort table by pct_ratio	
		if(rank.byCurr == "pct_ratio"){
			currDataOrig<-currDataOrig[
 				with(currDataOrig,order(pct_ratio,p_adj_FDR)),
			]
		}
		
		#Sort table by FDR	
		if(rank.byCurr == "FDR"){
			currDataOrig<-currDataOrig[
 				with(currDataOrig,order(p_adj_FDR,pct_ratio)),
			]
		}

	   #Get combined rank for pct_ratio and FDR and sort by this
	   if(rank.byCurr == "FDR_pctRatio_combined"){
 	  		currDataOrig<-cbind(currDataOrig,"rank_pctFDR"=apply(data.frame(rank(currDataOrig$pct_ratio),rank(currDataOrig$p_adj_FDR)),1,mean))
   			currDataOrig<-currDataOrig[
 				with(currDataOrig,order(rank_pctFDR,p_adj_FDR)),
			]
		}
	
		#Apply filter to gene set
		currData<-currDataOrig[which(currDataOrig$p_adj_FDR < p_adj_FDR & 
			currDataOrig$avg_logFC > avg_logFC & 
			currDataOrig$pct.1 > pct.1 &
			currDataOrig$pct_ratio < pct_ratio),]
	
		#If we don't get the minimum number of genes, try adding in genes representing the best combo of FDR and pct_ratio
		if(nrow(currData) < minGeneCurr){
			currData<-rbind(currData,currDataOrig[which(currDataOrig$p_adj_FDR < median(currDataOrig$p_adj_FDR) & 
				currDataOrig$avg_logFC > avg_logFC & 
				currDataOrig$pct.1 > pct.1 &
				currDataOrig$pct_ratio < median(currDataOrig$pct_ratio)),])
			currData<-currData[!duplicated(currData),]
		}	
		
		#If we still don't get the minimum number of genes, just add in genes, best FDR first, until we have enough or until out of genes, whatever comes first
		if(nrow(currData) < minGeneCurr){
			currData<-rbind(currData,currDataOrig[-which(currDataOrig$gene %in% currData$gene),])
		}
	
		#Make adjustments based on minGene and maxGene
		if(nrow(currData) > maxGeneCurr){
			currData<-currData[seq(maxGeneCurr),]
		}
		if(nrow(currData) < minGeneCurr){
			warning("The selection criteria are too stringent to get the specified minimum number of genes")
		}
	
		#Sort table by FDR	
		currData<-currData[
  			with(currData,order(p_adj_FDR,(1 - avg_logFC))),
		]
		
		#Add genes to list
		genesList[[length(genesList) + 1]]<-as.character(currData$gene)
	}
	
	names(genesList)<-paste("markers",unique(DEGsCons$comparison),sep="_")
	return(genesList)
}


























#New MASTDETest function that fully runs MAST (returning continuous, discrete, and hurdle results)
MASTDETest_expanded<-function(data.use=NULL,cells.1=NULL,cells.2=NULL,latent.vars=NULL,verbose=T){

#Prep input for MAST
  if (length(x = latent.vars) > 0) {
    latent.vars <- scale(x = latent.vars)
  }
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  "%||%" <- function(a, b) if (!is.null(a)) a else b
  latent.vars <- latent.vars %||% group.info
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  latent.vars.names <- c("condition", colnames(x = latent.vars))
  latent.vars <- cbind(latent.vars, group.info)
  latent.vars$wellKey <- rownames(x = latent.vars)
  fdat <- data.frame(rownames(x = data.use))
  colnames(x = fdat)[1] <- "primerid"
  rownames(x = fdat) <- fdat[, 1]

#Make expression object  
  sca <- MAST::FromMatrix(
    exprsArray = as.matrix(x = data.use),
    cData = latent.vars,
    fData = fdat
  )

#Do test
  cond <- factor(x = SummarizedExperiment::colData(sca)$group)
  cond <- relevel(x = cond, ref = "Group2")
  SummarizedExperiment::colData(sca)$condition <- cond
  fmla <- as.formula(
    object = paste0(" ~ ", paste(latent.vars.names, collapse = "+"))
  )
  zlmCond <- MAST::zlm(formula = fmla, sca = sca)

#Pull out summary results
  summaryCond <- summary(object = zlmCond, doLRT = 'conditionGroup1')
  summaryDt <- summaryCond$datatable
  #Discrete
  fcDiscrete <-summaryDt[contrast=='conditionGroup1' & component=='D', .(primerid, `Pr(>Chisq)`, coef, ci.hi, ci.lo)]
  fcDiscrete[,FDR:=p.adjust(`Pr(>Chisq)`, 'BH')]
  #Continuous
  fcContinuous <-summaryDt[contrast=='conditionGroup1' & component=='C', .(primerid, `Pr(>Chisq)`, coef, ci.hi, ci.lo)]
  fcContinuous[,FDR:=p.adjust(`Pr(>Chisq)`, 'BH')]
  #Hurdle
  fcHurdle <- merge(
     summaryDt[contrast=='conditionGroup1' & component=='H', .(primerid, `Pr(>Chisq)`)],
     summaryDt[contrast=='conditionGroup1' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid'
  )
  fcHurdle[,FDR:=p.adjust(`Pr(>Chisq)`, 'BH')]
  
  #Merge and get naming and ordering right
  fcAll<-merge(fcHurdle,merge(fcContinuous,fcDiscrete,by='primerid',all=T,
  	suffixes=c("_cont","_disc")),by='primerid',all=T)
  colnames(fcAll)[1]<-"gene"
  colnames(fcAll)[2]<-"pval_h"
  colnames(fcAll)[3]<-"LFC"
  colnames(fcAll)[4]<-"LFC_ci_hi"
  colnames(fcAll)[5]<-"LFC_ci_lo"
  colnames(fcAll)[6]<-"FDR_h"
  colnames(fcAll)[7]<-"pval_c"
  colnames(fcAll)[8]<-"coef_c"
  colnames(fcAll)[9]<-"coef_ci_hi_c"
  colnames(fcAll)[10]<-"coef_ci_lo_c"
  colnames(fcAll)[11]<-"FDR_c"
  colnames(fcAll)[12]<-"pval_d"
  colnames(fcAll)[13]<-"coef_d"
  colnames(fcAll)[14]<-"coef_ci_hi_d"
  colnames(fcAll)[15]<-"coef_ci_lo_d"
  colnames(fcAll)[16]<-"FDR_d"
  fcAll<-fcAll[,c(1,3:5,2,6,8:10,7,11,13:15,12,16)]
  return(as.data.frame(fcAll))
}




























####### Do Fisher test of enrichment for each cluster. Genes is a vector of genes for which you want to check if they are enriched.
#in the genes listed in the first column of the gene2module. The second column gives the groups/categories to separately test. 
#gene2module colnames must be c("Gene","Module"). 
DoEnrichment<-function(genes,gene2module,background=20441){

	moduleColors<-gene2module$Module
	
	#Create vectors to store results in for each module
	enrichment_hyp_p<-data.frame("pval"=rep(NA, length(unique(moduleColors))),
		"module"=rep(1, length(unique(moduleColors)))) #store hypergeometric pvalues
	enrichment_fish_p<-data.frame("pval"=rep(NA, length(unique(moduleColors))),
		"module"=rep(1, length(unique(moduleColors))))#store fisher pvalues
	enrichment_fish_o<-data.frame("pval"=rep(NA, length(unique(moduleColors))),
		"module"=rep(1, length(unique(moduleColors)))) #store fisher odds ratios

	##Run loop through each module
	for(i in 1:length(unique(moduleColors))){
		
		#Relevant variables
#		A=nrow(gene2module) #number of background genes (either the number of genes in the human hg_38 genome 
					   #(20,441 from biomart) or the number of DEGs (length(degs)).
		A=background #protein coding genes in hg_38 genome
		B=length(gene2module$Gene[gene2module$Module==unique(moduleColors)[i]])## number of genes sampled (
					   #i.e., num genes in this module)
		C=length(genes) #total number of genes in the target list (e.g., GWAS list)
		D=length(intersect(genes,gene2module$Gene[gene2module$Module==unique(moduleColors)[i]])) #number of genes 
		                 #in the target list that intersect with the sample (i.e., module).
		
		####Hypergeometric Test
		#Here are the inputs for the hypergeometric test
		targetsInSample = D #The number of target genes found in the set of selected genes
		targetsInBackground = C #The number of target genes overall (in the subsample and background together)
		backgroundMinusTargets = A-targetsInBackground #The number of genes in the background NOT in the target list
		sampleSize = B #The number of selected genes (i.e., number of genes in the module)
		
		#Now do the hypergeometric test (lower.tail = TRUE if testing for depletion; FALSE if testing for enrichment)
		#Note that when lower.tail is TRUE (default), probabilities are P[X ≤ x]. When lower.tail is FALSE, probabilities
		#are P[X > x], and thus we need to subtract x by 1 to get P[X ≥ x].
		enrichment_hyp_p[i,1]<-phyper((targetsInSample - 1),targetsInBackground,backgroundMinusTargets,
			sampleSize,lower.tail = FALSE)
		
		####Fisher's exact test
		#Here are the inputs for the Fisher's exact test
		targetsInSample = D #The number of target genes found in the set of selected genes
		targetsInBackgroundMinusTargetsInSample = targetsInBackground - targetsInSample #The number of target genes 
			#in the background minus those in the sample
		sampleMinusTargets = sampleSize - targetsInSample #The number of genes in the background NOT in the target list
		backgroundWithOnlyTargetsInSample = backgroundMinusTargets - sampleSize + targetsInSample #background minus the 
			#target genes and minus the genes in the sample. But then add in the number of target genes that are in the sample
		
		#Now fit these into a contingency table
		countsFisher<-matrix(data=c(targetsInSample,targetsInBackgroundMinusTargetsInSample,sampleMinusTargets,
			backgroundWithOnlyTargetsInSample),nrow=2)
		row.names(countsFisher)=c("InTargetList","NotInTargetList")
		colnames(countsFisher)=c("SampleDataset","BackgroundDataset")
	
		#Do Fisher's test (for enrichment do: alternative = 'greater') #store both p-value and odds ratio
		enrichment_fish_p[i,1]<-fisher.test(countsFisher,alternative='greater')$p.value
		enrichment_fish_o[i,1]<-unname(fisher.test(countsFisher)$estimate) 
			#odds ratio is (targetsInSample/targetsInBackgroundMinusTargetsInSample) / 
				#(sampleMinusTargets/backgroundWithOnlyTargetsInSample)
		
		#Add in module info
		enrichment_hyp_p[i,2]<-unique(moduleColors)[i]
		enrichment_fish_p[i,2]<-unique(moduleColors)[i]
		enrichment_fish_o[i,2]<-unique(moduleColors)[i]
	}	
	
#	#Sort
#	enrichment_hyp_p<-enrichment_hyp_p[order(enrichment_hyp_p$pval),]
#	enrichment_fish_p<-enrichment_hyp_p[order(enrichment_hyp_p$pval),]
#	enrichment_fish_o<-enrichment_hyp_p[order(enrichment_hyp_p$pval),]
	
	#Calculate FDR
	enrichment_hyp_p<-cbind(enrichment_hyp_p,"padj"=as.numeric(p.adjust(enrichment_hyp_p$pval,method="fdr")))
	enrichment_hyp_p<-enrichment_hyp_p[,c(2,1,3)]
	enrichment_fish_p<-cbind(enrichment_fish_p,"padj"=as.numeric(p.adjust(enrichment_fish_p$pval,method="fdr")))
	enrichment_hyp_p<-enrichment_fish_p[,c(2,1,3)]
		
	return(enrichment_hyp_p)
}
































#########Function for making dot plots for two compared populations using a Seurat DEG table including a set of genes 

#degs_tab = the DEG output table from Seurat
#colExtremes = extreme colors to use for scaling of fold change
#expMinSet = the minimum fold change to consider (anything less is set to this value)
#expMaxSet = the maximum fold change to consider (anything larger is set to this value)
#pctMinSize = the minimum point size to plot
#pctMaxSize = the maximum point size to plot
#maiVec = the plot margin settings (in inches)
#height = height of plot
#width = width of plot
make2WayDotplot<-function(degs_tab,colExtremes=c("#FFD5CE","#FF0000FF"),expMinSet=0.25,expMaxSet=2,
	pctMinSize=0.1,pctMaxSize=6,maiVec=c(1.02,0.1,0.82,0.2),height=7,width=3){
	
	#Load necessary libraries
	if(!require(plotrix)){
    	install.packages("plotrix")
    	library(plotrix)
	}
	if(!require(scales)){
    	install.packages("scales")
    	library(scales)
	}
	dev.new(height=height,width=width) #new window
	degs_tab<-degs_tab[seq(dim(degs_tab)[1],1),] #flip row order
	par(mai=maiVec) #set margins
	cexVec1<-scales::rescale(degs_tab$pct.1,to=c(pctMinSize,pctMaxSize),from=c(0,1)) #set point 1 size
	cexVec2<-scales::rescale(degs_tab$pct.2,to=c(pctMinSize,pctMaxSize),from=c(0,1)) #set point 2 size
	degs_tab$avg_logFC[which(degs_tab$avg_logFC > expMaxSet)]<-expMaxSet #break at expMaxSet
	degs_tab$avg_logFC[which(degs_tab$avg_logFC < expMinSet)]<-expMinSet #break at expMinSet
	colVec<-color.scale(scales::rescale(degs_tab$avg_logFC,to=c(0,1),from=c(expMinSet,expMaxSet)),extremes=colExtremes) #set color scale
	barplot(rep(nrow(degs_tab)+1,2), col = NA, border = NA, axes = FALSE,space=0.05) #Empty plot
	for(i in 1:nrow(degs_tab)){ #for each gene
		points(c(1,2),rep(i,2),pch=16,cex=c(cexVec1[i],cexVec2[i]),col=colVec[i]) #plot points
		text(0.3,i,labels=degs_tab$gene[i]) #add labels
	}
}

#EXAMPLE: Create plot
#Here's a function for creating dot plots from a DEG table
#degs<-read.table("SF_clusterDEGs_ALL.txt",sep="\t",header=T,stringsAsFactors=F)
#degs<-degs[1:10,]
#make2WayDotplot(degs_tab=degs)

































######################################################################################## MONOCLE FUNCTIONS


#################################### Function for calculating binned values for monocle
getBinnedValues<-function(pData,metadata_col_index,quantile=F,k=5){
	#Sort metadata by pseudotime
	pData_byPseud<-pData[order(pData$Pseudotime),]
	#Establish dataframe
	metadataTab<-data.frame(matrix(nrow=100,ncol=2))
	colnames(metadataTab)<-c(colnames(pData_byPseud)[metadata_col_index],"pseudotime")
	
	#Create bins for categorical variable
	d<-pData_byPseud[,metadata_col_index]
	p<-pData_byPseud$Pseudotime
	bins<-round(seq(1,length(d),by=(length(d) / 100)))
	bins_p<-seq(min(p),max(p),length.out = 100)
	for(i in 1:100){
		#For each bin, get pseudotime midpoint
		if(quantile){ #if you want to create bins based on the data (equal number of cells per bin)
			metadataTab[i,2]<-mean(range(p[bins[i]:(round(bins[i] + (length(d) / 100)) - 1)]))
		}else{ #if you want to bin the way Monocle does, i.e., bins at equal intervals across pseudotime and takes the midpoint
		metadataTab[i,2]<-mean(c(bins_p[i],(bins_p[i] + (bins_p[2] - bins_p[1]))))
		}
		#For each bin, get dominant categorial value
		if(quantile){#if you want to create bins based on the data (equal number of cells per bin)
			cellTab<-sort(table(d[bins[i]:(round(bins[i] + (length(d) / 100)) - 1)]),decreasing = T)
			metadataTab[i,1]<-sample(names(cellTab)[which(cellTab == max(cellTab))],1)
			#if you want to bin the way Monocle does, i.e., bins at equal intervals across pseudotime and takes the 
			#most representative of k nearest neighbor cells
		}else{
			cellTab<-sort(table(pData_byPseud[order(abs(pData_byPseud$Pseudotime - metadataTab[i,2])),metadata_col_index][1:k]))
			metadataTab[i,1]<-sample(names(cellTab)[which(cellTab == max(cellTab))],1)
		}
	}
	return(metadataTab)
}







######################## Deconstructing the plot_genes_in_pseudotime function 

#This function takes the plot_genes_in_pseudotime and gives it more options that I want
#Most importantly, it takes a list of datasets (cds_subsets) so that points and curves can be overlaid onto the same plot.
#I can also take a set of colors for points from each dataset and category 
#	(there need to be as many colors as there are number of categories to color x number of datasets.
#Print_points specifies whether or not you want points printed at all
#Curve_cols gives colors used for each curve (in order of the datasets in the list).
#Fix_y gives whether or not the y-axis should be fixed across all genes or allowed to adjust to the expression of each gene
#overlay gives whether you want to overlay the same genes from different datasets ("datasets") or different genes from the same dataset ("genes") - in the case
	#of different genes from the same dataset, you can only plot one dataset at a time. So cds_subsets should just be a list of the same dataset with different genes subsetted
	#that will all be overlaid onto a single plot. Also, I don't really have print_points = TRUE working for this case. 
#Finally, zeroToOneScale = TRUE will scale all curves between 0 and 1. If false, this is not done
prepPseudotimePlotDatasets<-function(cds_subsets,color_by = NULL,color_vec = NULL,nrow = 4,ncol = 4,print_points=TRUE,
	curve_cols = NULL,fix_y = T,overlay="datasets",min_expr = NULL,cell_size = 0.75,panel_order = NULL,trend_formula = "~ sm.ns(Pseudotime, df=3)",
	label_by_short_name = TRUE,relative_expr = FALSE,vertical_jitter = NULL, horizontal_jitter = NULL,zeroToOneScale=T){
	
	#Set of up list of processed datasets
	cds_expr_final<-list()
	
	#For each dataset to overlay, get ready for plotting
	for(i in 1:length(cds_subsets)){
		cds_subset<-cds_subsets[[i]]
    	f_id <- NA
    	Cell <- NA
    	if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", 
    	    "negbinomial.size")) {
    	    integer_expression <- TRUE
    	}else{
    	    integer_expression <- FALSE
    	    relative_expr <- TRUE
    	}
    	
    	#Get size factor adjusted expression
    	if (integer_expression) {
    	    cds_exprs <- exprs(cds_subset)
    	    if (relative_expr) {
    	        if (is.null(sizeFactors(cds_subset))) {
    	            stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
    	        }
    	        cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
    	    }
    	    cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
    	}else{
    	    cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
    	}
    	
    	#Set min expression to the lowest detection limit
    	if (is.null(min_expr)) {
    	    min_expr <- cds_subset@lowerDetectionLimit
    	}
    	
    	#Add in meta.data
    	colnames(cds_exprs) <- c("f_id", "Cell", "expression")
    	cds_pData <- pData(cds_subset)
    	cds_fData <- fData(cds_subset)
    	cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
    	cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
    	#Set adjusted_expression
    	if (integer_expression) {
    	    cds_exprs$adjusted_expression <- cds_exprs$expression
    	}else{
    	    cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
    	}
    	#Set feature_label
    	if (label_by_short_name == TRUE) {
    	    if (is.null(cds_exprs$gene_short_name) == FALSE) {
    	        cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
    	        cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    	    }else{
    	        cds_exprs$feature_label <- cds_exprs$f_id
    	    }
    	}else{
    	    cds_exprs$feature_label <- cds_exprs$f_id
    	}
    	cds_expr_final[[length(cds_expr_final) + 1]]<-cds_exprs
    }
    	    	
    #Get smooth curve for expression of each gene across pseudotime (but note cells not yet ordered by psudotime)
    #and merge these with the expression + metadata table	
    for(i in 1:length(cds_subsets)){
    	cds_subset<-cds_subsets[[i]]
    	cds_expr_final[[i]]$f_id <- as.character(cds_expr_final[[i]]$f_id) 
    	cds_expr_final[[i]]$feature_label <- factor(cds_expr_final[[i]]$feature_label)
    	new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
    	model_expectation <- genSmoothCurves(cds_subset, cores = 1, 
    	    trend_formula = trend_formula, relative_expr = relative_expr, new_data = new_data)
    	colnames(model_expectation) <- colnames(cds_subset)
    	expectation <- ddply(cds_expr_final[[i]], .(f_id, Cell), function(x) data.frame(expectation = model_expectation[x$f_id, 
    	    x$Cell]))
    	
    	#If overlaying genes, move the curves to be between 0 and 1
    	if(overlay == "genes"){
    		if(zeroToOneScale == T){
    			expectation$expectation<-(expectation$expectation - min(expectation$expectation)) / 
    				(max(expectation$expectation) - min(expectation$expectation))
    		}else{
    			expectation$expectation<-log(expectation$expectation)
    		}
    		#Another option to try is to scale and center each curve	
 #   		expectation$expectation=expectation$expectation - mean(expectation$expectation)
#  			expectation$expectation=expectation$expectation / sd(expectation$expectation)
    	}
    	cds_expr_final[[i]] <- merge(cds_expr_final[[i]], expectation)
    	
    	#Set expression below the specified minimum expression value to the minimum expression value
    	cds_expr_final[[i]]$expression[cds_expr_final[[i]]$expression < min_expr] <- min_expr
    	
    	#Do the same for the expected expression value
    	cds_expr_final[[i]]$expectation[cds_expr_final[[i]]$expectation < min_expr] <- min_expr
    	if (is.null(panel_order) == FALSE) {
    	    cds_expr_final[[i]]$feature_label <- factor(cds_expr_final[[i]]$feature_label, 
    	        levels = panel_order)
    	}
    }
 

    #########Now that we have the list of overlaid datasets prepped, plot each of them using this:
    
    #If there is more than one dataset, make the categorical variable specific to each dataset
    for(i in 1:length(cds_subsets)){
    	categorical_index<-grep(color_by,colnames(cds_expr_final[[i]]))
    	cds_expr_final[[i]][,categorical_index]<-paste(cds_expr_final[[i]][,categorical_index],"_",i,sep="")
    }
    
    #Then establish an empty plot with y axis fitting the global expression range (across all genes) and with the
    #x axis fitting the pseudotime range (across all datasets)
    appendedData<-data.frame()
    for(i in 1:length(cds_subsets)){
    	appendedData<-rbind(appendedData,cds_expr_final[[i]])
    }
 
    #Now make a ggplot
    q <- ggplot(aes(Pseudotime, expression), data = appendedData)

	#Make scatter plot of expression against pseudotime across all genes 
	#cell_size gives the size of the points
	#color_by gives the categories to color by
	#I've added point color specification here
	#Can skip this step if you don't want to include expression points
	if(print_points){
	    if(is.null(color_by) == FALSE){
		     q <- q + geom_point(aes_string(color = color_by), size = I(cell_size), alpha = 0.3, shape = 16,
		        position = position_jitter(horizontal_jitter, vertical_jitter),data = appendedData) +
		        scale_color_manual(values=color_vec[order(unique(appendedData[,categorical_index]))])
		}else{
		    q <- q + geom_point(size = I(cell_size), position = position_jitter(horizontal_jitter, 
	    	  	vertical_jitter),data = appendedData)
	    }
	}
    
    #Plot the expression expectation curves across the combined plot, one for each dataset
 	for(i in 1:length(cds_expr_final)){
	    q <- q + geom_line(aes(x = Pseudotime, y = expectation), color = curve_cols[i], #size = 1,
        	data = cds_expr_final[[i]])
    }
        
    #Now break up the plots based on gene (if overlaying onto different datasets)
    #scale_y_log10 places expression on a log scale
    #scales = "free_y" (default) lets the y-axis shift to fit the data
    #scales = "fixed" keeps the x and y axes the same across plots
    if(fix_y){
    	if(overlay == "datasets"){
 		   	q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow, 
        		ncol = ncol, scales = "fixed")
        }
    }else{
    	if(overlay == "datasets"){
 		   	q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow, 
        		ncol = ncol, scales = "free_y")
        }
    }
    
    if (min_expr < 1) {
        q <- q + expand_limits(y = c(min_expr, 1))
    }
    
    #Add labels
    if (relative_expr) {
        q <- q + ylab("Relative Expression")
    }else{
        q <- q + ylab("Absolute Expression")
    }
    if(overlay == "genes"){
    	q <- q + ylab("Scaled Expression")
    }
    q <- q + xlab("Pseudotime")
 
    #Get rid of ggplot background crap
    q <- q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_blank(), axis.line = element_line(colour = "black"))
	q
}

















#######Function for plotting separate curves for each gene with confidence intervals
plotIndividualCurves<-function(cds_subset,relative_expr=T,min_expr=0.01,colorVec="black",plot_points=F,zeroToOneScale=T,plot_interval=T){
	
	#Get scaled expression
	cds_exprs <- exprs(cds_subset)
    if(relative_expr) {
    	if(is.null(sizeFactors(cds_subset))) {
                stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
        }
        cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
    }
    cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
    
    #Set min expression to the lowest detection limit
    if (is.null(min_expr)) {
        min_expr <- cds_subset@lowerDetectionLimit
    }
    
    #Get fitted values
    f_expression<-cds_exprs$value
    fittedVals <- VGAM::vglm(as.formula(paste("f_expression","~sm.ns(cds_subset$Pseudotime,df=3)",sep="")),family = "negbinomial", epsilon=1e-1)
    seVals<-predictvglm(fittedVals,se.fit=T)
	
	#Get values scaled between zero and 1 if desired
	#Note I'm scaling the points and curves separately since I can't easily estimate the curves using the negative binomial after scaling
	#as scaling introduces non-integers
	if(zeroToOneScale==T){
		#Get scaled values
   		meanVals<-(fittedVals@fitted.values - min(fittedVals@fitted.values)) /
   			(max(fittedVals@fitted.values) - min(fittedVals@fitted.values))
   		lowerVals<-((fittedVals@fitted.values - seVals$se.fit[,1]) - min((fittedVals@fitted.values - seVals$se.fit[,1]))) /
   			(max((fittedVals@fitted.values - seVals$se.fit[,1])) - min((fittedVals@fitted.values - seVals$se.fit[,1])))
   		upperVals<-((fittedVals@fitted.values + seVals$se.fit[,1]) - min((fittedVals@fitted.values + seVals$se.fit[,1]))) /
   			(max((fittedVals@fitted.values + seVals$se.fit[,1])) - min((fittedVals@fitted.values + seVals$se.fit[,1])))
   		pointVals<-(cds_exprs$value - min(cds_exprs$value)) /
   			(max(cds_exprs$value) - min(cds_exprs$value))
   }else{
		meanVals<-fittedVals@fitted.values
		lowerVals<-fittedVals@fitted.values - seVals$se.fit[,1]
		upperVals<-fittedVals@fitted.values + seVals$se.fit[,1]
		pointVals<-cds_exprs$value
	}
   
    #Apply min_expr
   	meanVals[which(meanVals < min_expr)]<-min_expr
   	lowerVals[which(lowerVals < min_expr)]<-lowerVals
   	upperVals[which(upperVals < min_expr)]<-upperVals
   
    #Plot 
    #sqrt of non-scaled expression with points
    if(plot_points==T & zeroToOneScale==F){
    	plot(pData(cds_subset)$Pseudotime,sqrt(pointVals),las=1,ylab="",xlab="",bty="l",main=cds_exprs[1,1],cex.main=2,font.main=3,
    		col=colorVec,cex.axis=1.5,ylim=c(min(sqrt(pointVals)),(max(sqrt(pointVals)) + (max(sqrt(pointVals)) * 0.1))))	
    	points(x=cds_subset$Pseudotime, y=sqrt(meanVals),cex=0.1,pch=20)
    	if(plot_interval){
			points(x=cds_subset$Pseudotime, y=sqrt(upperVals),cex=0.1,col="grey")
			points(x=cds_subset$Pseudotime, y=sqrt(lowerVals),cex=0.1,col="grey")
		}
	}
	#Scaled expression with points
	if(plot_points==T & zeroToOneScale==T){
    	plot(pData(cds_subset)$Pseudotime,pointVals,las=1,ylab="Scaled expression",xlab="Pseudotime",bty="l",main=cds_exprs[1,1],
    		col=colorVec,cex.axis=0.9,ylim=c(min(pointVals),(max(pointVals) + (max(pointVals) * 0.1))))	
    	points(x=cds_subset$Pseudotime, y=meanVals,cex=0.1,las=1,ylab="Relative expression",xlab="Pseudotime",bty="l",main=cds_exprs[1,1])
    	if(plot_interval){
			points(x=cds_subset$Pseudotime, y=upperVals,cex=0.1,col="grey")
			points(x=cds_subset$Pseudotime, y=lowerVals,cex=0.1,col="grey")
		}
	}
	#Non-scaled expression without points
	if(plot_points==F & zeroToOneScale==F){
	    plot(x=cds_subset$Pseudotime, y=meanVals,cex=0.1,las=1,ylab="Relative expression",xlab="Pseudotime",bty="l",
    		main=cds_exprs[1,1],cex.axis=0.9)
    	if(plot_interval){
			points(x=cds_subset$Pseudotime, y=v,cex=0.1,col="grey")
			points(x=cds_subset$Pseudotime, y=lowerVals,cex=0.1,col="grey")
		}
	}
	#Scaled expression without points
	if(plot_points==F & zeroToOneScale==T){
		plot(x=cds_subset$Pseudotime, y=meanVals,cex=0.1,las=1,ylab="Scaled expression",xlab="Pseudotime",bty="l",
			main=cds_exprs[1,1],cex.axis=0.9)
		if(plot_interval){
			points(x=cds_subset$Pseudotime, y=upperVals,cex=0.1,col="grey")
			points(x=cds_subset$Pseudotime, y=lowerVals,cex=0.1,col="grey")
		}
	}  
}

















############# Getting positions of gene clusters in the heat map for use in affinity designer
#So when plotting heat maps in Monocle without clusters, there is no break between the clusters
#So, the below function returns the y coordinates (in points) to place the divider between each cluster in the heat map
#To use this, input the pseudotime DEG table, with clusters identified in a column called "cluster".
#Secondly, position the heat map in an Affinity Designer page such that the top is at points = 0 and the bottom of the
#heat map is at point "bottomLocation".

getCoordinateVec<-function(DEG_table, bottomLocation){
	position<-c(table(DEG_table$cluster)[1])
	position_prop<-c(table(DEG_table$cluster)[1] / nrow(DEG_table) * bottomLocation)
	for(i in 2:length(unique(DEG_table$cluster))){
		position<-append(position,
			(position[length(position)] + table(DEG_table$cluster)[i]))
		position_prop<-append(position_prop,
			(position[(length(position) - 1)] + table(DEG_table$cluster)[i]) / nrow(DEG_table) * bottomLocation)
		names(position)[length(position)]<-names(table(DEG_table$cluster))[i]
		names(position_prop)[(length(position))]<-names(table(DEG_table$cluster))[i]
	}
	return(position_prop)
}