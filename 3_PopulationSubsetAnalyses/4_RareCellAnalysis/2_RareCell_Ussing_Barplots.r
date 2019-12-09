############## Re-making barplots for Kate's qPCR and Ussing data

library(ggplot2)
library(openxlsx)
source("../CommonFunctions.r")

#Bring in all of Kate's data
dat<-readExcelTabs("RareCellData_forStats.xlsx")





############################ qPCR data #############################
#Isolate qpcr
qpcr<-dat[[3]]

#We want to exclude T133 from all the plots
#qpcr<-qpcr[-which(qpcr$Donor == "T133"),]

#Change '2^-(goi-hk)' to value
names(qpcr)[4]<-"value"

#Simplify the sample column so just lists inserts 1-3
qpcr$sample<-sapply(strsplit(qpcr$sample," "),function(x)x[length(x)])



library(lme4)
qpcr_norm<-data.frame()
#Now, for each gene
for(i in 1:length(unique(qpcr$gene))){
	#For each KO status
	for(j in 2:length(unique(qpcr$KO_status))){
		
		#Get estimated values and standard error, accounting for different donors
		#Treat donor as random effect, since the three inserts are non-independent
		#KO test
		currKO<-qpcr[which(qpcr$gene == unique(qpcr$gene)[i] & qpcr$KO_status == unique(qpcr$KO_status)[j]),]
		currlmer<-lmer(currKO$value~(1 | currKO$Donor))
		estimate_KO<-coef(summary(currlmer))[1]
		SE_KO<-coef(summary(currlmer))[2]
		
		#Scramble test
		currSCBL<-qpcr[which(qpcr$gene == unique(qpcr$gene)[i] & qpcr$KO_status == unique(qpcr$KO_status)[1]),]
		currlmerSCBL<-lmer(currSCBL$value~(1 | currSCBL$Donor))
		estimate_SCBL<-coef(summary(currlmerSCBL))[1]
		SE_SCBL<-coef(summary(currlmerSCBL))[2]
		
		#Now normalize the KO values to the Scramble values
		normEstimate_KO<-estimate_KO/estimate_SCBL
		normSE_KO<-SE_KO/estimate_SCBL
		
		#Now add these values to master table
		currNorm<-data.frame("KO_status"=currKO$KO_status[1],"gene"=currKO$gene[1],"value"=normEstimate_KO,"SE"=normSE_KO)
		qpcr_norm<-rbind(qpcr_norm,currNorm)
		
#		################Extra, if want to do the old way, by calculating the mean of the means
#		#For each donor
#		for(k in 1:length(unique(qpcr$Donor))){
#			#Get data for KO
#			currData<-qpcr[which(qpcr$gene == unique(qpcr$gene)[i] & qpcr$KO_status == unique(qpcr$KO_status)[j] &
#				qpcr$Donor == unique(qpcr$Donor)[k]),]
#			#Get data from scramble
#			currDataSCBL<-qpcr[which(qpcr$gene == unique(qpcr$gene)[i] & qpcr$KO_status == unique(qpcr$KO_status)[1] &
#				qpcr$Donor == unique(qpcr$Donor)[k]),]
#			#Get ratio of each KO treatment to the scramble
#			valueRatio<-mean(currData$value / mean(currDataSCBL$value))
#			#Let's keep track of the standard error across the three reps, also scaled to the mean of the three scramble values

	}		
}

### Make barplot
#library(gplots)
#pdf("RareCell_qPCR_barplots.pdf")
dev.new(height=5,width=5)
barplot2(rev(qpcr_norm$value - 1),beside = T,las=1,col=c("gold","darkolivegreen3"),ylim=c(1,22),names.arg="",
	xlab="mRNA expression normalized to baseline",horiz=T, xlim=c(-1,1.3),ci.l = rev((qpcr_norm$value - 1) - qpcr_norm$SE), 
	ci.u = rev((qpcr_norm$value - 1) + qpcr_norm$SE), plot.ci=T, ci.color="grey45",xaxt="n",
	space=rep(c(1,0),(length(qpcr_norm$value) / 2)),width=0.7)
abline(v=0)
axis(1,c(-1,-0.5,0,0.5,1),labels=c(0,0.5,1,1.5,2))
mtext(as.character(qpcr_norm$gene)[1],3,line= -2, at=-1.2)
mtext(as.character(qpcr_norm$gene)[3],3,line= -3.5, at=-1.2)
mtext(as.character(qpcr_norm$gene)[5],3,line= -5, at=-1.2)
mtext(as.character(qpcr_norm$gene)[7],3,line= -6.5, at=-1.2)
mtext(as.character(qpcr_norm$gene)[9],3,line= -8, at=-1.2)
mtext(as.character(qpcr_norm$gene)[11],3,line= -9.5, at=-1.2)
mtext(as.character(qpcr_norm$gene)[13],3,line= -11, at=-1.2)
mtext(as.character(qpcr_norm$gene)[15],3,line= -12.5, at=-1.2)
mtext(as.character(qpcr_norm$gene)[17],3,line= -14, at=-1.2)
mtext(as.character(qpcr_norm$gene)[19],3,line= -15.5, at=-1.2)




















############################ Ussing data #############################
################# First, run each stim, 1 at a time
#Isolate qpcr
ussing<-dat[[4]]

#We want to exclude T133 from all the plots
#ussing<-ussing[-which(ussing$Donor == "T133"),]

#Change '2^-(goi-hk)' to value
names(ussing)[5]<-"value"

#Now ussing needs to be a subset of one of the Ussing parameters
#Re-run for each of the four using parameters
ussing<-ussing[which(ussing$Ussing_parameter == "Voltage"),]
ussing<-ussing[which(ussing$Ussing_parameter == "Current"),]
ussing<-ussing[which(ussing$Ussing_parameter == "Resistance"),]
ussing<-ussing[which(ussing$Ussing_parameter == "Conductance"),]

library(lme4)
ussing_norm<-data.frame()
#Now, for each stim
for(i in 1:length(unique(ussing$stim))){
	#For each KO status
	for(j in 2:length(unique(ussing$KO_status))){
		
		#Get estimated values and standard error, accounting for different donors
		#Treat donor as random effect, since the three inserts are non-independent
		#KO test
		currKO<-ussing[which(ussing$stim == unique(ussing$stim)[i] & ussing$KO_status == unique(ussing$KO_status)[j]),]
		currlmer<-lmer(currKO$value~(1 | currKO$Donor))
		estimate_KO<-coef(summary(currlmer))[1]
		SE_KO<-coef(summary(currlmer))[2]
		
		#Scramble test
		currSCBL<-ussing[which(ussing$stim == unique(ussing$stim)[i] & ussing$KO_status == unique(ussing$KO_status)[1]),]
		currlmerSCBL<-lmer(currSCBL$value~(1 | currSCBL$Donor))
		estimate_SCBL<-coef(summary(currlmerSCBL))[1]
		SE_SCBL<-coef(summary(currlmerSCBL))[2]
		
		#Now normalize the KO values to the Scramble values
		normEstimate_KO<-estimate_KO/estimate_SCBL
		normSE_KO<-abs(SE_KO/estimate_SCBL)
		
		#Now add these values to master table
		currNorm<-data.frame("KO_status"=currKO$KO_status[1],"stim"=currKO$stim[1],"value"=normEstimate_KO,"SE"=normSE_KO)
		ussing_norm<-rbind(ussing_norm,currNorm)
	}		
}

#library(gplots)
#pdf("RareCell_ussing_barplots_voltage.pdf")
dev.new(height=3.5,width=5)
barplot2(rev(ussing_norm$value - 1),beside = T,las=1,col=c("gold","darkolivegreen3"),ylim=c(1,11),names.arg="",
	xlab="Ussing value normalized to baseline",horiz=T, xlim=c(-1,2),ci.l = rev((ussing_norm$value - 1) - ussing_norm$SE), 
	ci.u = rev((ussing_norm$value - 1) + ussing_norm$SE), plot.ci=T, ci.color="grey45",xaxt="n",
	space=rep(c(1,0),(length(ussing_norm$value) / 2)),width=0.7,main=ussing$Ussing_parameter[1])
abline(v=0)
axis(1,c(-1,-0.5,0,0.5,1,1.5,2),labels=c(0,0.5,1,1.5,2,2.5,3))
mtext(as.character(ussing_norm$stim)[1],3,line= -2, at=-1)
mtext(as.character(ussing_norm$stim)[3],3,line= -3.5, at=-1)
mtext(as.character(ussing_norm$stim)[5],3,line= -5, at=-1)
mtext(as.character(ussing_norm$stim)[7],3,line= -6.5, at=-1)
mtext(as.character(ussing_norm$stim)[9],3,line= -8, at=-1)





################# Second, plot all four stims on the same plot on the same plot
#Isolate qpcr
ussing<-dat[[4]]

##We want to exclude T133 from all the plots
#ussing<-ussing[-which(ussing$Donor == "T133"),]

#We want to exclude current from all the plots
#ussing<-ussing[-which(ussing$Ussing_parameter == "Current"),]

#Change '2^-(goi-hk)' to value
names(ussing)[5]<-"value"

library(lme4)
ussing_norm<-data.frame()
#For each ussing paramter
for(h in 1:length(unique(ussing$Ussing_parameter))){
	#Now, for each stim
	for(i in 1:length(unique(ussing$stim))){
		#For each KO status
		for(j in 2:length(unique(ussing$KO_status))){
			
			#Get estimated values and standard error, accounting for different donors
			#Treat donor as random effect, since the three inserts are non-independent
			#KO test
			currKO<-ussing[which(ussing$Ussing_parameter == unique(ussing$Ussing_parameter)[h] & 
				ussing$stim == unique(ussing$stim)[i] & ussing$KO_status == unique(ussing$KO_status)[j]),]
			currlmer<-lmer(currKO$value~(1 | currKO$Donor))
			estimate_KO<-coef(summary(currlmer))[1]
			SE_KO<-coef(summary(currlmer))[2]
			
			#Scramble test
			currSCBL<-ussing[which(ussing$Ussing_parameter == unique(ussing$Ussing_parameter)[h] &
				ussing$stim == unique(ussing$stim)[i] & ussing$KO_status == unique(ussing$KO_status)[1]),]
			currlmerSCBL<-lmer(currSCBL$value~(1 | currSCBL$Donor))
			estimate_SCBL<-coef(summary(currlmerSCBL))[1]
			SE_SCBL<-coef(summary(currlmerSCBL))[2]
			
			#Now normalize the KO values to the Scramble values
			normEstimate_KO<-estimate_KO/estimate_SCBL
			normSE_KO<-abs(SE_KO/estimate_SCBL)
			
			#Now add these values to master table
			currNorm<-data.frame("Ussing_parameter"=currKO$Ussing_parameter[1],"KO_status"=currKO$KO_status[1],
				"stim"=currKO$stim[1],"value"=normEstimate_KO,"SE"=normSE_KO)
			ussing_norm<-rbind(ussing_norm,currNorm)
		}		
	}
}

#library(gplots)
pdf("RareCell_ussing_barplots_allParameters.pdf",height=8,width=5)
#dev.new(height=8,width=5)
barplot2(rev(ussing_norm$value - 1),beside = T,las=1,col=c("gold","darkolivegreen3"),ylim=c(1,42),names.arg="",
	xlab="Ussing value normalized to baseline",horiz=T, xlim=c(-0.5,1.5),ci.l = rev((ussing_norm$value - 1) - ussing_norm$SE), 
	ci.u = rev((ussing_norm$value - 1) + ussing_norm$SE), plot.ci=T, ci.color="grey45",xaxt="n",
	space=rep(c(1,0),(length(ussing_norm$value) / 2)),width=0.7,main=ussing$Ussing_parameter[1])
abline(v=0)
axis(1,c(-0.5,0,0.5,1,1.5),labels=c(0.5,1,1.5,2,2.5))
dev.off()






















############################ Ussing clampcurrent data #############################
################# First, run just the APT stim (for the paper figure)
#Isolate qpcr
ussingClamp<-dat[[5]]

#We want to exclude T133 from all the plots
#ussingClamp<-ussingClamp[-which(ussingClamp$Donor == "T133"),]

#We only want the ATP stim, for now
ussingClamp<-ussingClamp[which(ussingClamp$stim == "ATP"),]

#Change '2^-(goi-hk)' to value
names(ussingClamp)[5]<-"value"


library(lme4)
ussingClamp_norm<-data.frame()
#Now, for each stim
for(i in 1:length(unique(ussingClamp$stim))){
	#For each KO status
	for(j in 2:length(unique(ussingClamp$KO_status))){
		
		#Get estimated values and standard error, accounting for different donors
		#Treat donor as random effect, since the three inserts are non-independent
		#KO test
		currKO<-ussingClamp[which(ussingClamp$KO_status == unique(ussingClamp$KO_status)[j]),]
		#Scramble test
		currSCBL<-ussingClamp[which(ussingClamp$KO_status == unique(ussingClamp$KO_status)[1]),]
		#This time, I want to normalize ahead of time, since each value in the KO should be normalized to the value from the same batch in the scramble
		#Each scramble measure is paired with a KO measure
		currKO$value<-currKO$value/currSCBL$value
		currlmer<-lmer(currKO$value~1 + (1 | currKO$Donor))
		normEstimate_KO<-coef(summary(currlmer))[1]
		normSE_KO<-coef(summary(currlmer))[2]
		
		#Now add these values to master table
		currNorm<-data.frame("KO_status"=currKO$KO_status[1],"stim"=currKO$stim[1],"value"=normEstimate_KO,"SE"=normSE_KO)
		ussingClamp_norm<-rbind(ussingClamp_norm,currNorm)
	}		
}

#library(gplots)
#pdf("RareCell_ussingClampedATP_barplots.pdf")
dev.new(height=2,width=2)
barplot2(rev(ussingClamp_norm$value - 1),las=1,col=c("gold","darkolivegreen3"),ylim=c(1,2),names.arg="",
	xlab="Clamped current spike (ATP) normalized to baseline",horiz=T, xlim=c(-1,2),ci.l = rev((ussingClamp_norm$value - 1) - ussingClamp_norm$SE), 
	ci.u = rev((ussingClamp_norm$value - 1) + ussingClamp_norm$SE), plot.ci=T, ci.color="grey45",xaxt="n",space=c(0,2),
	width=0.7,main=ussingClamp$ussingClamp_parameter[1])
abline(v=0)
axis(1,c(-1,-0.5,0,0.5,1,1.5),labels=c(0,0.5,1,1.5,2,2.5))






################# Now, just do all the stims together
#Isolate qpcr
ussingClamp<-dat[[5]]

#We want to exclude T133 from all the plots
#ussingClamp<-ussingClamp[-which(ussingClamp$Donor == "T133"),]

#Change '2^-(goi-hk)' to value
names(ussingClamp)[5]<-"value"


library(lme4)
ussingClamp_norm<-data.frame()
#Now, for each stim
for(i in 1:length(unique(ussingClamp$stim))){
	#For each KO status
	for(j in 2:length(unique(ussingClamp$KO_status))){
		
		#Get estimated values and standard error, accounting for different donors
		#Treat donor as random effect, since the three inserts are non-independent
		#KO test
		currKO<-ussingClamp[which(ussingClamp$stim == unique(ussingClamp$stim)[i] & ussingClamp$KO_status == unique(ussingClamp$KO_status)[j]),]
		#Scramble test
		currSCBL<-ussingClamp[which(ussingClamp$stim == unique(ussingClamp$stim)[i] & ussingClamp$KO_status == unique(ussingClamp$KO_status)[1]),]
		#This time, I want to normalize ahead of time, since each value in the KO should be normalized to the value from the same batch in the scramble
		#Each scramble measure is paired with a KO measure
		currKO$value<-currKO$value/currSCBL$value
		currlmer<-lmer(currKO$value~1 + (1 | currKO$Donor))
		normEstimate_KO<-coef(summary(currlmer))[1]
		normSE_KO<-coef(summary(currlmer))[2]
		
		#Now add these values to master table
		currNorm<-data.frame("KO_status"=currKO$KO_status[1],"stim"=currKO$stim[1],"value"=normEstimate_KO,"SE"=normSE_KO)
		ussingClamp_norm<-rbind(ussingClamp_norm,currNorm)
	}		
}

#library(gplots)
#pdf("RareCell_ussingClampedALL_barplots.pdf")
dev.new(height=3.5,width=3)
barplot2(rev(ussingClamp_norm$value - 1),las=1,col=c("gold","darkolivegreen3"),ylim=c(1,11),names.arg="",
	xlab="Clamped current spike normalized to baseline",horiz=T, xlim=c(-1,5),ci.l = rev((ussingClamp_norm$value - 1) - ussingClamp_norm$SE), 
	ci.u = rev((ussingClamp_norm$value - 1) + ussingClamp_norm$SE), plot.ci=T, ci.color="grey45",xaxt="n",space=rev(c(0,1,0,1,0,1,0,1,0,1)),
	width=0.7,main=ussingClamp$ussingClamp_parameter[1])
abline(v=0)
axis(1,c(-1,0,1,2,3,4),labels=c(0,1,2,3,4,5))





























####################### Finally, do for quantification of GRP and FOXI1 imaging

#######GRP
#Isolate data
imaging<-dat[[1]]

##We want to exclude T133 from all the plots
#imaging<-imaging[-which(imaging$Donor == "T133"),]

#Change '2^-(goi-hk)' to value
names(imaging)[3]<-"value"

library(lme4)
imaging_norm<-data.frame()
#For each KO status
for(j in 2:length(unique(imaging$KO_status))){
	
	#Get estimated values and standard error, accounting for different donors
	#Treat donor as random effect, since the three inserts are non-independent
	#KO test
	currKO<-imaging[which(imaging$KO_status == unique(imaging$KO_status)[j]),]
	currlmer<-lmer(currKO$value~(1 | currKO$Donor))
	estimate_KO<-coef(summary(currlmer))[1]
	SE_KO<-coef(summary(currlmer))[2]
	
	#Scramble test
	currSCBL<-imaging[which(imaging$KO_status == unique(imaging$KO_status)[1]),]
	currlmerSCBL<-lmer(currSCBL$value~(1 | currSCBL$Donor))
	estimate_SCBL<-coef(summary(currlmerSCBL))[1]
	SE_SCBL<-coef(summary(currlmerSCBL))[2]
	
	#Now normalize the KO values to the Scramble values
	normEstimate_KO<-estimate_KO/estimate_SCBL
	normSE_KO<-abs(SE_KO/estimate_SCBL)
	
	#Now add these values to master table
	currNorm<-data.frame("KO_status"=currKO$KO_status[1],"value"=normEstimate_KO,"SE"=normSE_KO)
	imaging_norm<-rbind(imaging_norm,currNorm)
}
GRP_norm<-imaging_norm	

#library(gplots)
#pdf("RareCell_imagingGRP_barplots.pdf",height=8,width=5)
dev.new(height=2,width=2.5)
barplot2(rev(imaging_norm$value - 1),las=1,col=c("gold","darkolivegreen3"),ylim=c(1,2),names.arg="",
	xlab="No. GRP+ cells/field",horiz=T, xlim=c(-1,2.2),ci.l = rev((imaging_norm$value - 1) - imaging_norm$SE), 
	ci.u = rev((imaging_norm$value - 1) + imaging_norm$SE), plot.ci=T, ci.color="grey45",xaxt="n",space=c(0,2),
	width=0.7,main=imaging$imaging_parameter[1])
abline(v=0)
axis(1,c(-1,0,1),labels=c(0,1,2))





###########FOXI1
#Isolate data
imaging<-dat[[2]]

##We want to exclude T133 from all the plots
#imaging<-imaging[-which(imaging$Donor == "T133"),]

#Change '2^-(goi-hk)' to value
names(imaging)[3]<-"value"

library(lme4)
imaging_norm<-data.frame()
#For each KO status
for(j in 2:length(unique(imaging$KO_status))){
	
	#Get estimated values and standard error, accounting for different donors
	#Treat donor as random effect, since the three inserts are non-independent
	#KO test
	currKO<-imaging[which(imaging$KO_status == unique(imaging$KO_status)[j]),]
	currlmer<-lmer(currKO$value~(1 | currKO$Donor))
	estimate_KO<-coef(summary(currlmer))[1]
	SE_KO<-coef(summary(currlmer))[2]
	
	#Scramble test
	currSCBL<-imaging[which(imaging$KO_status == unique(imaging$KO_status)[1]),]
	currlmerSCBL<-lmer(currSCBL$value~(1 | currSCBL$Donor))
	estimate_SCBL<-coef(summary(currlmerSCBL))[1]
	SE_SCBL<-coef(summary(currlmerSCBL))[2]
	
	#Now normalize the KO values to the Scramble values
	normEstimate_KO<-estimate_KO/estimate_SCBL
	normSE_KO<-abs(SE_KO/estimate_SCBL)
	
	#Now add these values to master table
	currNorm<-data.frame("KO_status"=currKO$KO_status[1],"value"=normEstimate_KO,"SE"=normSE_KO)
	imaging_norm<-rbind(imaging_norm,currNorm)
}	
FOXI1_norm<-imaging_norm		

#library(gplots)
#pdf("RareCell_imagingFOXI1_barplots.pdf",height=8,width=5)
dev.new(height=2,width=2.5)
barplot2(rev(imaging_norm$value - 1),las=1,col=c("gold","darkolivegreen3"),ylim=c(1,2),names.arg="",
	xlab="No. GRP+ cells/field",horiz=T, xlim=c(-1,2.2),ci.l = rev((imaging_norm$value - 1) - imaging_norm$SE), 
	ci.u = rev((imaging_norm$value - 1) + imaging_norm$SE), plot.ci=T, ci.color="grey45",xaxt="n",space=c(0,2),
	width=0.7,main=imaging$imaging_parameter[1])
abline(v=0)
axis(1,c(-1,-1.5,0,1),labels=c(0,0.5,1,2))


















#Print all the raw results to a table for Kate
imaging_norm<-rbind(cbind(data.frame("KO_status"=GRP_norm$KO_status,stringsAsFactors=F),"gene"="GRP",GRP_norm[,2:3]),
	cbind(data.frame("KO_status"=FOXI1_norm$KO_status,stringsAsFactors=F),"gene"="FOXI1",FOXI1_norm[,2:3]))
tablesToPrint<-list(qpcr_norm,ussing_norm,ussingClamp_norm,imaging_norm)
write.xlsx(tablesToPrint,file="RareCell_barplotData.xlsx",sheetName=c("qPCR","Ussing","UssingClamp","Imaging"),firstRow=T,colWidths="auto")































############################ New Ussing conductance-only data from Kate - 11/20/19 #############################
################# First, run each stim, 1 at a time
#Isolate qpcr
ussing_cond<-dat[[6]]

#We want to exclude T133 from all the plots
#ussing_cond<-ussing_cond[-which(ussing_cond$Donor == "T133"),]

#Change '2^-(goi-hk)' to value
names(ussing_cond)[5]<-"value"

library(lme4)
ussing_cond_norm<-data.frame()
#Now, for each stim
for(i in 1:length(unique(ussing_cond$treatment))){
	#For each KO status
	for(j in 2:length(unique(ussing_cond$KO_status))){
		
		#Get estimated values and standard error, accounting for different donors
		#Treat donor as random effect, since the three inserts are non-independent
		#KO test
		currKO<-ussing_cond[which(ussing_cond$treatment == unique(ussing_cond$treatment)[i] & ussing_cond$KO_status == unique(ussing_cond$KO_status)[j]),]
		currlmer<-lmer(currKO$value~(1 | currKO$donor))
		estimate_KO<-coef(summary(currlmer))[1]
		SE_KO<-coef(summary(currlmer))[2]
		
		#Scramble test
		currSCBL<-ussing_cond[which(ussing_cond$treatment == unique(ussing_cond$treatment)[i] & ussing_cond$KO_status == unique(ussing_cond$KO_status)[1]),]
		currlmerSCBL<-lmer(currSCBL$value~(1 | currSCBL$donor))
		estimate_SCBL<-coef(summary(currlmerSCBL))[1]
		SE_SCBL<-coef(summary(currlmerSCBL))[2]
		
		#Now normalize the KO values to the Scramble values
		normEstimate_KO<-estimate_KO/estimate_SCBL
		normSE_KO<-abs(SE_KO/estimate_SCBL)
		
		#Now add these values to master table
		currNorm<-data.frame("KO_status"=currKO$KO_status[1],"treatment"=currKO$treatment[1],"value"=normEstimate_KO,"SE"=normSE_KO)
		ussing_cond_norm<-rbind(ussing_cond_norm,currNorm)
	}		
}


#library(gplots)
#pdf("RareCell_ussing_cond_barplots.pdf")
dev.new(height=3.5,width=5)
barplot2(rev(ussing_cond_norm$value - 1),beside = T,las=1,col=c("gold","darkolivegreen3"),ylim=c(1,11),names.arg="",
	xlab="ussing_cond value normalized to baseline",horiz=T, xlim=c(-0.5,0.5),ci.l = rev((ussing_cond_norm$value - 1) - ussing_cond_norm$SE), 
	ci.u = rev((ussing_cond_norm$value - 1) + ussing_cond_norm$SE), plot.ci=T, ci.color="grey45",xaxt="n",
	space=rep(c(1,0),(length(ussing_cond_norm$value) / 2)),width=0.7,main=ussing_cond$ussing_cond_parameter[1])
abline(v=0)
axis(1,c(-0.5,0,0.5),labels=c(0.5,1.0,1.5))
mtext(as.character(ussing_cond_norm$treatment)[1],3,line= -2, at=-1)
mtext(as.character(ussing_cond_norm$treatment)[3],3,line= -3.5, at=-1)
mtext(as.character(ussing_cond_norm$treatment)[5],3,line= -5, at=-1)
mtext(as.character(ussing_cond_norm$treatment)[7],3,line= -6.5, at=-1)
mtext(as.character(ussing_cond_norm$treatment)[9],3,line= -8, at=-1)








############# Running a test of difference between KO and Scramble
library(lmerTest) #Changes the lmer function to output Satterthwaite approximated p-values (can also do ANOVA comparing full model (with KO status) and partial model (without)
#Also, test for differences between scramble and KO for each treatment
ussing_cond_test<-data.frame()
for(i in 1:length(unique(ussing_cond$treatment))){
	#For each KO status
	for(j in 2:length(unique(ussing_cond$KO_status))){
		
		#Get estimated values and standard error, accounting for different donors
		#Treat donor as random effect, since the three inserts are non-independent
		#KO test
		currKO<-ussing_cond[which(ussing_cond$treatment == unique(ussing_cond$treatment)[i] & (ussing_cond$KO_status == unique(ussing_cond$KO_status)[j] |
			ussing_cond$KO_status == unique(ussing_cond$KO_status)[1])),]
		currlmer<-lmer(currKO$value~currKO$KO_status + (1 | currKO$donor))
		estimate_diffSCRB<-summary(currlmer)$coefficients[2,1]
		pvalue_diffSCRB<-coef(summary(currlmer))[2,5]
		
		#Now add these values to master table
		currNorm<-data.frame("KO_status"=currKO$KO_status[2],"treatment"=currKO$treatment[1],"estimate_diffSCRB"=estimate_diffSCRB,"pvalue_diffSCRB"=pvalue_diffSCRB)
		ussing_cond_test<-rbind(ussing_cond_test,currNorm)
	}		
}





############# Running a second test looking at decreasing differences between scramble and KO with subsequent tests
library(lmerTest) #Changes the lmer function to output Satterthwaite approximated p-values (can also do ANOVA comparing full model (with KO status) and partial model (without)
#Also, test for differences between scramble and KO for each treatment

#Here, I need to change the treatment to be an ordered factor (1:5)
ussing_cond$treatment<-factor(ussing_cond$treatment,levels=c(unique(ussing_cond$treatment)))
ussing_cond$treatment<-as.numeric(ussing_cond$treatment)

#Now, for each KO in turn, calculate shift between KO and scramble for each treatment and then test for effect on this shift as function of incremental treatments
#for(j in 2:length(unique(ussing_cond$KO_status))){
j=2 #for POU2F3
j=3 #for FOXI1		
		#Get estimated values and standard error, accounting for different donors
		#Treat donor as random effect, since the three inserts are non-independent
		#KO test
		#1
		currKO_1<-ussing_cond[which(ussing_cond$treatment == unique(ussing_cond$treatment)[1] & ussing_cond$KO_status == unique(ussing_cond$KO_status)[j]),]
		currSCBL<-ussing_cond[which(ussing_cond$treatment == unique(ussing_cond$treatment)[1] & ussing_cond$KO_status == unique(ussing_cond$KO_status)[1]),]
		currKO_1$value<-sapply(currKO_1$value,function(x)(x * -1) + mean(currSCBL$value))
		#2
		currKO_2<-ussing_cond[which(ussing_cond$treatment == unique(ussing_cond$treatment)[2] & ussing_cond$KO_status == unique(ussing_cond$KO_status)[j]),]
		currSCBL<-ussing_cond[which(ussing_cond$treatment == unique(ussing_cond$treatment)[2] & ussing_cond$KO_status == unique(ussing_cond$KO_status)[1]),]
		currKO_2$value<-sapply(currKO_2$value,function(x)(x * -1) + mean(currSCBL$value))
		#3
		currKO_3<-ussing_cond[which(ussing_cond$treatment == unique(ussing_cond$treatment)[3] & ussing_cond$KO_status == unique(ussing_cond$KO_status)[j]),]
		currSCBL<-ussing_cond[which(ussing_cond$treatment == unique(ussing_cond$treatment)[3] & ussing_cond$KO_status == unique(ussing_cond$KO_status)[1]),]
		currKO_3$value<-sapply(currKO_3$value,function(x)(x * -1) + mean(currSCBL$value))
		#4
		currKO_4<-ussing_cond[which(ussing_cond$treatment == unique(ussing_cond$treatment)[4] & ussing_cond$KO_status == unique(ussing_cond$KO_status)[j]),]
		currSCBL<-ussing_cond[which(ussing_cond$treatment == unique(ussing_cond$treatment)[4] & ussing_cond$KO_status == unique(ussing_cond$KO_status)[1]),]
		currKO_4$value<-sapply(currKO_4$value,function(x)(x * -1) + mean(currSCBL$value))
		#5
		currKO_5<-ussing_cond[which(ussing_cond$treatment == unique(ussing_cond$treatment)[5] & ussing_cond$KO_status == unique(ussing_cond$KO_status)[j]),]
		currSCBL<-ussing_cond[which(ussing_cond$treatment == unique(ussing_cond$treatment)[5] & ussing_cond$KO_status == unique(ussing_cond$KO_status)[1]),]
		currKO_5$value<-sapply(currKO_5$value,function(x)(x * -1) + mean(currSCBL$value))
		
		#Combind all these together
		currKO<-rbind(currKO_1,currKO_2,currKO_3,currKO_4,currKO_5)
		
		#Now test whether shifts in values between KO and scramble decrease as a function of subsequent tests
		currlmer<-lmer(currKO$value~currKO$treatment + (1 | currKO$donor))
		estimate_diffSCRB<-summary(currlmer)$coefficients[2,1]
		pvalue_diffSCRB<-coef(summary(currlmer))[2,5]
		
		#Now add these values to master table
		currNorm<-data.frame("KO_status"=currKO$KO_status[2],"treatment"=currKO$treatment[1],"estimate_diffSCRB"=estimate_diffSCRB,"pvalue_diffSCRB"=pvalue_diffSCRB)
		ussing_cond_test<-rbind(ussing_cond_test,currNorm)
	}		
}

