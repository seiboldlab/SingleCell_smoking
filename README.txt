#The scripts within this directory contain code used to produce the figures within the paper "Dissecting the cellular specificity of smoking effects and reconstructing lineages in the human airway epithelium". They are organized into three major directories: 1) preliminary analyses, 2) main population analyses, and 3) population subset analyses. All scripts contain R code. The necessary packages to load can be found at the top of each script. We used R version 3.5.1 for all analyses. As for the two most important R packages we used for analysis, we used Seurat 3.0.3.9015 and Monocle 2.8.0. All analysis was done
either in a MAC OS or linux environment.

The scripts are organized within these directories as follows:

I. 1_PreliminaryAnalyses	
	a. 1_QCingDataset.r - filters cells and genes from the original 10X count matrices for 15 donors 
		from the in vivo scRNA-seq dataset.
	b. 2_IntegrationAndClustering.r - batch integrates the dataset across 15 donors and then infers main 
		clusters (Fig. 1b, S1b).
	
II. 2_MainPopulationAnalyses
	a. 1_MainClusterAnalyses.r - carries out analyses across the main populations, including differential 
		expression, pathway analysis, Dot plot production (Fig. 1c), assessment of cell proportion shifts 
		between smokers and nonsmokers (Fig. 2c, S3d) marker inference, marker heat map production 
		(Fig. S1c), and plotting of SMG versus epithelial genes (Fig. S1d).
	b. 2_AnalysisOfSmokingEffects.r - carries out smoking effect inference for the main populations, including 
		differential expression, assessment of core versus unique responses (Fig. 2ab, S3b), plotting of bulk 
		smoking genes (Fig. S3a), and creating smoking response gene enrichment plots (Fig. 2d, 7h).
	c. 3_AnalysisOfSMGCells.r - makes dot plots of SMG markers (Fig. 4d), does differential expression, pathway 
		analysis, and heat maps of genes that define SMG secretory versus surface secretory populations 
		(Fig. 4ab, S6ab), makes pie chart (Fig 4c) and bar plots (Fig. S6c) for different mucin panels for SMG
		secretory cells, and makes dot plots of SMG smoking response genes (Fig. S6d).
	d. 4_AnalysisOfCiliatedCells.r - makes box plots in Fig. 5e and Fig. 9d.
	e. 5_AnalysisOfRareCells.r - makes marker violin plots (Fig. 6ac), the CFTR feature plot and  
		CFTR cell-type proportion table (Fig. 7b), makes the heat map in Fig. 6b, makes the FOXI1
		violin plots/pie charts in Fig. 6c, the box plots in Fig. S10, the violin plots in both Fig. S10a and
		in Table S3.
	f. 6_FGNet_allSmokingDEGs.r - code for making the functional gene network plots in Fig. 8.
	g. 7_FGNET_core_up_smokingDEGs.r - code for making the functional gene network plot in Fig. S3c.
	h. 8_FGNET_ciliated.r - code for making the functional gene network plots in Fig. 5gh.

III. 3_PopulationSubsetAnalyses
	a. 1_SecretorySubclusterAnalysis
		1. 1_Subset_integration.r - carries out Seurat integration for the mucus secretory and KRT8-high cells.
		2. 2_ClusteringAndDifferentialExpression.r - carries out clustering and differential expression for
			this subset of cells (Fig. S5a).
		3. 3_Monocle_pseudotime_trajectory.r - reconstructs the pseudotime trajectory in Fig. S5b, produces the
			transcription factor plots in Fig. 3b, and then subsets a "mature" population based on pseudotime.
		4. 4_Monocle_pseudotime_heatmap.r - produces the heat map in Fig. 3a.
		5. 5_AnalysisOfMUCACandMUC5Bcells.r - makes pie charts of mucin profiles (Fig. 3c) and bar graphs of smoking 
		effects (S5e), calculates co-expressed genes for MUC5AC and MUC5B (Fig. 3e), makes box plots of how top
		co-expressed genes for these mucins shift due to smoking (Fig. 3f), and produces subcluster UMAPs (Fig. S5cd).
	b. 2_SMGSubclusterAnalysis 
		1. 1_Subset_integration.r - carries out Seurat integration for SMG cells.
		2. 2_ClusteringAndDifferentialExpression.r - subclusters SMG basal cells (Fig. 4e), carries out
			differential expression, and makes the heat map in Fig. S6f.
		3. 3_Monocle_pseudotime_trajectory.R - reconstructs the pseudotime trajectory in Fig. S6g and the
			transcription factor plot in Fig. 4f.
		4. 4_Monocle_pseudotime_heatmap.r - creates the pseudotime gene heat map in Fig. S6g.
	c. 3_CiliatedSubclusterAnalysis
		1_ClusteringAndDifferentialExpression.r - carries out subclustering and differential expression for
			ciliated cells.
	d. 4_RareCellAnalysis
		1_ClusteringAndDifferentialExpression.r - carries out subclustering and differential expression for
			the rare cells.
		2_RareCell_Ussing_Barplots.r - makes the qPCR and Ussing bar plots (Fig. 6ghi, 7aef, S11bde).

In addition to these scripts, the "CommonFunctions.r" script in the main directory contains functions commonly
called throughout, and should thus always be loaded into the environment.

For running the preliminary analyses, one will need to download the raw, donor-specific, count tables (for
the first script) or the processed combined count table (second script) as specified in the scripts. All
data are housed in GEO here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134174. For the second script,
the "Metadata.txt" table (also included in the main directory) will need to be used.

For all the subsequent scripts, one can generally start with the Seurat object ("Processed_invitro_seurat.Rdata") 
we have also housed at GEO. Once downloaded, this file can be loaded into R (as explained in the scripts)
and used throughout.

The Gene_lists directory contains gene lists used throughout the scripts.
