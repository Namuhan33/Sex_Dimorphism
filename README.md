methylation
	
	extract_paired_sample
		# Retrieve all the sample barcodes available on TCGA for each cancer
		# kept the barcodes for samples in pair (used in DM_analysis)
	
	DM_analysis
		# DM analysis on all/female/male group
	
	cgID_convert
		# based on the ID of each DM CpG, to match all the potential targets for each DM CpG
	
	sex-specific_DMC
		# overlap the female-based DM CpGs and male-based DM CpGs to obtain female-specific and male-specific DM CpGs
	
	CG_classification_prep
		# find the nearest gene of each CpG site and find the longest transcript of the nearest gene, generate an annotation file;
		# generating the annotation file contain the region of the CpG sitting in, inside or out side a CpG island
	
	DM_target
		# retrieve the potential targets of each DM CpG based on the annotation file generated previously
		# classify the DMC with linked DE targets by their position, methylation status, inside or outside island




miRNA
	
	TCGA_file_preprocessing
		# filter the paired samples and their raw counts, generate expression matrix
	
	DE_analysis
		# Differentially expression analysis on miRNA counts
		# DV analysis

	DEmiR_target
		# retrieve the potential targets of all the DEmiRs

	sex-specific_DMC
		# overlap the female-based and male-based DM CpG lists and identify female-specific and male-specific DM CpGs
		# overlap the female-based and male-based DM CpG targets and identify female-specific and male-specific DM CpG targets




miRNA-methylation

	ref_prepare.R
		# prepare the annotation file contain miRNA-CpG pairs
	
	dm_demir.R
		# link all potential miRNA (as target) to DM CpG and do starburst plot




mRNA

	TCGA_data_preprocessing
		# retrieve all the samples counts and merge them into one counts matrix, keep the pair samples

	DE
		# differential expression analysis under three linear regression model
		# differential expression analysis in female-only samples and male-only samples

	DEcomparison
		# comparison on Design 1 & Design 2, Design 1 & Design 3
		# sex-specific DE genes

	DV
		# DV analysis on all/female/male group, including Levene\'92s test and MAD calculation

	DE_DV
		# overlap between DE gene and DV gene and Fisher exact test
	
	enrichment_analysis
		# sex-specific DE genes enrichment analysis against Hallmark pathway, MsigDB (gene ontology)




mRNA-methylation

	starburst.R
		# compare the potential target list of the DM CpGs to the DE genes and generate starburst plot
	
	true target.R
		# select out the up-regulated target gene & hypo-methylated CpG and down-regulated target gene & hyper-methylated DM CpG pairs to identify true targets





mRNA-miRNA
	
	mRNA-miRNA_target_sex-specific.R
		# identify the targets of sex-specific DE miRNA

	starburst_mRNA-miRNA.R
		# do starburst plot to reveal the true targets of the DE miRNAs

	upset_plot.R
		# identify the common gene/miRNA/CpG/gene sets enriched pathway present in most cancers

	intersect.R
		# identify universal gene/miRNA/CpG/gene sets enriched pathway in 6 cancers




Single-cell
	
	celltype-sex-tissue.R
		# the statistics of normal/tumor cell numbers in female/male group for each cell type

	seurat_colorectal.R
		# filter, normalize, DE analysis, DV analysis on cells
		# DE analysis based on each cell type and DE analysis based on all cells

	compare_to_bulk.R
		# find overlap between DE genes of female epithelial cell and DE genes of female bulk RNA in colon cancer; 
		#  find overlap between DE genes of male epithelial cell and DE genes of male bulk RNA in colon cancer; 
		# find overlap between DE genes of all cells of female and DE genes of female bulk RNA in colon cancer;
		# find overlap between DE genes of all cells of male and DE genes of male bulk RNA in colon cancer




		
