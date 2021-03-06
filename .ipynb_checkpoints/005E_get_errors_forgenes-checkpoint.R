date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

input_data_dir = "/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/INPUT_DATA/"

#Load in pval and qval dt
window_pval_dt = fread(pff("/data/005C_pvaldt_100KBerrorwindows.csv"))
window_qval_dt = fread(pff("/data/005C_qvaldt_100KBerrorwindows.csv"))

windows_gr = GRanges(window_pval_dt$chr, 
					 IRanges(window_pval_dt$start, 
							 window_pval_dt$start + 99999))

#Only keep Chromosomes 1-22,X
chr_to_keep = paste("chr", c(1:22, "X"), sep = "")

#Load in gene datable GENCODE
GENCODE = fread(paste0(input_data_dir, "GENCODE_hg38_PROCESSED.txt"))[gene_type == "protein_coding"][chr %in% chr_to_keep]
GENCODE_gr = GRanges(GENCODE$chr, IRanges(GENCODE$start, GENCODE$end))

get_gene_pvals = function(gene_index, p_val_dt, cohorts_to_keep){
        
	gene_name = GENCODE[gene_index]$gene_name

	gene_coordinates = GENCODE_gr[gene_index]
	
	pvals = p_val_dt[, .SD, .SDcols = cohorts_to_keep]
	
	#Get overlapping windows
	overlapping_windows = subjectHits(findOverlaps(gene_coordinates, windows_gr))

	if(length(overlapping_windows) == 0){
		gene_pvals = rep(NA, length(cohorts_to_keep))
	}

	if(length(overlapping_windows) == 1){
		gene_pvals = as.numeric(pvals[overlapping_windows])
	}

	if(length(overlapping_windows) > 1){

		#Get window with most overlap
		widths = unlist(lapply(overlapping_windows, 
							   function(x) width(intersect(gene_coordinates, 
														   windows_gr[x]))))

		top_window = overlapping_windows[which.max(widths)]

		gene_pvals = as.numeric(pvals[top_window])                     
	}

    
    
return(gene_pvals)}

								   
#Get pval for each gene
#Keep only core 14 cancer types
cancer_types_to_keep = c("Breast-AdenoCa", "Kidney-RCC", "Liver-HCC", "ColoRect-AdenoCA", "Eso-AdenoCa", 
						"CNS-GBM", "Stomach-AdenoCA", "Biliary-AdenoCA", "Lung-AdenoCA", "Lung-SCC", 
						 "Head-SCC", "Lymph-BNHL", "Lymph-CLL", "Skin-Melanoma")
							   
							   
gene_pval_dt = as.data.table(do.call("rbind.data.frame", 
									 mclapply(1:nrow(GENCODE), 
											  get_gene_pvals, 
											  window_pval_dt, 
											  cancer_types_to_keep, mc.cores = 8)))
colnames(gene_pval_dt) = cancer_types_to_keep                             
    
gene_pval_dt_full = as.data.table(cbind.data.frame(gene_name = GENCODE$gene_name, gene_pval_dt))
                                 
fwrite(gene_pval_dt_full, pff("/data/005E_gene_pvaldt_100KBerrorwindows.csv"))    
								   
								   
#Get qval for each gene								   
gene_qval_dt = as.data.table(do.call("cbind.data.frame", 
									 lapply(1:26, get_gene_pvals, window_qval_dt)))
colnames(gene_qval_dt) = colnames(gene_qval_dt)[-c(1, 2)]                              
    
gene_qval_dt_full = as.data.table(cbind.data.frame(gene_name = GENCODE$gene_name, gene_qval_dt))
                                 
fwrite(gene_qval_dt_full, pff("/data/005E_gene_qvaldt_100KBerrorwindows.csv"))                                 								   