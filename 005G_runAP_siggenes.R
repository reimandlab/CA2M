date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

input_data_dir = "/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/INPUT_DATA/"

gmt = paste0(input_data_dir, "GOBP_REAC.gmt")

#Load in gene pval dt and remove genes with NAs
gene_pval_dt = na.omit(fread(pff("data/005E_gene_pvaldt_100KBerrorwindows.csv")))


cancer_types_to_keep = c("Breast-AdenoCa", "Kidney-RCC", "Liver-HCC", "ColoRect-AdenoCA", "Eso-AdenoCa", 
						"CNS-GBM", "Stomach-AdenoCA", "Biliary-AdenoCA", "Lung-AdenoCA", "Lung-SCC", 
						 "Head-SCC", "Lymph-BNHL", "Lymph-CLL", "Skin-Melanoma")	
gene_pval_dt_core = gene_pval_dt[, .SD, .SDcols = cancer_types_to_keep]

gene_pval_dt_core_matrix = data.matrix(gene_pval_dt_core)
rownames(gene_pval_dt_core_matrix) =  gene_pval_dt[[1]]

AP_integrated = ActivePathways(unique(gene_pval_dt_core_matrix), 
							 gmt, 
							 correction.method="fdr", 
							 cytoscape.file.tag = pff("data/005G_Integrated_AP_cytofiles/"))

fwrite(AP_integrated, pff("data/005G_integratedAP_results.csv"))

nrow(AP_integrated)
# [1] 220

#Get how many cancer types each pathway takes evidence from
num_evidence = sapply(1:nrow(AP_integrated), 
					  function(x) length(unlist(strsplit(unlist(AP_integrated$evidence[x]), 
														 split = ",", fixed = T))))
num_combined = sum(unlist(AP_integrated$evidence) == "combined")
					  
sum(num_evidence > 1) + num_combined				  
# [1] 162
					  
#Edit subroups file to add PCAWG colours
library(mgsub)

subgroups = fread(pff("data/005G_Integrated_AP_cytofiles/subgroups.txt"))

change_cols = function(file, old_cols, new_cols){
    
    file$instruct = unlist(lapply(file$instruct, function(x) mgsub(x, old_cols, new_cols)))
                                
    return(file)                            
    
    
}
PCAWG_colours = readRDS(paste0(input_data_dir, "PCAWG_colour_palette.RDS"))

new_cols = c(as.character(PCAWG_colours)[match(tolower(colnames(subgroups)[-c(1,  
																			  ncol(subgroups))]), 
											   names(PCAWG_colours))])
new_cols[2] = "#FF0000"								  
								  
old_cols = c("#FF0000", 
			 "#FF6D00", 
			 "#FFDB00", 
			 "#B6FF00",  
			 "#49FF00", 
			 "#00FF24", 
			 "#00FF92", 
			 "#00FFFF", 
			 "#0092FF", 
			 "#0024FF",  
			 "#4900FF",  
			 "#B600FF", 
			 "#FF00DB", 
			 "#FF006D")
                                
subgroups_PCAWG = change_cols(subgroups,old_cols,new_cols)
write.table(subgroups_PCAWG, 
			pff("data/005G_Integrated_AP_cytofiles/subgroups_PCAWG.txt"), 
			quote = FALSE, 
			sep = '\t', 
		    row.names = F)                                
