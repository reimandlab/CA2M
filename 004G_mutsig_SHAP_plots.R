date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

library(reticulate)
shap = import("shap")
plt = import('matplotlib.pyplot')
backend=import("matplotlib.backends.backend_pdf")

input_data_dir = "/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/INPUT_DATA/"

importance_paths = list.files(pff("data/004A_Sig_RF_Results"), 
					   full.names = T)
cohort_names = unlist(lapply(list.files(pff("data/004A_Sig_RF_Results")),
                            function(x) unlist(strsplit(x ,split = ".csv"))[1]))

#Permutation test paths and cancer types							 
permutation_paths = list.files(pff("data/004B_Sig_importance_permRF_results/"), full.name = T)							 
permutation_cohort_names = list.files(pff("data/004B_Sig_importance_permRF_results/"))							 
							 
#Bootstrap paths and cancer types
bootstrap_paths = list.files(pff("data/004C_Sig_importance_bootstrapRF_results"), full.name = T)           
bootstrap_cohort_names = list.files(pff("data/004C_Sig_importance_bootstrapRF_results"))
							 
#Load in predictor supplementary table
predictor_supp = fread(paste0(input_data_dir, "predictor_supp_table_4.csv"))           

#Load predictor matching list
matching_list = readRDS(paste0(input_data_dir, "Matching_predictors_dictionary.RDS"))

#Get paths of SHAP results
SHAP_paths = list.files(pff("data/004D_Sig_RF_SHAPvalues"), full.name = T)    							 
SHAP_cohort_names = list.files(pff("data/004D_Sig_RF_SHAPvalues"))							 

#Load in predictors
ATAC_seq = fread(pff("/data/001B_TCGA_ATACSeq_1MBwindow_processed.csv"))[, -c(1,2)]
Roadmap = fread(pff("/data/001A_Roadmap_DNaseSeq_1MBwindow_processed.csv"))[,-c(1, 2)]
Vierstra = fread(pff("/data/001E_vierstra_dnaseseq_1MB.csv"))[, -c(1,2)]
RT = fread(pff("/data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))[, -c(1, 2)]

Preds = as.data.table(cbind.data.frame(ATAC_seq, Roadmap, Vierstra, RT))							 
							 
#Get plot data							 
plot_SHAP_ct = function(cohort_index, top_N_predictors){
	
	cohort_name = cohort_names[cohort_index]
    importances = fread(importance_paths[cohort_index])
    predictors = importances[[1]][-1]
    
	
    signatures = unlist(lapply(list.files(paste0(pff("data/001I_Sig_RMVs/"), cohort_name)), 
            function(x) unlist(strsplit(x, split = ".csv"))[1]))
    
	#Get permutation test paths
    permutation_paths_ct = list.files(permutation_paths[which(permutation_cohort_names == cohort_name)], 
									full.name = T)
    permutation_sigs_ct = unlist(lapply(list.files(permutation_paths[which(permutation_cohort_names == cohort_name)]), 
								 function(x) unlist(strsplit(x, split = ".csv"))[1]))
	#Get bootstrap paths								  
	bootstrap_paths_ct = list.files(bootstrap_paths[which(bootstrap_cohort_names == cohort_name)], 
									full.name = T)
    bootstrap_sigs_ct = unlist(lapply(list.files(bootstrap_paths[which(bootstrap_cohort_names == cohort_name)]), 
								 function(x) unlist(strsplit(x, split = ".csv"))[1]))
	
	#Get SHAP paths
	SHAP_paths_ct = list.files(SHAP_paths[which(SHAP_cohort_names == cohort_name)], 
									full.name = T)								  
	SHAP_sigs_ct = unlist(lapply(list.files(SHAP_paths[which(SHAP_cohort_names == cohort_name)]), 
								 function(x) unlist(strsplit(x, split = ".csv"))[1]))
									  
	#Remove hypermutated windows from lymphoma and leukemia
    if(cohort_name %in% c("Lymph-CLL", "Lymph-BNHL")){
        Preds = Preds[-c(292, 2438)]
    }	
								 
	plot_SHAP_sig = function(signature){
		
		#Get permutation test p-values for each predictor
		data = importances[[signature]][-1]
		permutation_importances_dt = fread(permutation_paths_ct[which(permutation_sigs_ct == signature)])
		
		predictor_pvals = unlist(lapply(1:length(predictors), 
										function(x)  sum(data[x] < permutation_importances_dt[[x]])/1000))

		significance = ifelse(predictor_pvals == 0, "*", "")
		names(significance) = predictors
		top_N_predictors = min(top_N_predictors, sum(significance == "*"))

		#Get predictor importance means and sd from bootstrap experiment for significant predictor                     
		Bootstrap_dt = fread(bootstrap_paths_ct[which(bootstrap_sigs_ct == signature)])
										
		sig_predictor_means = unlist(lapply(which(significance == "*"), 
											function(x) mean(Bootstrap_dt[[x]])))                                   
		top_predictors = names(sig_predictor_means[order(sig_predictor_means, 
														 decreasing = T)][1:top_N_predictors])
											
		#Get predictor categories
		top_predictor_descriptions = predictor_supp$Description[match(top_predictors, predictor_supp$Predictor)]
											
		#Load in SHAP values
		SHAP_dt = fread(SHAP_paths_ct[which(SHAP_sigs_ct == signature)])
		
		#Keep SHAP values for top predictors
		SHAP_dt_top = SHAP_dt[,.SD,.SDcols = top_predictors]
		colnames(SHAP_dt_top) = top_predictor_descriptions
		Input_top = Preds[,.SD,.SDcols = top_predictors]                                  
		colnames(Input_top) = top_predictor_descriptions

		#Plot SHAP summary plot                                  
		fig = plt$figure(figsize = c(3, 3))
		shap$summary_plot(data.matrix(SHAP_dt_top), 
						  Input_top, 
						  show = F,
						  sort = F, 
						  max_display = 5L)
		plt$title(paste(cohort_name, signature))
		pdf$savefig(fig, bbox_inches = "tight")}
											
	lapply(signatures, plot_SHAP_sig)}

cohorts_to_keep = c(1, 4, 5, 8, 9, 10, 12, 13, 14, 15, 16, 17, 23, 24)											

plt$close()                           											
pdf = backend$PdfPages(pff("data/004G_mutsig_SHAP_plots.pdf"))                          
lapply(cohorts_to_keep, plot_SHAP_ct, 5)     
pdf$close()   											