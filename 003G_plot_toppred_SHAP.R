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

#Load in top preds
importance = fread(pff("/data/003A_feature_importances_allpreds.csv"))

#Load in predictor supplementary table
predictor_supp = fread(paste0(input_data_dir, "predictor_supp_table_4.csv"))           

#Load predictor matching list
matching_list = readRDS(paste0(input_data_dir, "Matching_predictors_dictionary.RDS"))

#Load in permutation test paths
p_val_paths = list.files(pff("/data/003B_Feature_importance_permutations/") ,full.name = T)
p_val_cancer_types = unlist(lapply(list.files(pff("/data/003B_Feature_importance_permutations/")), 
											  function(x) unlist(strsplit(x, split = ".csv"))[1]))
                                 
#Load in bootstrap test paths                                 
bootstrap_paths = list.files(pff("/data/003C_Feature_importance_bootstrapping/"), full.name = T)                
bootstrap_cancer_types = unlist(lapply(list.files(pff("/data/003C_Feature_importance_bootstrapping/")), 
												  function(x) unlist(strsplit(x, split = ".csv"))[1]))
      
#Load in SHAP dataset paths                                     
SHAP_paths = list.files(pff("data/003F_RF_SHAP_Results/"), full.name = T)
SHAP_cancer_types = unlist(lapply(list.files(pff("data/003F_RF_SHAP_Results/")), 
								function(x) unlist(strsplit(x, split = ".csv"))[1]))

#Load in predictors
ATAC_seq = fread(pff("/data/001B_TCGA_ATACSeq_1MBwindow_processed.csv"))[, -c(1,2)]
Roadmap = fread(pff("/data/001A_Roadmap_DNaseSeq_1MBwindow_processed.csv"))[,-c(1, 2)]
Vierstra = fread(pff("/data/001E_vierstra_dnaseseq_1MB.csv"))[, -c(1,2)]
RT = fread(pff("/data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))[, -c(1, 2)]

Preds = as.data.table(cbind.data.frame(ATAC_seq, Roadmap, Vierstra, RT))
                         
plot_SHAP_plot = function(cohort_index, top_N_predictors){
   	
	cohort_name = colnames(importance)[-1][cohort_index]
	predictors = importance[[1]]
	
	#Remove hypermutated windows from lymphoma and leukemia
    if(cohort_name %in% c("Lymph-CLL", "Lymph-BNHL")){
        Preds = Preds[-c(292, 2438)]
    }
	
	#Get permutation test p-values for each predictor
    data = importance[[cohort_index + 1]]
    permutation_importances_dt = as.data.table(do.call("rbind.data.frame", 
									 lapply(list.files(p_val_paths[which(p_val_cancer_types == cohort_name)], 
													   full.name = T), fread)))
    predictor_pvals = unlist(lapply(1:length(predictors), 
									function(x)  sum(data[x] < permutation_importances_dt[[x]])/1000))
									
    significance = ifelse(predictor_pvals == 0, "*", "")
    names(significance) = predictors
	top_N_predictors = min(top_N_predictors, sum(significance == "*"))
	
	#Get predictor importance means and sd from bootstrap experiment for significant predictor                     
    Bootstrap_dt = as.data.table(do.call("rbind.data.frame", 
										 lapply(list.files(bootstrap_paths[which(bootstrap_cancer_types == cohort_name)], 
														   full.name = T), fread)))
    sig_predictor_means = unlist(lapply(which(significance == "*"), 
										function(x) mean(Bootstrap_dt[[x]])))                                   
    top_predictors = names(sig_predictor_means[order(sig_predictor_means, 
													 decreasing = T)][1:top_N_predictors])
    top_predictor_means = as.numeric(sig_predictor_means[order(sig_predictor_means, 
															   decreasing = T)][1:top_N_predictors])
    top_predictor_SDs = unlist(lapply(top_predictors, function(x) sd(Bootstrap_dt[[x]])))								

    #Get predictor categories
    top_predictor_categories = paste0(predictor_supp$`General Category`[match(top_predictors, 
																	   predictor_supp$Predictor)], " CA")
	top_predictor_categories = ifelse(top_predictor_categories=="RT CA", 
									  "Replication timing", top_predictor_categories)
									  
    top_predictor_matching = ifelse(top_predictors %in% matching_list[[cohort_name]], 
									"Matching", 
									"Non-Matching")
    top_predictor_fill = unlist(lapply(1:length(top_predictors), 
									   function(x) ifelse(top_predictor_matching[x]=="Matching", 
														  paste(top_predictor_categories[x], "from matching tissue", sep = " "), 
														  paste(top_predictor_categories[x], " from non-matching tissue", sep = ""))))	
									   
	top_predictor_descriptions = predictor_supp$Description[match(top_predictors, predictor_supp$Predictor)]	
    
    #Load in SHAP values
    SHAP_dt = fread(SHAP_paths[which(SHAP_cancer_types == cohort_name)])
                                      
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
    plt$title(cohort_name)
    pdf$savefig(fig, bbox_inches = "tight")                                  
                                      
}
                                      
                                      
plt$close()                           
pdf = backend$PdfPages(pff("data/003G_toppred_SHAP_plots.pdf"))                          
lapply(1:26, plot_SHAP_plot, 5)     
pdf$close()                                                               