date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

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

							 
#Get plot data							 
get_plot_data_ct = function(cohort_index, top_N_predictors){
	
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
									  
									  
	get_plot_data_sig = function(signature){
		
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

		#Create data table for plotting                         
		results.m = as.data.table(cbind.data.frame(value = top_predictor_means,
											   SD = top_predictor_SDs, 
											   variable = top_predictors,
											   significance = rep("*", length(top_predictors)), 
											   fill = top_predictor_fill))

		results.m$variable = factor(results.m$variable, levels = top_predictors)
		results.m$fill = factor(results.m$fill, levels = c("Primary tumour CA from matching tissue", 
														   "Normal cells CA from matching tissue", 
														   "Cancer cell line CA from matching tissue",
														   "Replication timing from matching tissue",
														   "Primary tumour CA from non-matching tissue", 
														   "Normal cells CA from non-matching tissue", 
														   "Cancer cell line CA from non-matching tissue",
														   "Replication timing from non-matching tissue")) 									
											
							

										   
		p = ggplot(results.m, aes(x = variable, 
									   y = value,
									   label = significance, 
									   fill = fill))+
			geom_bar(stat = "identity", color="black")+
			geom_text(aes(y=value+SD+SD/8), colour = "black", size = 6)+
			geom_errorbar(aes(ymin = value - SD, ymax = value + SD), width=.2,
				 position=position_dodge(.9))+
			theme_bw()+
			theme(legend.text = element_text(size = 6, colour = "black"),
				legend.title = element_text(size = 6, face = "bold"),
				legend.key.size = unit(0.4, "cm"))+
			scale_fill_manual(drop = F, values = fill_colours)+
			labs(x = "Predictor", y = "IncMSE", fill = "Predictor category",
				title = paste(cohort_name, signature))      		
			}

										   
	lapply(signatures, get_plot_data_sig)}

#Set fill colors								 
fill_colours = c("#FFC300", "#58D68D", "#2980B9", "#FF3352",
				 "#9A7d0A", "#186A3B", "#6C3483", "#922B21")
pdf(pff("data/004H_mutsig_predimp_braplots.pdf"), width = 13, height = 6)
lapply(1:26, get_plot_data_ct, 5)                                         
dev.off()										   
										   
							 