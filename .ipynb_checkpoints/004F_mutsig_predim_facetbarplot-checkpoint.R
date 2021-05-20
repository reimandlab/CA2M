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
		results.m = as.data.table(cbind.data.frame(value = rep(1, length(top_predictors)), 
											   variable = top_predictors, 
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
											
							

		#Keep barplot facet order                                 
		results.m$cohort_name = rep(cohort_name, top_N_predictors)
		results.m$signature = rep(signature, top_N_predictors)								   
		results.m$order_infacet = match(results.m$variable, top_predictors)
		i <<- i + top_N_predictors
										   
		return(results.m)								   
		
	}

										   
	result_dt=as.data.table(do.call("rbind.data.frame", lapply(signatures, get_plot_data_sig)))
    return(result_dt)}

i = 0
result_dt = as.data.table(do.call("rbind.data.frame", 
								lapply(1:26, get_plot_data_ct, 5)))                                         
										   
#Plot boxplot
cancer_types_to_keep = c("Breast-AdenoCa", "Kidney-RCC", "Liver-HCC", "ColoRect-AdenoCA", "Eso-AdenoCa", 
						"CNS-GBM", "Stomach-AdenoCA", "Biliary-AdenoCA", "Lung-AdenoCA", "Lung-SCC", 
						 "Head-SCC", "Lymph-BNHL", "Lymph-CLL", "Skin-Melanoma")
										   
										   
plot_dt_select = result_dt[cohort_name %in% cancer_types_to_keep]	
										   
plot_dt_select$cohort_name = factor(plot_dt_select$cohort_name, 
									levels = cancer_types_to_keep)

SBS_order = c("ALL", "SBS1", "SBS2", "SBS13", "SBS3", 
			  "SBS4", "SBS7a", "SBS7b", "SBS22", 
			  "SBS5, SBS40", 
			  "SBS9", "SBS12", "SBS17a", "SBS17b", "SBS18")										   
plot_dt_select$signature = factor(plot_dt_select$signature,
								 levels = SBS_order)
										   
#Set fill colors								 
fill_colours = c("#FFC300", "#58D68D", "#2980B9", "#FF3352",
				 "#9A7d0A", "#186A3B", "#6C3483", "#922B21")
										   
#Plot barplot
p = ggplot(plot_dt_select, aes(x = order_infacet, 
							   y = value, 
							   fill = fill))+
    geom_bar(stat = "identity", color="black")+
    facet_grid(vars(signature), vars(cohort_name), scales = "free")+                                     
    theme_bw()+
    theme(axis.text = element_blank(),
		  axis.ticks = element_blank(),
         strip.text.y.right = element_text(angle = 0, size = 6),
         strip.text.x.top = element_text(size = 5),
        legend.text = element_text(size = 6, colour = "black"),
        legend.title = element_text(size = 6, face = "bold"),
        legend.key.size = unit(0.4, "cm"))+
	scale_fill_manual(drop = F, values = fill_colours)+
    labs(x = "", y = "", fill = "Predictor category")                                     

pdf(pff("data/004F_mutsig_predimp_facetbarplot.pdf"), width = 13, height = 6)
p
dev.off()										   
										   
							 