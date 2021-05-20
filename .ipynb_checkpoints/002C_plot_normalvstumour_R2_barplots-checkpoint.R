date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

#Load in cohort names in analysis
cohort_names = unlist(lapply(list.files(pff("data/002A_primarytumour_RF_Results/")), 
										function(x) unlist(strsplit(x, 
																	split = ".csv", 
																	fixed = T))[1]))



#Load in tumor CA results
data_paths_tumor = list.files(pff("data/002A_primarytumour_RF_Results"), full.names = T)
accuracies_tumor = unlist(lapply(data_paths_tumor, 
								 function(x) as.numeric(fread(x)[["Test_adj_R2"]])))


#Load in normal CA results
data_paths_normal = list.files(pff("data/002B_normalcell_RF_Results"), full.names = T)
accuracies_normal = unlist(lapply(data_paths_normal, 
								 function(x) as.numeric(fread(x)[["Test_adj_R2"]])))

#Combine results into data table
accuracies.dt = as.data.table(cbind.data.frame(cohort_name = unlist(lapply(cohort_names,rep, 1000)), 
											   Tumor_adjR2 = accuracies_tumor, 
											   Normal_adjR2 = accuracies_normal))
								  
accuracies.dt$diff_adjR2 = accuracies.dt$Tumor_adjR2 - accuracies.dt$Normal_adjR2
                               

plot_dt = accuracies.dt[,.(median_tumor_adjR2 = median(Tumor_adjR2), 
						   median_normal_adjR2 = median(Normal_adjR2), 
						   mean_diff_adjR2 = mean(diff_adjR2), 
						   low_CI_adjR2 = quantile(diff_adjR2, 0.025), 
						   high_CI_adjR2 = quantile(diff_adjR2, 0.975)), 
						by = cohort_name]                             
                
order = plot_dt$cohort_name[order(plot_dt$mean_diff_adjR2)]

#Get empirical p-val
get_emp_pval = function(cohort_name_in){
	data = accuracies.dt[cohort_name == cohort_name_in]
	mean_diff = mean(data$diff_adjR2)
	p_val  = sum(sign(data$diff_adjR2) != sign(mean_diff))/1000
	return(p_val)
	
}
emp_pvals = sapply(plot_dt$cohort_name, get_emp_pval)			   
			   
#Save results
plot_dt$emp_pvals = emp_pvals
plot_dt$emp_FDR = p.adjust(emp_pvals, method="fdr")							  

plot_dt$bar_label = ifelse(plot_dt$emp_pvals < 0.001, "***",
						ifelse(plot_dt$emp_pvals < 0.01, "**",
							ifelse(plot_dt$emp_pvals < 0.05, "*",
								ifelse(plot_dt$emp_pvals < 0.1, "\u2022", ""))))		
							  
cairo_pdf(pff("data/002C_CA_tumourvsnormal_barplots.pdf"), width = 9, height = 5)                       
ggplot(plot_dt, aes(x = factor(cohort_name, 
					levels = rev(order)), 
					y = mean_diff_adjR2, 
					fill = median_tumor_adjR2)) + 
  geom_bar(stat = "identity", colour = 'black')+
  geom_errorbar(aes(ymin = low_CI_adjR2, ymax = high_CI_adjR2), width = .2)+
  geom_text(aes(label = bar_label, y = high_CI_adjR2 + 0.02), fontface = "bold", size = 3)+
    theme_bw()+
    scale_fill_gradient2(low = "white", mid = "white", high = "red", midpoint = 0.3)+
    labs(y = "Change in model accuracy,\ntumour vs. normal tissue chromatin (adjusted R2)", 
		 x = "Cohort name", 
		 fill = "Model accuracy using\ntumour chromatin\n(adjusted R2)")+
    theme(axis.title = element_text(size=9),
        axis.text.x = element_text(size = 8, 
								   angle = 90, 
								   hjust = 1, 
								   vjust = 0.5,
								   colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
         legend.text = element_text(size = 8),
         legend.title =element_text(size = 10))
dev.off()                               
							  
							  
#Create CSV of statistical results
pvals = sapply(unique(accuracies.dt$cohort_name), 
			   function(x) wilcox.test(accuracies.dt[cohort_name==x]$diff_adjR2)$p.val)	
plot_dt$p_val = as.numeric(pvals[match(plot_dt$cohort_name, names(pvals))])
plot_dt$FDR = p.adjust(plot_dt$p_val, method="fdr")
			   
			   

			   
fwrite(accuracies.dt, pff("data/002C_resultsdata.csv"))
fwrite(plot_dt, pff("data/002C_results_summary.csv"))			   