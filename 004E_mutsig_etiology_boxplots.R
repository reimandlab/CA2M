date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))



filepaths = list.files(pff("data/004A_Sig_RF_Results"), 
					   full.names = T)
cohort_names = unlist(lapply(list.files(pff("data/004A_Sig_RF_Results")),
                            function(x) unlist(strsplit(x ,split = ".csv"))[1]))

PCAWG_mutations_dt_hg38_SNP = fread(pff("data/001G_PCAWG_MAF_SNP_withSBSSIG.csv"))
PCAWG_sample_table = PCAWG_mutations_dt_hg38_SNP[, 
												 .(Cohort_size = uniqueN(Donor_ID)), 
												 by = Project_Code]

n_mut_sig_project_SBS = PCAWG_mutations_dt_hg38_SNP[, .N, by = .(Project_Code, SBS_Signature)]
n_mut_sig_project = PCAWG_mutations_dt_hg38_SNP[, .N, by = .(Project_Code)]
n_samples = PCAWG_mutations_dt_hg38_SNP[, length(unique(Donor_ID)), by = Project_Code]
						 
							 
get_plot_data = function(cohort_index){
    
    cohort_name = cohort_names[cohort_index]
	print(cohort_name)

    data = fread(filepaths[which(cohort_names == cohort_name)])
                            
    signatures = colnames(data)[-c(1, 2)]
    
    Adj_R2 = as.numeric(data[1, -c(1, 2)])
	
    get_mut_fraction = function(SBS){
        if(cohort_name == "PANCAN"){
                n_mut = sum(n_mut_sig_project_SBS[SBS_Signature == SBS]$N)

                n_fraction = n_mut/nrow(PCAWG_mutations_dt_hg38_SNP)
            }
            else{
				mut_nums = as.numeric(fread(paste0(pff("/data/001I_Sig_RMVs/"), 
									 cohort_name, 
									 "/", 
									 SBS, 
									 ".csv"))[["Mut_counts"]])
				
				
                n_fraction = mean(mut_nums)/n_samples[Project_Code == cohort_name]$V1
				
            }
            return(n_fraction)}
    mut_fractions = unlist(lapply(signatures, get_mut_fraction))
    
    get_mut_num = function(SBS){
        if(cohort_name=="PANCAN"){
            n_mut = sum(n_mut_sig_project_SBS[SBS_Signature == SBS]$N)
            }
        else{
            n_mut = n_mut_sig_project_SBS[Project_Code == cohort_name & SBS_Signature == SBS]$N
            }
            return(n_mut)}
    mut_nums = unlist(lapply(signatures, get_mut_num))
	
    plot_dt = cbind.data.frame(cohort_name = rep(cohort_name, length(signatures)), 
							   signatures, 
							   mut_nums, 
							   Adj_R2, 
							   mut_fractions)
    
    return(plot_dt)
                            
                            
}
full_plot_dt = as.data.table(do.call("rbind.data.frame", 
									 lapply(1:24, get_plot_data)))  
                            
# etiology_dt = as.data.table(cbind.data.frame(
#     etiology = c(rep("APOBEC enzyme activity", 2), 
# 				 rep("Defective DNA repair", 4), 
# 				 rep("Age-related", 1), 
# 				 rep("Exogenous/Carcinogen", 5), 
# 				 rep("Unknown/Other", 10)), 
#     signature = c(c("SBS13", "SBS2"), 
# 				  c("SBS3", "SBS26", "SBS36", "SBS44"), 
# 				  c("SBS1"), 
# 				  c("SBS4", "SBS7a", "SBS7b", "SBS29", "SBS22"),  
# 				  c("SBS8", "SBS9", "SBS12", "SBS17a", "SBS17b", "SBS18", "SBS40", "SBS41", "SBS60","SBS5"))))
							 
							 
etiology_dt = as.data.table(cbind.data.frame(
    etiology = c(rep("APOBEC enzyme activity", 2), 
				 rep("Defective DNA repair", 4), 
				 rep("SBS1", 1), 
				 rep("Exogenous/Carcinogen", 5),
				 rep("SBS5, SBS40", 2),
				 rep("Unknown/Other", 8)), 
    signature = c(c("SBS13", "SBS2"), 
				  c("SBS3", "SBS26", "SBS36", "SBS44"), 
				  c("SBS1"), 
				  c("SBS4", "SBS7a", "SBS7b", "SBS29", "SBS22"),
				  c("SBS5", "SBS40"),
				  c("SBS8", "SBS9", "SBS12", "SBS17a", "SBS17b", "SBS18", "SBS41", "SBS60"))))							 
                            
full_plot_dt$etiologies = etiology_dt$etiology[match(full_plot_dt$signatures, 
											  etiology_dt$signature)]
                            
#Create boxplots
full_plot_dt$etiologies = factor(full_plot_dt$etiologies, 
								 levels = c("SBS1", 
											"APOBEC enzyme activity", 
											"Defective DNA repair", 
											"Exogenous/Carcinogen",
											"SBS5, SBS40",
											"Unknown/Other"))

full_plot_dt$labels = paste0(full_plot_dt$cohort_name, 
							"\n", 
							full_plot_dt$signatures)                        


full_plot_dt$labels_to_keep = ifelse(full_plot_dt$labels %in% labels_to_keep, 
									 full_plot_dt$labels, 
									 "")                      
#Plot boxplot
cancer_types_to_keep = c("Breast-AdenoCa", "Kidney-RCC", "Liver-HCC", "ColoRect-AdenoCA", "Eso-AdenoCa", 
						"CNS-GBM", "Stomach-AdenoCA", "Biliary-AdenoCA", "Lung-AdenoCA", "Lung-SCC", 
						 "Head-SCC", "Lymph-BNHL", "Lymph-CLL", "Skin-Melanoma")
							 
full_plot_dt_select = full_plot_dt[cohort_name %in% cancer_types_to_keep]							 

endogenous_etiologies = c("APOBEC enzyme activity", "Defective DNA repair", "SBS1")							 
p_val_1 = wilcox.test(full_plot_dt_select[etiologies %in% endogenous_etiologies]$Adj_R2, 
					  full_plot_dt_select[etiologies == "Exogenous/Carcinogen"]$Adj_R2)$p.val
p_val_2 = wilcox.test(full_plot_dt_select[etiologies %in% endogenous_etiologies]$Adj_R2, 
					  full_plot_dt_select[etiologies == "SBS5, SBS40"]$Adj_R2)$p.val	
p_val_3 = wilcox.test(full_plot_dt_select[etiologies %in% endogenous_etiologies]$Adj_R2, 
					  full_plot_dt_select[etiologies == "Unknown/Other"]$Adj_R2)$p.val	
							 
							 
							 
#Run ANOVA tests to account for mutation burden of signatures
full_plot_dt_select$class = factor(ifelse(full_plot_dt_select$etiologies %in% endogenous_etiologies, 
								   "endogenous", as.character(full_plot_dt_select$etiologies)))

#Run for carcinogenic vs. endogenous signatures							 
H0 = lm(Adj_R2 ~ mut_fractions, 
		data = full_plot_dt_select[class %in% c("endogenous", "Exogenous/Carcinogen")])
							 
H1 = lm(Adj_R2 ~ mut_fractions + class, 
		data = full_plot_dt_select[class %in% c("endogenous", "Exogenous/Carcinogen")])	
							 
pval_anova_1 = anova(H0, H1)$`Pr(>F)`[2]							 
							 							 
#Run for SBS5/40 vs. endogenous signatures							 
H0 = lm(Adj_R2 ~ mut_fractions, 
		data = full_plot_dt_select[class %in% c("endogenous", "SBS5, SBS40")])
							 
H1 = lm(Adj_R2 ~ mut_fractions + class, 
		data = full_plot_dt_select[class %in% c("endogenous", "SBS5, SBS40")])	
							 
pval_anova_2 = anova(H0, H1)$`Pr(>F)`[2]							 
							 
#Run for unknown vs. endogenous signatures							 
H0 = lm(Adj_R2 ~ mut_fractions, 
		data = full_plot_dt_select[class %in% c("endogenous", "Unknown/Other")])
							 
H1 = lm(Adj_R2 ~ mut_fractions + class, 
		data = full_plot_dt_select[class %in% c("endogenous", "Unknown/Other")])	
							 
pval_anova_3 = anova(H0, H1)$`Pr(>F)`[2]									 
							 
#Keep certain labels
labels_to_keep = c("Skin-Melanoma\nSBS7a", 
				   "Lung-AdenoCA\nSBS4", 
				   "Lung-AdenoCA\nSBS40", 
				   "Breast-AdenoCa\nSBS5", 
				   "Breast-AdenoCa\nSBS1", 
				   "Breast-AdenoCa\nSBS13", 
				   "Breast-AdenoCa\nSBS3", 
				   "Eso-AdenoCa\nSBS17b", 
				   "CNS-GBM\nSBS40", 
				   "CNS-GBM\nSBS1", 
				   "Kidney-RCC\nSBS1", 
				   "Kidney-RCC\nSBS5", 
				   "Kidney-RCC\nSBS22",
				   "ColoRect-AdenoCA\nSBS40",
				   "Lung-SCC\nSBS13",
				   "Eso-AdenoCa\nSBS1",
				   "Stomach-AdenoCA\nSBS3", 
				   "Liver-HCC\nSBS4",
				   "Lung-SCC\nSBS4",
				   "Lymph-BNHL\nSBS9") 
							 							 
pos = position_jitter(width = 0.3, seed = 4)                        
pdf(pff("data/004E_etiology_boxplots.pdf"), width = 6, height = 4)                        
ggboxplot(full_plot_dt_select, 
		  x = "etiologies", 
		  y = "Adj_R2", 
		  color = "etiologies", 
		  outlier.shape = NA)+
    geom_jitter(shape = 21, 
				aes(fill = etiologies, 
					size = mut_fractions), 
				colour = "black", 
				position = pos)+
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
         axis.title = element_text(size = 7),
         legend.text = element_text(size = 6),
         legend.title = element_text(size = 7))+              
    labs(x = "Signature etiology", 
		 y = "Adj. R2", 
		 size = "% of total mutations in cohort")+
    scale_colour_d3(drop = F)+
    scale_fill_d3(drop = F)+
    geom_text_repel(aes(label = labels_to_keep), 
					size = 2.5, 
					box.padding = 0.5, 
					min.segment.length = 0, 
					position = pos, 
					segment.size = 0.2)+
    scale_radius(range = c(1, 3))+
    guides(color = FALSE, fill = FALSE)+
	ylim(c(0,1))
dev.off()                        

