date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))
library(chemometrics)
error_dt_paths = list.files(pff("/data/005A_100KB_RF_errorDT_results"), full.names = T)
error_dt_cancertypes = unlist(lapply(list.files(pff("/data/005A_100KB_RF_errorDT_results")), 
											  function(x) unlist(strsplit(x, split = ".csv"))[1]))
                           
                                   
get_pvals = function(cohort_index, log){
    print(cohort_index)
    project_code = error_dt_cancertypes[cohort_index]
    
    error_dt = fread(error_dt_paths[cohort_index])
    
    error_dt$errors = as.numeric(error_dt$observed - error_dt$predicted)
	#Trim extreme values for mean and sd
	p_values = pnorm(error_dt$errors, 
					    mean = mean(error_dt$errors),
						sd = sd(error_dt$errors),
						lower.tail = F, 
						log.p = log) #Natural log

    if(cohort_index == 1){
        return(cbind.data.frame(chr = error_dt$chr, 
								start = error_dt$start, 
								p_values))
    }else{
        return(p_values)
    }
}

#Save p-val of errors									 
p_val_dt = as.data.table(do.call("cbind.data.frame", lapply(c(1:26), get_pvals, log = F)))                                   
colnames(p_val_dt)[-c(1, 2)] = error_dt_cancertypes[c(1:26)] 
fwrite(p_val_dt, pff("/data/005C_pvaldt_100KBerrorwindows.csv"))									 

									 
									 
#Save log p-val of errors									 
log_p_val_dt = as.data.table(do.call("cbind.data.frame", lapply(c(1:26), get_pvals, log = T)))                                   
colnames(log_p_val_dt)[-c(1, 2)] = error_dt_cancertypes[c(1:26)]                                   
fwrite(log_p_val_dt, pff("/data/005C_logpvaldt_100KBerrorwindows.csv"))

#Save q-val of errors
q_val_dt = as.data.table(cbind.data.frame(p_val_dt[,c(1, 2)], 
										  apply(p_val_dt[,-c(1, 2)], 
												2, 
												p.adjust, 
												method = "fdr")))         									 
fwrite(q_val_dt, pff("/data/005C_qvaldt_100KBerrorwindows.csv"))

#Get dt for 12 core cancer types
qval_dt_core = q_val_dt[, c(3, 6, 7, 10, 11, 12, 14, 15, 16, 17, 18, 19, 25, 26)]

#Get number of unique windows significant in at least one cancer type
sum(rowSums(qval_dt_core < 0.05, na.rm = T) > 0, na.rm = T)
# 1330

#Get number of windows significant for each cancer type (overlapping)
num_sig = apply(qval_dt_core, 2, function(x) sum(x < 0.05, na.rm = T))
num_sig				
#  Biliary-AdenoCA   Breast-AdenoCa          CNS-GBM ColoRect-AdenoCA 
#              110               62              217              114 
#      Eso-AdenoCa         Head-SCC       Kidney-RCC        Liver-HCC 
#              166              108               26               48 
#     Lung-AdenoCA         Lung-SCC       Lymph-BNHL        Lymph-CLL 
#              155              182               55              110 
#    Skin-Melanoma  Stomach-AdenoCA 
#               95              128 
				
#For each window get number of significant cancer types
num_sig_cancertypes = table(apply(qval_dt_core, 1, function(x) sum(x < 0.05, na.rm = T)))
								  
num_sig_cancertypes								  
#     0     1     2     3     4     5     9    11 
# 23110  1147   149    25     5     1     1     2		