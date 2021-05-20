date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

#Load in predictors
ATAC_seq = fread(pff("/data/001B_TCGA_ATACSeq_1MBwindow_processed.csv"))[,-c(1, 2)]
Roadmap = fread(pff("/data/001A_Roadmap_DNaseSeq_1MBwindow_processed.csv"))[,-c(1, 2)]
Vierstra = fread(pff("/data/001E_vierstra_dnaseseq_1MB.csv"))[,-c(1, 2)]
RT = fread(pff("/data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))[,-c(1, 2)]

Preds = as.data.table(cbind.data.frame(ATAC_seq, Roadmap_DNase, Vierstra, RT))

#Load 1 MB mutation count dataset
mutation_counts_dt = fread(pff("data/001D_PCAWG_mutation_counts_1MBwindow_processed.csv"))[,-c(1, 2)]

get_importances = function(cohort_index){

    project_code = colnames(mutation_counts_dt)[cohort_index] 

    output = mutation_counts_dt[[cohort_index]]

    #Remove hypermutated windows from lymphoma and leukemia
    if(project_code %in% c("Lymph-CLL", "Lymph-BNHL")){
        Preds = Preds[-c(292, 2438)]
        output = output[-c(292, 2438)]
    }
    
    #Train randomForest
    rf = randomForest(x = Preds, 
					  y = output, 
					  keep.forest = T,
					  ntree = 1000, 
					  do.trace = F,
					  importance = T)
    
    #Save predictor importances
    importances=as.numeric(importance(rf, type = 1, scale = F))
    
    return(importances)
}


RF_results = as.data.table(do.call("cbind.data.frame", 
								   mclapply(seq(26), 
											get_importances, 
											mc.cores = 16)))

RF_importances = as.data.table(cbind.data.frame(colnames(Preds), RF_results))
colnames(RF_importances) = c("Predictor", colnames(mutation_counts_dt))
fwrite(RF_importances, pff("/data/003A_feature_importances_allpreds.csv"))

#Get top 5 features
apply(RF_importances[,-1], 2, function(x) RF_importances[,1][order(x, decreasing = T)[1:5]] )


                                                                 
