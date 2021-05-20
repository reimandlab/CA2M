date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

source(pff("/bin/999_run_RF_withSHAP.R"))

#Load in predictors
ATAC_seq = fread(pff("/data/001B_TCGA_ATACSeq_1MBwindow_processed.csv"))[, -c(1,2)]
Roadmap = fread(pff("/data/001A_Roadmap_DNaseSeq_1MBwindow_processed.csv"))[,-c(1, 2)]
Vierstra = fread(pff("/data/001E_vierstra_dnaseseq_1MB.csv"))[, -c(1,2)]
RT = fread(pff("/data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))[, -c(1, 2)]

Preds = as.data.table(cbind.data.frame(ATAC_seq, Roadmap, Vierstra, RT))

#Load 1 MB mutation count dataset
mutation_counts_dt = fread(pff("data/001D_PCAWG_mutation_counts_1MBwindow_processed.csv"))[ , -c(1, 2)] 

#Get SHAP values in each window for each predictor in each cancer type
get_SHAP_values = function(cohort_index){

    print(cohort_index)
    
    project_code = colnames(mutation_counts_dt)[cohort_index] 

    output = mutation_counts_dt[[cohort_index]]

    #Remove hypermutated windows from lymphoma and leukemia
    if(project_code %in% c("Lymph-CLL", "Lymph-BNHL")){
        Preds = Preds[-c(292, 2438)]
        output = output[-c(292, 2438)]
    }
    

    #Train randomForest and get SHAP values
    SHAP_values = run_RF_and_get_SHAP(Preds, output)
    
    fwrite(SHAP_values, paste0(pff("data/003F_RF_SHAP_Results/"), project_code, ".csv"))
}


mclapply(seq(26), get_SHAP_values, mc.cores = 13)



                                                                 
