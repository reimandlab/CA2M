#Get cohort number as argument
args = commandArgs(trailingOnly = TRUE)
cohort_index = as.integer(args[1])

date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

source(pff("/bin/999_run_randomforest_experiment.R"))

#Load 1 MB mutation count dataset
mutation_counts_dt = fread(pff("data/001D_PCAWG_mutation_counts_1MBwindow_processed.csv"))[ ,.SD, .SDcols = -c(1, 2)] 

Output = mutation_counts_dt[[cohort_index]]                                    
                   
project_code = colnames(mutation_counts_dt)[cohort_index]

print(project_code)

#Load in predictors
ATAC_seq = fread(pff("/data/001B_TCGA_ATACSeq_1MBwindow_processed.csv"))[,.SD,.SDcols=-c(1,2)]
Vierstra = fread(pff("/data/001E_vierstra_dnaseseq_1MB.csv"))[,.SD,.SDcols=-c(1,2)]
RT = fread(pff("/data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))[,-c(1, 2)]

#Combine primary tumour predictor datasets
Preds = as.data.table(cbind.data.frame(ATAC_seq, 
									   Vierstra[, .SD, .SDcols = c("h.renal.cell.carcinoma-DS26693", 
																	 "h.renal.cell.carcinoma-DS37973")],
									   RT))

#Remove hypermutated windows from lymphoma and leukemia
if(project_code %in% c("Lymph-CLL", "Lymph-BNHL")){
    Preds = Preds[-c(292, 2438)]
    Output = Output[-c(292, 2438)]
}

#Run random forest with monte carlo cross-validation
RF_result = run_RF_and_get_pred_importances(Preds, Output, 1000, n_tree = 1000, cores = 16, train_split = 0.8)

#Save results
fwrite(RF_result, paste0(pff("/data/002D_tumorCAplusRT_RF_Results/"), project_code, ".csv"))

print("Finished")