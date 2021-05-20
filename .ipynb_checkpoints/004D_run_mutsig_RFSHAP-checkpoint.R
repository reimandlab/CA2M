#Get cohort number as argument
args = commandArgs(trailingOnly = TRUE)
cohort_index = as.integer(args[1])

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

project_codes = list.files(pff("/data/001I_Sig_RMVs/"))


run_SHAP = function(cohort_index){
project_code = project_codes[cohort_index]

print(project_code)

Sigs = unlist(lapply(list.files(paste0(pff("/data/001I_Sig_RMVs/"), project_code)), 
            function(x) unlist(strsplit(x, split = ".csv"))[1]))
                   
#Function to convert R2 to adjusted R2
adjust_R2 = function(R2, n, k){
	adjusted_R2 = 1 - (1 - R2)*(n - 1)/(n - k - 1)
	return(adjusted_R2)
}
    
dir.create(paste0(pff("data/004D_Sig_RF_SHAPvalues/"), project_code, "/"))
                   
run_sig_RF = function(signature){
    
    #Process MAF into RMV
    Output = as.numeric(fread(paste0(pff("/data/001I_Sig_RMVs/"), 
									 project_code, 
									 "/", 
									 signature, 
									 ".csv"))[["Mut_counts"]])
    
    if(project_code %in% c("Lymph-CLL", "Lymph-BNHL")){
        Preds = Preds[-c(292,2438)]
        Output = Output[-c(292,2438)]}
    
    SHAP_values = run_RF_and_get_SHAP(Preds, Output)

    fwrite(SHAP_values, paste0(pff("data/004D_Sig_RF_SHAPvalues/"), 
							   project_code, "/", signature, ".csv"))

}

mclapply(Sigs, run_sig_RF, mc.cores = 9)

print("Finished")}
					 
cohorts_to_run = c(8, 9, 10, 12, 13, 14, 15, 16, 17, 23, 24)
					 
lapply(cohorts_to_run, run_SHAP)					 
					 
					 