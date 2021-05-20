#Get cohort number as argument
args = commandArgs(trailingOnly = TRUE)
cohort_index = as.integer(args[1])

date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

project_codes = list.files(pff("/data/001I_Sig_RMVs/"))

project_code = project_codes[cohort_index]

print(project_code)

#Load in predictors
ATAC_seq = fread(pff("/data/001B_TCGA_ATACSeq_1MBwindow_processed.csv"))[, -c(1,2)]
Roadmap = fread(pff("/data/001A_Roadmap_DNaseSeq_1MBwindow_processed.csv"))[,-c(1, 2)]
Vierstra = fread(pff("/data/001E_vierstra_dnaseseq_1MB.csv"))[, -c(1,2)]
RT = fread(pff("/data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))[, -c(1, 2)]

Preds = as.data.table(cbind.data.frame(ATAC_seq, Roadmap, Vierstra, RT))

Sigs = unlist(lapply(list.files(paste0(pff("/data/001I_Sig_RMVs/"), project_code)), 
            function(x) unlist(strsplit(x, split = ".csv"))[1]))
                   
#Function to convert R2 to adjusted R2
adjust_R2 = function(R2, n, k){
	adjusted_R2 = 1 - (1 - R2)*(n - 1)/(n - k - 1)
	return(adjusted_R2)
}

                   
#Create function to run RF
run_RF = function(signature){

    #Process MAF into RMV
    Output = as.numeric(fread(paste0(pff("/data/001I_Sig_RMVs/"), 
									 project_code, 
									 "/", 
									 signature, 
									 ".csv"))[["Mut_counts"]])
    
    if(project_code %in% c("Lymph-CLL", "Lymph-BNHL")){
        Preds = Preds[-c(292,2438)]
        Output = Output[-c(292,2438)]}
    
    #Train randomForest
    rf = randomForest(x = Preds, 
					  y = Output, 
					  keep.forest = T, 
					  ntree = 1000, 
					  do.trace = F, 
					  importance = T)
    
    importances = as.numeric(importance(rf, type = 1, scale = F))
    
    R2 = (cor(rf$predicted, Output))**2
    
    adj_R2 = adjust_R2(R2, nrow(Preds), ncol(Preds))
    
    results = c(adj_R2, importances)
    
    return(results)}

#Run for all mutations

result_dt = as.data.table(do.call("cbind.data.frame", 
								  mclapply(Sigs, run_RF, mc.cores = 8)))
colnames(result_dt) = Sigs
                   
result_dt = as.data.table(cbind.data.frame(Pred = c("Adj_R2", colnames(Preds)), 
										   result_dt))                   
                   
print("Finished")
                   
fwrite(result_dt, 
	   paste0(pff("data/004A_Sig_RF_Results/"), project_code, ".csv"))