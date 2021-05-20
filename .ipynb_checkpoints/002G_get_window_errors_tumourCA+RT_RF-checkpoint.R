date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

source(pff("/bin/999_run_randomforest_experiment.R"))

#Load 1 MB mutation count dataset
mutation_counts_dt = fread(pff("data/001D_PCAWG_mutation_counts_1MBwindow_processed.csv"))[ , -c(1, 2)] 

Output = mutation_counts_dt[[cohort_index]]                                    
                   
project_code = colnames(mutation_counts_dt)[cohort_index]

print(project_code)

#Load in predictors
ATAC_seq = fread(pff("/data/001B_TCGA_ATACSeq_1MBwindow_processed.csv"))
Vierstra = fread(pff("/data/001E_vierstra_dnaseseq_1MB.csv"))[, -c(1,2)]
RT = fread(pff("/data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))[, -c(1, 2)]

#Combine primary tumour predictor datasets
Preds = as.data.table(cbind.data.frame(ATAC_seq[, -c(1,2)], 
									   Vierstra[, .SD, .SDcols = c("h.renal.cell.carcinoma-DS26693", 
																	 "h.renal.cell.carcinoma-DS37973")],
									   RT))

#Get observed vs predicted data tables
get_dt = function(cohort_index){
    cohort_name = colnames(mutation_counts_dt)[cohort_index]
    
    output = mutation_counts_dt[[cohort_index]]
    
    #Remove hypermutated windows from lymphoma and leukemia
    if(cohort_name %in% c("Lymph-CLL", "Lymph-BNHL")){
        ATAC_seq = ATAC_seq[-c(292,2438)]
        Preds = Preds[-c(292,2438)]
        output = output[-c(292,2438)]
    }
    
    rf=randomForest(x = Preds, 
					y = output, 
					keep.forest = T, 
					ntree = 1000, 
					do.trace = F, 
					importance = F)
    
    #Save observed vs. predicted values
    dt = as.data.table(cbind.data.frame(chr = ATAC_seq$chr, 
										start = ATAC_seq$start, 
										observed = output, 
										predicted = rf$predicted))
    
    fwrite(dt, paste0(pff("data/002G_tumourCAplusRT_RF_obsvsexpected/"), cohort_name, ".csv"))}

mclapply(1:26, get_dt, mc.cores = 16)