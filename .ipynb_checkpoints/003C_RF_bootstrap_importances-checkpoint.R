#Get cohort number as argument
args = commandArgs(trailingOnly = TRUE)
cohort_index = as.integer(args[1])
block = as.integer(args[2])

date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

#Load in predictors
ATAC_seq = fread(pff("/data/001B_TCGA_ATACSeq_1MBwindow_processed.csv"))[,-c(1, 2)]
Roadmap_DNase=fread("/.mounts/labs/reimandlab/private/users/oocsenas/ChrAcc2MutData/Roadmap_DNaseSeq_1MBwindow_processed.csv")[,-c(1, 2)]
Vierstra = fread(pff("/data/001E_vierstra_dnaseseq_1MB.csv"))[,-c(1, 2)]
RT = fread(pff("/data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))[,-c(1, 2)]

Preds = as.data.table(cbind.data.frame(ATAC_seq,Roadmap_DNase,Vierstra, RT))

#Load 1 MB mutation count dataset
mutation_counts_dt = fread(pff("data/001D_PCAWG_mutation_counts_1MBwindow_processed.csv"))[,-c(1, 2)]

project_code = colnames(mutation_counts_dt)[cohort_index] 

output = mutation_counts_dt[[cohort_index]]

#Remove hypermutated windows from lymphoma and leukemia
if(project_code %in% c("Lymph-CLL", "Lymph-BNHL")){
    Preds = Preds[-c(292, 2438)]
    output = output[-c(292, 2438)]
}

#Bootstrap sample dataset and train random forest
get_RF_dt=function(seed){
    print(seed)
    set.seed(seed)

    #Convert input and output to bootsrapped data.table
    bootstrap_sample = sample(nrow(Preds), replace = T)
    input_bootstrap = Preds[bootstrap_sample]
    output_bootstrap = output[bootstrap_sample]

    #Train randomForest
    rf = randomForest(x = input_bootstrap, 
					  y = output_bootstrap, 
					  keep.forest = T, 
					  ntree = 1000, 
					  do.trace = F,
					  importance = T)
    
    return(as.numeric(importance(rf, type = 1, scale = F)))}

blocks=seq(100, 1000, 100)           
           
RF_results=as.data.table(do.call("rbind.data.frame", 
								 mclapply(seq(blocks[block]-99, blocks[block]), 
										  get_RF_dt, mc.cores = 10)))
colnames(RF_results) = colnames(Preds)
           
dir.create(paste0(pff("/data/003C_Feature_importance_bootstrapping/"), project_code, "/"))

fwrite(RF_results, paste0(pff("/data/003C_Feature_importance_bootstrapping/"), project_code, "/", block, ".csv"))

print("Finished")
                                                                   
