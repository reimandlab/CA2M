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
Roadmap = fread(pff("/data/001A_Roadmap_DNaseSeq_1MBwindow_processed.csv"))[,-c(1, 2)]
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

#Permute model output to run permutation test
get_RF_dt = function(seed){
    print(seed)
    set.seed(seed)

    #Convert output to data.table
    output_scrambled = sample(as.numeric(output))

    #Train randomForest
    rf=randomForest(x = Preds, y = output_scrambled, keep.forest = T, 
					ntree = 1000, do.trace = F, importance = T)
    
    return(as.numeric(importance(rf, type = 1, scale = F)))}

blocks = seq(100, 1000, 100)           
           
RF_results = as.data.table(do.call("rbind.data.frame", 
								   mclapply(seq(blocks[block]-99, blocks[block]), 
											get_RF_dt, mc.cores = 10)))
colnames(RF_results) = colnames(Preds)
           
dir.create(paste0(pff("/data/003B_Feature_importance_permutations/"), project_code, "/"))

fwrite(RF_results, paste0(pff("/data/003B_Feature_importance_permutations/"), project_code, "/", block, ".csv"))

print("Finished")
                                                                