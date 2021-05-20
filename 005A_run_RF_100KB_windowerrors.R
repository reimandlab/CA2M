#Get cohort number as argument
args = commandArgs(trailingOnly = TRUE)
cohort_index = as.integer(args[1])

date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

mutation_counts_dt = fread(pff("data/001D_PCAWG_mutation_counts_100KBwindow_processed.csv"))[ , -c(1, 2)] 

#Load in predictors
ATAC_seq = fread(pff("/data/001B_TCGA_ATACSeq_100KBwindow.csv"))
Roadmap = fread(pff("/data/001A_Roadmap_DNaseSeq_100KBwindow_processed.csv"))[, -c(1, 2)]
Vierstra = fread(pff("/data/001E_vierstra_dnaseseq_100KB.csv"))[, -c(1, 2)]
RT = fread(pff("/data/001F_ENCODE_repliseq_100KBwindow_processed.csv"))[, -c(1, 2)]

Preds = as.data.table(cbind.data.frame(ATAC_seq[, .SD, .SDcols = -c(1,2)], 
									   Roadmap, 
									   Vierstra,
									   RT))

cohort_name = colnames(mutation_counts_dt)[cohort_index]

output = mutation_counts_dt[[cohort_index]]

#Remove hypermutated regions of the genome
ATAC_seq = ATAC_seq[-c(2858:2866, 20254, 24212:24219)]
Preds = Preds[-c(2858:2866, 20254, 24212:24219)]
output = output[-c(2858:2866, 20254, 24212:24219)]

dim(Preds)
# [1] 24440   773

rf = randomForest(x = Preds, 
				  y = output, 
				  keep.forest = T, 
				  ntree = 1000, 
				  do.trace = F, 
				  importance = F)

dt = as.data.table(cbind.data.frame(chr = ATAC_seq$chr,
									start = ATAC_seq$start, 
									observed = output, 
									predicted = rf$predicted))

fwrite(dt, paste0(pff("/data/005A_100KB_RF_errorDT_results/"), cohort_name, ".csv"))