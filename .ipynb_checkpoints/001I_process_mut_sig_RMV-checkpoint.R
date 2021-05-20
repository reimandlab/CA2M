#Get cohort number as argument
args = commandArgs(trailingOnly = TRUE)
cohort_index = as.integer(args[1])

date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

source(pff("bin/999_process_data.R"))

#Load in mutations file
MAF = fread(pff("/data/001G_PCAWG_MAF_SNP_withSBSSIG.csv"))[,.SD,.SDcols=c(2,3,4,42,43,46)]

cancer_type_sample_sizes = MAF[, .(number_of_distinct_donors = length(unique(Donor_ID))), 
							   by = Project_Code]

project_code = cancer_type_sample_sizes[number_of_distinct_donors > 25]$Project_Code[cohort_index]

print(project_code)

#Create directory
dir.create(paste0(pff("/data/001I_Sig_RMVs/"), project_code))

MAF = MAF[Project_Code==project_code]

#Get number of mutations for each signature
sigs_table = table(MAF$SBS_Signature)
sigs_table_percent = sigs_table*100/nrow(MAF)

#Create function to run RF
get_RMV = function(signature){

    #Process MAF into RMV
    RMV = create_RMV(MAF[SBS_Signature==signature], 
					 1000000)
    RMV_mappable = filter_windows_with_UMAP(RMV, 1000000)

    fwrite(RMV_mappable,paste0(pff("/data/001I_Sig_RMVs/"), project_code, "/", signature, ".csv"))
}

#Run for all mutations with >10000 mutations
lapply(names(which(sigs_table > 10000 & sigs_table_percent>5)), get_RMV)

RMV = create_RMV(MAF, 1000000)

RMV_mappable = filter_windows_with_UMAP(RMV, 1000000)

fwrite(RMV_mappable, paste(pff("/data/001I_Sig_RMVs/"), project_code, "/ALL.csv", sep = ""))



