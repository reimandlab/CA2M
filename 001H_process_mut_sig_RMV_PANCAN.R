date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

source(pff("bin/999_process_data.R"))

#Load in mutations file
MAF = fread(pff("/data/001G_PCAWG_MAF_SNP_withSBSSIG.csv"))[,.SD,.SDcols=c(2,3,4,46)]

#Create directory
dir.create(pff("/data/Sig_RMVs/PANCAN"))

#Get number of mutations for each signature
sigs_table = table(MAF$SBS_Signature)
sigs_table_percent = sigs_table*100/nrow(MAF)

#Create function to run RF
get_RMV = function(signature){

    #Process MAF into RMV
    RMV = create_RMV(MAF[SBS_Signature == signature], 1000000)
    RMV_mappable = filter_windows_with_UMAP(RMV, 1000000)

    fwrite(RMV_mappable, paste(pff("/data/Sig_RMVs/PANCAN/"), signature, ".csv", sep = ""))
}

#Run for all mutations with minimum 10,000 mutations and >5% of total mutations in cohort
lapply(names(which(sigs_table > 10000 & sigs_table_percent>5)), get_RMV)

RMV = create_RMV(MAF,1000000)

RMV_mappable = filter_windows_with_UMAP(RMV, 1000000)

fwrite(RMV_mappable,pff("/data/Sig_RMVs/PANCAN/ALL.csv"))



