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
Roadmap = fread(pff("/data/001A_Roadmap_DNaseSeq_1MBwindow_processed.csv"))[,.SD,.SDcols=-c(1,2)]
Vierstra = fread(pff("/data/001E_vierstra_dnaseseq_1MB.csv"))[,.SD,.SDcols=-c(1,2)]
RT = fread(pff("/data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))[,-c(1, 2)]


#Keep only healthy cell lines
Roadmap = Roadmap[, .SD, .SDcols = -c(40,42,43,48)]
vierstra_cancer_cell_lines = c("HeLa_S3-DS24790", "K562-DS15363", "K562-DS16924", "M059J-DS20493", "NB4-DS12543", "RPMI.7951-DS20909", "SK-N-DZ-DS21405",
                      "SK-N-DZ-DS21387", "SkNSH-DS8482", "h.786-O-DS27192", "h.786-O-DS37199", "h.A172-DS24744", "h.A549-DS14289",
                      "h.A673-DS27751", "h.A673-DS27667", "h.ACHN-DS24547", "h.ACHN-DS24471", "h.Caki-2-DS24211", "h.Daoy-DS24567",
                      "h.G-401-DS24147", "h.HT29-DS25113", "h.HT29-DS34646", "h.HT29-DS25119", "h.HT29-DS38587", "h.HepG2-DS7764",
                      "h.HepG2-DS24838", "h.K562-DS52908", "h.Karpas-422-DS27648", "h.Karpas-422-DS27641", "h.MG_63-DS24640",
                      "h.NAMALWA-DS22823", "h.NAMALWA-DS22816", "h.NCI-H226-DS21975", "h.NCI-H226-DS21935", "h.Oci-Ly-7-DS27687",
                      "h.Oci-Ly-7-DS29219", "h.PC9-DS27166", "h.RKO-DS40362", "h.RKO-DS44414", "h.RPMI.8226-DS22919", "h.RPMI.8226-DS22913",
                      "h.SJCRH30-DS23123", "h.SJSA-1-DS24235", "h.SK-MEL-5-DS23131", "h.SW-480-DS38592", "h.b.lymphoblast.mm.1s-DS27522",
                      "h.renal.cell.carcinoma-DS37973", "h.renal.cell.carcinoma-DS26693")
Vierstra = Vierstra[, .SD, .SDcols = -(vierstra_cancer_cell_lines)]

Preds = as.data.table(cbind.data.frame(Roadmap, Vierstra, RT))


#Remove hypermutated windows from lymphoma and leukemia
if(project_code %in% c("Lymph-CLL", "Lymph-BNHL")){
    Preds = Preds[-c(292, 2438)]
    Output = Output[-c(292, 2438)]
}

#Run random forest with monte carlo cross-validation
RF_result = run_RF_and_get_pred_importances(Preds, Output, 1000, n_tree = 1000, cores = 16, train_split = 0.8)

#Save results
fwrite(RF_result, paste0(pff("/data/002E_normalCAplusRT_RF_Results/"), project_code, ".csv"))

print("Finished")