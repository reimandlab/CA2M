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

#Load 1 MB mutation count dataset
mutation_counts_dt = fread(pff("data/001D_PCAWG_mutation_counts_1MBwindow_processed.csv"))[,-c(1, 2)]

Preds_1 = as.data.table(cbind.data.frame(ATAC_seq, Roadmap_DNase, Vierstra))
Preds_2 = as.data.table(cbind.data.frame(ATAC_seq, Roadmap_DNase, Vierstra, RT))

adjust_R2 = function(R2, n, k){
	#Function to convert R2 to adjusted R2
	#R2: R-squared accuracy of model
	#n: Number of samples in training set
	#k: Number of predictors in model
	
	adjusted_R2 = 1 - (1 - R2)*(n - 1)/(n - k - 1)
	
	return(adjusted_R2)
}

get_adj_R2s = function(cohort_index, preds){
	
	output = mutation_counts_dt[[cohort_index]]
	
	rf = randomForest(preds, output, n_tree = 500, do.trace = F)
	
	R2 = (cor(rf$predicted, output))^2
	
	adj_R2 = adjust_R2(R2, nrow(preds), ncol(preds))
	
	return(adj_R2)
}

adj_R2_1 = mclapply(1:26, get_adj_R2s, Preds_1, mc.cores = 8)
adj_R2_2 = mclapply(1:26, get_adj_R2s, Preds_2, mc.cores = 8)

plot_dt = as.data.table(cbind.data.frame(Cancer_type = colnames(mutation_counts_dt), 
									   Adj_R2_diff = unlist(adj_R2_2)- unlist(adj_R2_1)))

plot_dt$Cancer_type = factor(plot_dt$Cancer_type, 
							 levels = plot_dt$Cancer_type[order(plot_dt$Adj_R2_diff, decreasing = T)])

pdf(pff("/data/001M_CAvsCART_barplots.pdf"))
ggplot(plot_dt, aes(x = Cancer_type, y = Adj_R2_diff))+
	geom_bar(stat = "identity", fill = "steelblue")+
	theme(axis.text.x = element_text(angle = 45, hjust = 1))+
	labs( x = "Cancer Type", y = "Diff. Adj. R2\n(CA+RT model - CA model)")
	
dev.off()