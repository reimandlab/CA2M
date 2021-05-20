date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

#Load 1 MB mutation count dataset
mutation_counts_dt_MB = fread(pff("data/001D_PCAWG_mutation_counts_1MBwindow_processed.csv"))[ ,.SD, .SDcols = -c(1, 2)] 

mut_cor = cor(mutation_counts_dt_MB, method = "spearman")
pdf(pff("/data/001K_1MBmuttrack_cormap.pdf"), width = 10, height = 10)
heatmap.2(x = mut_cor, 
    Colv = T,
	Rowv = T,	  
    dendrogram = "both",
    col = "bluered",
    trace = "none",
    xlab = "Cancer Type",
	ylab = "Cancer Type",
	margins = c(13,13),
	cexRow = 1.3, cexCol = 1.3)
dev.off()

#Load 100KB mutation count dataset
mutation_counts_dt_100KB = fread(pff("data/001D_PCAWG_mutation_counts_100KBwindow_processed.csv"))[ ,.SD, .SDcols = -c(1, 2)] 

mut_cor = cor(mutation_counts_dt_100KB, method = "spearman")
pdf(pff("/data/001K_100KBmuttrack_cormap.pdf"), width = 10, height = 10)
heatmap.2(x = mut_cor, 
    Colv = T,
	Rowv = T,	  
    dendrogram = "both",
    col = "bluered",
    trace = "none",
    xlab = "Cancer Type",
	ylab = "Cancer Type",
	margins = c(13,13),
	cexRow = 1.3, cexCol = 1.3)
dev.off()