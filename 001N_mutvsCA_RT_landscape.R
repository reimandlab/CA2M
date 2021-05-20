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

Preds_2 = as.data.table(cbind.data.frame(ATAC_seq, Roadmap_DNase, Vierstra, RT))


#Function to make plot

make_plot = function(CA, mut, CA_label, mut_label){
	plot_dt = cbind.data.frame(chr = mut$chr, start = mut$start, CA = -1*scale(CA), mut = scale(mut[[3]]))
	plot_dt.m = melt(plot_dt, id.vars = c("chr", "start"))
	
	#Calculate spearman correlation
	rho_cor = signif(cor(plot_dt[plot_dt$chr == "chr2","CA"], 
			  plot_dt[plot_dt$chr == "chr2","mut"], 
			  method = "spearman"), digits = 3)

	p = ggplot(plot_dt.m[plot_dt.m$chr == "chr2",], aes(x = start/(1e6), y = value, color = variable))+
	geom_smooth(span=0.05, se = FALSE)+
	theme_bw()+
	labs(x = "Chromosomal coordinate (mb)", 
		 y = "Scaled value", 
		 title = "Inverse chromatin accessibility vs. mutation density (chr2)",
		 color = "Data type")+
	scale_color_manual(labels = c(paste(CA_label,"chromatin accessibility\n(inverse scale)") , paste(mut_label,"mutation density")), 
					   values = c("Cyan", "Red"))+
	annotate("text", label = paste0("r=", rho_cor), color = "red", x = 10, y = 1.9)
	
	return(p)
}

####Melanoma 7a vs melanocyte

#Load in mutations file
mel_sig = fread("/.mounts/labs/reimandlab/private/users/oocsenas/ATACSEQ_MUT_THESIS/Random_Forest/data/240520/Sig_RMVs/Skin-Melanoma/SBS7a.csv")

#Compare with top predictor
top_pred=names(sort(sapply(colnames(Preds), function(x) cor(Preds[[x]], mel_sig[[3]], method="spearman"))))[1]
melanocyte_CA = Preds[[top_pred]]


pdf("/.mounts/labs/reimandlab/private/users/oocsenas/ChrAcc2Mut/data/001H_melanocyteCA_vs_melanoma7Amut.pdf",width=10,height=4)
make_plot(melanocyte_CA, mel_sig, "Melanocyte", "Melanoma SBS7a")
dev.off()

####Breast AdenoCA vs BRCA
#Load in mutations file
breast_mut = fread("/.mounts/labs/reimandlab/private/users/oocsenas/ATACSEQ_MUT_THESIS/Random_Forest/data/240520/Sig_RMVs/Breast-AdenoCa/ALL.csv")

#Compare with top predictor
top_pred=names(sort(sapply(colnames(Preds), function(x) cor(Preds[[x]], breast_mut[[3]], method="spearman"))))[1]
BRCA_CA = Preds[[top_pred]]


pdf("/.mounts/labs/reimandlab/private/users/oocsenas/ChrAcc2Mut/data/001H_BRCACA_vs_Breastmut.pdf",width=10,height=4)
make_plot(BRCA_CA, breast_mut, "BRCA", "Breast-AdenoCa")
dev.off()


####GBM vs LGG
#Load in mutations file
GBM_mut = fread("/.mounts/labs/reimandlab/private/users/oocsenas/ATACSEQ_MUT_THESIS/Random_Forest/data/240520/Sig_RMVs/CNS-GBM/ALL.csv")

#Compare with top predictor
top_pred=names(sort(sapply(colnames(Preds), function(x) cor(Preds[[x]], GBM_mut[[3]], method="spearman"))))[1]
LGG_CA = Preds[[top_pred]]


pdf("/.mounts/labs/reimandlab/private/users/oocsenas/ChrAcc2Mut/data/001H_LGGCA_vs_GBMmut.pdf",width=10,height=4)
make_plot(LGG_CA, GBM_mut, "LGG", "GBM")
dev.off()
						   
####### STAD vs STAD
STAD_mut = fread("/.mounts/labs/reimandlab/private/users/oocsenas/ATACSEQ_MUT_THESIS/Random_Forest/data/240520/Sig_RMVs/Stomach-AdenoCA/ALL.csv")

#Compare with top predictor
top_pred=names(sort(sapply(colnames(Preds), function(x) cor(Preds[[x]], STAD_mut[[3]], method="spearman"))))[1]
top_CA = Preds[[top_pred]]


pdf("/.mounts/labs/reimandlab/private/users/oocsenas/ChrAcc2Mut/data/001H_STADCA_vs_STADmut.pdf",width=10,height=4)
make_plot(top_CA, STAD_mut, "STAD", "Stomach-AdenoCA")
dev.off()
						   
####### Lymphoma
BNHL_mut = fread("/.mounts/labs/reimandlab/private/users/oocsenas/ATACSEQ_MUT_THESIS/Random_Forest/data/240520/Sig_RMVs/Lymph-BNHL/ALL.csv")

#Compare with top predictor
top_pred=names(sort(sapply(colnames(Preds), function(x) cor(Preds[[x]], BNHL_mut[[3]], method="spearman"))))[1]
top_CA = Preds[[top_pred]]


pdf("/.mounts/labs/reimandlab/private/users/oocsenas/ChrAcc2Mut/data/001H_BcellCA_vs_BNHLmut.pdf",width=10,height=4)
make_plot(top_CA, STAD_mut, "Primary B cell", "Lymph-BNHL")
dev.off()
						   
						   
##############Make separate plots
####### STAD vs STAD
STAD_mut = fread("/.mounts/labs/reimandlab/private/users/oocsenas/ATACSEQ_MUT_THESIS/Random_Forest/data/240520/Sig_RMVs/Stomach-AdenoCA/ALL.csv")

#Compare with top predictor
top_pred=names(sort(sapply(colnames(Preds), function(x) cor(Preds[[x]], STAD_mut[[3]], method="spearman"))))[1]
top_CA = Preds[[top_pred]]
						   
plot_dt = cbind.data.frame(chr = STAD_mut$chr, start = STAD_mut$start, CA = (top_CA), mut = (STAD_mut[[3]]))
plot_dt.m = melt(plot_dt, id.vars = c("chr", "start"))

#Make CA line plot
pdf("ChrAcc2MutData/mutvsCA_landscape_plots/STAD_CA_lineplot.pdf",width=10,height=4)						   
ggplot(plot_dt[plot_dt$chr == "chr2",], aes(x = start/(1e6), y = CA))+
	geom_smooth(span=0.05, se = FALSE)+
	theme_bw()+
	labs(x = "Chromosomal coordinate (mb)", 
		 y = "Chromatin accessibility", 
		 title = "STAD chromatin accessibility (chr2)")
dev.off()
						   
#Make CA line plot
pdf("ChrAcc2MutData/mutvsCA_landscape_plots/STAD_mut_barplot.pdf",width=10,height=4)						   
ggplot(plot_dt[plot_dt$chr == "chr2",], aes(x = start/(1e6), y = mut))+
	geom_bar(stat = "identity", color = "black", fill = "red")+
	theme_bw()+
	labs(x = "Chromosomal coordinate (mb)", 
		 y = "Total mutations", 
		 title = "Stomach AdenoCa mutation density (chr2)")
dev.off()
						   
