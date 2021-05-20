date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

#Load in 1 MB RT and CA
RT = fread(pff("/data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))[,-c(1,2)]
CA = fread(pff("/data/001E_vierstra_dnaseseq_1MB.csv"))[,-c(1,2)]                                                         
                                          										
#Plot correlation
cor_mat = cor(CA, RT)
										
#Create cofactor bars
Cell_stage = rep("", 96)
Cell_stage[grep("G1", colnames(cor_mat))] = "G1"                                          
Cell_stage[grep("G2", colnames(cor_mat))] = "G2"                                          
Cell_stage[grep("S1", colnames(cor_mat))] = "S1"                                          
Cell_stage[grep("S2", colnames(cor_mat))] = "S2"                                          
Cell_stage[grep("S3", colnames(cor_mat))] = "S3"                                          
Cell_stage[grep("S4", colnames(cor_mat))] = "S4"                                          

Cell_type = rep("", 96)
Cell_type[1:6] = "Stem_cells"                                                                                  
Cell_type[7:18] = "Skin"
Cell_type[19:48] = "Blood"
Cell_type[49:54] = "Cervix"
Cell_type[55:60] = "Liver"
Cell_type[61:66] = "Endothelial"
Cell_type[67:72] = "Lung"
Cell_type[73:78] = "Blood"
Cell_type[79:84] = "Breast"
Cell_type[85:90] = "Skin"
Cell_type[91:96] = "Brain"
                                          
rowside_mat = cbind.data.frame(Cell_stage, Cell_type)                                        
rownames(rowside_mat) = colnames(cor_mat)                                    
                                          
my_palette = colorRampPalette(c("blue", "white", "red"))(n = 299)

# Specify colors
ann_colors = list(
    Cell_stage = c(G1 = "Red", G2 = "Orange", S1 = "Yellow", S2 = "Green", S3 = "Blue", S4 = "Purple"),
    Cell_type =c(Stem_cells = "#A6CEE3", Skin = "#1F78B4", Blood = "#B2DF8A", Cervix = "#33A02C", 
				 Liver = "#FB9A99", Endothelial = "#E31A1C", Lung = "#FDBF6F", Breast = "#FF7F00", Brain = "#CAB2D6")
)
                                                     								
										
pdf(pff("/data/001L_RT_CA_cor_heatmap2.pdf"), width = 20, height = 50)
pheatmap(cor_mat, color = my_palette, breaks = seq(-1, 1, length.out = 300), 
		 annotation_col = rowside_mat, annotation_colors = ann_colors)
dev.off()  										


#Get auto correlation heatmap of RT
cor_mat = cor(RT, RT)

pdf(pff("/data/001L_RT_RT_cor_heatmap2.pdf"), width = 35, height = 30)
pheatmap(cor_mat, color = my_palette, breaks = seq(-1, 1, length.out = 300), 
		 annotation_col = rowside_mat, annotation_row = rowside_mat, annotation_colors = ann_colors)
dev.off()  