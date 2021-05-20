date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

source(pff("/bin/999_process_data.R"))

input_data_dir = "/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/INPUT_DATA/"

#Load in hg38 PCAWG file
PCAWG_mutations_dt_hg38 = fread(pff("data/001C_PCAWG_mutations_hg38.csv"))
PCAWG_mutations_dt_hg38.gr = GRanges(PCAWG_mutations_dt_hg38$Chromosome, 
									 IRanges(PCAWG_mutations_dt_hg38$Start_position, 
																		 PCAWG_mutations_dt_hg38$End_position))


#Load in Cohort sizes for PCAWG
PCAWG_sample_table = PCAWG_mutations_dt_hg38[,.(Cohort_size = uniqueN(Donor_ID)), by = Project_Code]
cohorts_to_keep = sort(PCAWG_sample_table[Cohort_size>30]$Project_Code)


#Load in human genome hg38
genome = BSgenome.Hsapiens.UCSC.hg38

#Convert genome (chromosomes 1 through 22) into windows of specified size and cut last window in each chromosome to stay within chromosome
gr.windows = tileGenome(seqinfo(Hsapiens)[paste("chr", 1:22, sep = "")], 
						tilewidth = 1000000, 
						cut.last.tile.in.chrom=TRUE)

get_mean_VAF_tracks = function(cohort){
	print(cohort)
	mut = PCAWG_mutations_dt_hg38[Project_Code == cohort]
    mut.gr = GRanges(mut$Chromosome, 
						  IRanges(mut$Start_position, 
								 mut$End_position))
	
    #Get number of counts of mutations overlapping each genomic window
    overlaps = findOverlaps(gr.windows, mut.gr)
	means = unlist(lapply(1:length(seqnames(gr.windows)), 
						 function(x)  mean(mut$i_VAF[subjectHits(overlaps)[which(queryHits(overlaps) == x)]], na.rm = T)))
						  
						  
	return(means)}
						  
mutation_VAF_track_dt = as.data.table(do.call("cbind", lapply(cohorts_to_keep, get_mean_VAF_tracks)))
mutation_VAF_track_dt = as.data.table(cbind.data.frame(as.character(seqnames(gr.windows)), 
														  start(gr.windows),
														  mutation_VAF_track_dt))
colnames(mutation_VAF_track_dt) =  c("chr", "start", cohorts_to_keep)

#Remove windows with less than or equal to 80% mappability
mutation_VAF_track_dt_mappable = filter_windows_with_UMAP(mutation_VAF_track_dt, 1000000)

fwrite(mutation_VAF_track_dt_mappable, pff("data/001P_mut_VAF_track_mappable.csv"))