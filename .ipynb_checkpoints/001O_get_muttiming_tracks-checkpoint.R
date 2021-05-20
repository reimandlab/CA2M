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

#Load in PCAWG timing data
timing = fread(paste0(input_data_dir, "PCAWG_timing/2018-07-19-allSegmentsTimeRaw.icgc.controlled.txt"))

#Process PCAWG timing data
timing$chromosome = paste0("chr", timing$seqnames)

#Change samples to proper IDs
timing$Donor_ID = PCAWG_mutations_dt_hg38$Donor_ID[match(timing$sample, 
														 PCAWG_mutations_dt_hg38$Tumor_Sample_Barcode)]


#Add timing to each mutation in MAF by samples
get_timing = function(sample_ID){
	
	
	mut_sample = PCAWG_mutations_dt_hg38[which(PCAWG_mutations_dt_hg38$Tumor_Sample_Barcode == sample_ID)]
	mut_sample.gr = PCAWG_mutations_dt_hg38.gr[which(PCAWG_mutations_dt_hg38$Tumor_Sample_Barcode == sample_ID)]

	timing_sample = timing[sample == sample_ID][!(is.na(time))]
	timing_sample.gr = GRanges(timing_sample$chromosome, 
							   IRanges(timing_sample$start, timing_sample$end))
	
	Times_muts_sample = timing_sample$time[queryHits(findOverlaps(timing_sample.gr, 
												  mut_sample.gr))]
	Mut_data_with_overlap = mut_sample[subjectHits(findOverlaps(timing_sample.gr, 
												  mut_sample.gr))]
	
	cols_to_keep = c(2, 3, 4, 13, 42, 43)
	dt = as.data.table(cbind.data.frame(Mut_data_with_overlap[, .SD, .SDcols = cols_to_keep], 
										time = Times_muts_sample))
	
	return(dt)
		
}

mut_timing_dt = as.data.table(do.call("rbind.data.frame", 
									  mclapply(unique(timing$sample), get_timing, mc.cores = 8)))

fwrite(mut_timing_dt, pff("001O_muttiming_MAF.csv"))

#Load in human genome hg38
genome = BSgenome.Hsapiens.UCSC.hg38

#Convert genome (chromosomes 1 through 22) into windows of specified size and cut last window in each chromosome to stay within chromosome
gr.windows = tileGenome(seqinfo(Hsapiens)[paste("chr", 1:22, sep = "")], 
						tilewidth = 1000000, 
						cut.last.tile.in.chrom=TRUE)

get_mut_timing = function(cohort){
	
	mut = mut_timing_dt[Project_Code == cohort]
    mut.gr = GRanges(mut$Chromosome, 
						  IRanges(mut$Start_position, 
								 mut$End_position))
	
    #Get number of counts of mutations overlapping each genomic window
    overlaps = findOverlaps(gr.windows, mut.gr)
	means = unlist(mclapply(1:length(seqnames(gr.windows)), 
						 function(x)  mean(mut$time[subjectHits(overlaps)[which(queryHits(overlaps) == x)]]), 
							mc.cores = 8))
						  
						  
	return(means)}
						  
mutation_timing_track_dt = as.data.table(do.call("cbind", lapply(cohorts_to_keep, get_mut_timing)))
mutation_timing_track_dt = as.data.table(cbind.data.frame(as.character(seqnames(gr.windows)), 
														  start(gr.windows),
														  mutation_timing_track_dt))
colnames(mutation_timing_track_dt) =  c("chr", "start", cohorts_to_keep)

#Remove windows with less than or equal to 80% mappability
mutation_timing_track_dt_mappable = filter_windows_with_UMAP(mutation_timing_track_dt, 1000000)

fwrite(mutation_timing_track_dt_mappable, pff("data/001O_mut_timing_track_mappable.csv"))