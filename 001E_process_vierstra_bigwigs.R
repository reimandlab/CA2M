date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

source(pff("/bin/999_process_data.R"))
	   
input_data_dir = "/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/INPUT_DATA/"

#Get data paths
track_names = list.files(paste0(input_data_dir, "vierstra_dnaseseq_bigwigs/"))
track_paths = list.files(paste0(input_data_dir, "vierstra_dnaseseq_bigwigs/"), 
						 full.names = T)

genome = BSgenome.Hsapiens.UCSC.hg38

get_vierstra_dataset = function(window_size){
	
	gr.windows = tileGenome(seqinfo(Hsapiens)[paste("chr", 1:22, sep = "")], 
							tilewidth = window_size, 
							cut.last.tile.in.chrom = TRUE)

	merge_tracks = function(track_index){

		obs_path = paste(track_paths[track_index], "/interval.all.obs.bw", sep = "")
		exp_path = paste(track_paths[track_index], "/interval.all.exp.bw", sep = "")

		obs.gr.data = import(obs_path)
		exp.gr.data = import(exp_path)

		#Get unified ranges
		overlaps = findOverlaps(obs.gr.data,exp.gr.data)
		unified_ranges = obs.gr.data[queryHits(overlaps)]
		FC = obs.gr.data[queryHits(overlaps)]$score / exp.gr.data[subjectHits(overlaps)]$score

		unified_ranges$score = FC

		windowed_track = create_windowed_track(unified_ranges,window_size)[[3]] #Use function from script to convert to windowed track

	return(windowed_track)

    }

    DNaseSEQ_scores = as.data.table(do.call("cbind", 
											mclapply(1:length(track_names), 
													 merge_tracks, 
													 mc.cores=16)))
    colnames(DNaseSEQ_scores) = track_names
    DNase_SEQ_dt = as.data.table(cbind.data.frame(chr = as.character(seqnames(gr.windows)), 
												  start = start(gr.windows), 
												  DNaseSEQ_scores))

    #Filter out windows with less than or equal to 80% mappability
    DNase_SEQ_dt_mappable = filter_windows_with_UMAP(DNase_SEQ_dt, window_size)
	
    return(DNase_SEQ_dt_mappable)}

fwrite(get_vierstra_dataset(1000000), pff("data/001E_vierstra_dnaseseq_1MB.csv"))
fwrite(get_vierstra_dataset(100000), pff("/data/001E_vierstra_dnaseseq_100KB.csv"))

