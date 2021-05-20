date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

source(pff("/bin/999_process_data.R"))

input_data_dir = "/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/INPUT_DATA/"

#Import Roadmap supplementary data table
EID_table = fread(paste0(input_data_dir, 
		"/Roadmap_Consolidated_EpigenomeIDs_summary_Table.csv"))

#Import paths to Roadmap bigwigs
DNase_Seq_paths = list.files(paste0(input_data_dir, "Roadmap_DNAse/"), full.names = T)
DNase_Seq_names = unlist(lapply(list.files(paste0(input_data_dir, "Roadmap_DNAse/")),
	function(x) unlist(strsplit(x, split = "-"))[1]))
                              
#Import human genome build hg38                            
genome = BSgenome.Hsapiens.UCSC.hg38
                            
#Import liftover files
path = system.file(package = "liftOver", "extdata", "hg19ToHg38.over.chain")
ch = import.chain(paste0(input_data_dir, "/hg19ToHg38.over.chain"))

#Function to average individual tracks                           
merge_tracks = function(bigwig_path, window_size){

	gr.data = import(bigwig_path)

	#Liftover to hg38
	gr.data_hg38 = unlist(liftOver(gr.data, ch))
	
	#Use function from script to convert to windowed track
    windowed_track = create_windowed_track(gr.data_hg38, 
										   window_size)[[3]] 
	
    return(windowed_track)

}                            
                            
#Function to process and merge tracks into one data table                            
process_roadmap_data = function(paths, names, window_size){
    
    gr.windows = tileGenome(seqinfo(Hsapiens)[paste0("chr", 1:22)], 
							tilewidth = window_size, 
							cut.last.tile.in.chrom = TRUE)
    
    scores = as.data.table(do.call("cbind.data.frame", 
								   mclapply(paths, 
											function(x) merge_tracks(x, window_size), 
											mc.cores = 8)))
    colnames(scores) = names
    scores_dt = as.data.table(cbind.data.frame(chr = as.character(seqnames(gr.windows)), 
											   start = start(gr.windows), 
											   scores))
    
    scores_mappable = filter_windows_with_UMAP(scores_dt, 
											   window_size)

    colnames(scores_mappable)[-c(1,2)] = EID_table[[15]][match(colnames(scores_mappable)[-c(1,2)], 
															   EID_table[[2]])]
    
    return(scores_mappable)
}                            

fwrite(process_roadmap_data(DNase_Seq_paths, DNase_Seq_names, 1000000) , 
	   pff("/data/001A_Roadmap_DNaseSeq_1MBwindow_processed.csv"))
fwrite(process_roadmap_data(DNase_Seq_paths, DNase_Seq_names, 100000)  , 
	   pff("/data/001A_Roadmap_DNaseSeq_100KBwindow_processed.csv"))