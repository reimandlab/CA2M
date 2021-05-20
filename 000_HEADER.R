date_tag = "210317"
home_dir = paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# home_dir = paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag)
setwd(home_dir)

packages = c("data.table", "ActivePathways", "parallel", "gridExtra", "grid", 
			 "ggplot2", "mgsub", "gplots", "RColorBrewer", "GenomicRanges", 
			 "rtracklayer", "BSgenome.Hsapiens.UCSC.hg38", "plyr", "randomForest", 
			 "caret", "rCGH", "ggrastr", "cowplot", "pheatmap", "ggsci", "ggpubr",
			 "ggrepel", "gtools")

lapply(packages, require, character.only = TRUE)

pff = function(x = "") {
	paste0("./", paste0(x, collapse = ""), sep = "")
}

opff = function(x) {
	system(paste("open", pff(x)))
}




file_open_call = function(fname) {
	fname1 = tail(strsplit(fname, '/')[[1]], 1)
	cat(paste0("\n\nscp -C $hn:", getwd(), "/", fname, " ./ && open ", fname1, "\n\n"))
}

file_open_call2 = function(fname) {
	system(paste("cp ", fname, "~/tmp/"))
	fname1 = tail(strsplit(fname, '/')[[1]], 1)
	cat(paste0("\n\n scp -C $cw:~/tmp/", fname1, " ./ && open ", fname1, "\n\n"))
}