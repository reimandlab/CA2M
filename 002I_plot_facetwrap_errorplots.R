date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

source(pff("/bin/999_run_randomforest_experiment.R"))

#Load in predictor datasets
Roadmap = fread(pff("/data/001A_Roadmap_DNaseSeq_1MBwindow_processed.csv"))[,-c(1, 2)]
Vierstra = fread(pff("/data/001E_vierstra_dnaseseq_1MB.csv"))[, -c(1,2)]
RT = fread(pff("/data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))[, -c(1, 2)]
ATAC_seq = fread(pff("/data/001B_TCGA_ATACSeq_1MBwindow_processed.csv"))[, -c(1,2)]

#Get normal predictors + RT
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
Vierstra_normal = Vierstra[, .SD, .SDcols = -(vierstra_cancer_cell_lines)]

Normal_Preds = as.data.table(cbind.data.frame(Roadmap, Vierstra_normal, RT))

#Get tumour predictors + RT
Tumour_Preds = as.data.table(cbind.data.frame(ATAC_seq, 
									   Vierstra[, .SD, .SDcols = c("h.renal.cell.carcinoma-DS26693", 
																	 "h.renal.cell.carcinoma-DS37973")],
									   RT))

#Adjust R2 function
adjust_R2 = function(R2, n, k){
	adjusted_R2 = 1 - (1 - R2)*(n - 1)/(n - k - 1)
	return(adjusted_R2)
}

#Load in data paths
tumour_paths = list.files(pff("data/002G_tumourCAplusRT_RF_obsvsexpected"), full.name = T)
normal_paths = list.files(pff("data/002H_normalCAplusRT_RF_obsvsexpected"), full.name = T)

#Get cancer types
cohort_names = unlist(lapply(list.files(pff("data/002G_tumourCAplusRT_RF_obsvsexpected")),
                           function(x) unlist(strsplit(x, split = ".csv"))[1]))

                           
cohorts_to_examine = c("PANCAN", "Breast-AdenoCa", "Prost-AdenoCA", "Ovary-AdenoCA", "Skin-Melanoma")
                           
#Load in errors and put into DT
get_dt = function(cohort_name, paths, Preds, type){
        
    cohort_index = which(cohort_names == cohort_name)
    
    data = fread(paths[cohort_index])
    
    R2 = (cor(data$observed,data$predicted))**2
    
    adjR2=adjust_R2(R2, nrow(Preds), ncol(Preds))  
    
    dt = as.data.table(cbind.data.frame(observed = data$observed, 
										predicted = data$predicted, 
										cancer_type = cohort_name, 
										predictor_set = paste0(cohort_name, 
															  "\n", type, " CA + RT\n", 
															  "Adj.R2=", 
															  round(adjR2, 2))))
    
    return(dt)

}

tumour_dt = as.data.table(do.call("rbind.data.frame", 
								  lapply(cohorts_to_examine, get_dt, tumour_paths, Tumour_Preds, "Tumour")))
										
normal_dt = as.data.table(do.call("rbind.data.frame", 
								  lapply(cohorts_to_examine, get_dt, normal_paths, Normal_Preds, "Normal")))
                           
full_dt = as.data.table(rbind.data.frame(tumour_dt, normal_dt))                           
                           
                           
#Plot facet wrapped scatter plots
full_dt$predictor_set = factor(full_dt$predictor_set, 
							   levels = unique(full_dt$predictor_set))

###################################Functions to change scales on facets                           
scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}
                           
CustomFacetWrap <- ggproto(
  "CustomFacetWrap", FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)
                           
facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
    shrink = facet_super$shrink,
    params = facet_super$params
  )
}                           

################################################################################################################                           
                           
                           
pdf(pff("data/002I_facetwrap_errorplots_2.pdf"), width = 5.5, height = 3)
ggplot(full_dt, aes(x = observed, y = predicted))+
    rasterise(geom_point(alpha = 0.4, size = 0.3), dpi = 400)+
    geom_smooth(method = 'loess', span = 0.9, size = 0.5)+                      
    theme_bw()+
    labs(x = "Observed mutations per Mb", y = "Predicted mutations per Mb")+
    theme(axis.title = element_text(size = 7),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 5, colour = "black", angle = 30, hjust = 1),
         strip.text = element_text(size = 6))+
    facet_wrap(~predictor_set, nrow = 2, scales = "free")+                       
    facet_wrap_custom(~predictor_set, nrow = 2, scales = "free", scale_overrides = list(                      
        scale_override(1, scale_y_continuous(breaks = seq(4000, 20000, 4000), limits=c(4000, 20000))),
        scale_override(2, scale_y_continuous(breaks = seq(300, 900, 200), limits=c(300, 950))),
        scale_override(3, scale_y_continuous(breaks = seq(100, 600, 100), limits=c(100, 600))),
        scale_override(4, scale_y_continuous(breaks = seq(200, 800, 200), limits=c(200, 800))),                      
        scale_override(5, scale_y_continuous(breaks = seq(500, 2000, 500), limits=c(300, 2000))),                      
        scale_override(6, scale_y_continuous(breaks = seq(4000, 20000, 4000), limits=c(4000, 20000))),
        scale_override(7, scale_y_continuous(breaks = seq(300, 900, 200), limits=c(300, 950))),
        scale_override(8, scale_y_continuous(breaks = seq(100, 600, 100), limits=c(100, 600))),
        scale_override(9, scale_y_continuous(breaks = seq(200, 800, 200), limits=c(200, 800))),                      
        scale_override(10, scale_y_continuous(breaks = seq(500, 2000, 500), limits=c(300, 2000)))))                              
dev.off() 