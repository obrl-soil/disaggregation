####################################################################################################
#                                                                                                  #
#                        DSMART  - individual soil class probability maps                          #
#                                                                                                  #
#  This function writes probability surfaces for each predicted soil class. It requires a 'probs'  #
#  output from DSMART_tally_reals() (use option `keep_tallies = TRUE`).                            #
#                                                                                                  #
#  Usage                                                                                           #
#  class_probs(probstack = NULL, lookup = NULL, cpus = 1)                                          #
#                                                                                                  #
#  Arguments                                                                                       #
#  probstack    RasterStack; 'class_props.tif'; output from DSMART_tally_reals().                  #                                                              #
#                                                                                                  #
#  lookup       Data frame; contains the full range of predicted ID (num) and CLASS (chr) values.  #
#               Can be retrieved from one of the realisation RATs, or (more reliably) imported     #
#               from the .txt file in /dsmartOuts/rasters.                                         #
#                                                                                                  #
#  cpus         Numeric; number of compute nodes to use. Default is 1.                             #
#                                                                                                  #
####################################################################################################

class_probs <- function(probstack = NULL, lookup = NULL, cpus = 1) {
  
  if (!dir.exists('dsmartOuts/summaries/class_maps/')) {
    dir.create('dsmartOuts/summaries/class_maps/', showWarnings = F)
  }
  class_dir <- paste0(getwd(), '/dsmartOuts/summaries/class_maps/')

  beginCluster(cpus)
  suppressWarnings(
    probs_list <- unstack(probs)
    )
    
  for (i in 1:length(probs_list)) {
    probs_list[[i]]@data@names <- as.character(lookup[i, 2])
  }
    
  prob_tifs <- lapply(probs_list, 
                      FUN = function(x) {
                        writeRaster(x, 
                                    filename = paste0(class_dir, names(x), "_probability.tif"),
                                    format = "GTiff",
                                    NAflag = -9999,
                                    datatype = 'FLT4S',
                                    overwrite = TRUE)}
                      ) 
  
  endCluster()
}
