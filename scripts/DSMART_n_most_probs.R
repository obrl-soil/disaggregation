####################################################################################################
#                                                                                                  #
#                           DSMART  - make n most-probable soil maps                               #
#                                                                                                  #
#  This function takes the final outputs of DSMART_tally_reals() and unstacks the top n layers,    # 
#  writing them to file.                                                                           #
#                                                                                                  #
#  Usage                                                                                           #
#  n_most_probs(class_predictions = NULL, class_probabilities = NULL, lookup = NULL, nmaps = 1,    # 
#  cpus = 1)                                                                                       #
#                                                                                                  #
#  Arguments                                                                                       #
#  class_predictions    RasterStack; 'class_predictions_ranked.tif'; the tallied and sorted soil   #
#                       classes from DSMART_tally_reals().                                         #
#                                                                                                  #
#  class_probabilities  RasterStack; 'class_probabilities_ranked.tif'; the tallied and sorted soil #
#                       class probabilites from DSMART_tally_reals().                              #
#                                                                                                  #
#  lookup               Data frame; contains the full range of predicted ID (num) and CLASS (chr)  #
#                       values. Can be retrieved from one of the realisation RATs, or (more        #
#                       reliably) imported from the .txt file in /dsmartOuts/rasters.              #
#                                                                                                  #
#  nmaps                Integer; the number of most-probable soil maps desired.                    #
#                                                                                                  #
#  cpus                 Numeric; number of compute nodes to use. Default is 1.                     #
#                                                                                                  #
####################################################################################################

n_most_probs <- function(class_predictions = NULL, class_probabilities = NULL, lookup = NULL, nmaps = 1, cpus = 1) {
  
  if (!dir.exists('dsmartOuts/summaries/most_probable_maps/')) { 
    dir.create('dsmartOuts/summaries/most_probable_maps/', showWarnings = F)
  }
  mp_dir <- file.path(getwd(), 'dsmartOuts', 'summaries', 'most_probable_maps') 

  beginCluster(cpus)

  if (!is.null(class_predictions)) {
     # only unstack the number of maps requested
      suppressWarnings(
        map_list <- unstack(class_predictions[[1:nmaps]])
      )
  } else {
      stop('class_predictions raster stack has not been specified')
  }
  
  for (i in 1:nmaps) {
    ratify(map_list[[i]])
    levels(map_list[[i]])[[1]] <- dplyr::left_join(levels(map_list[[i]])[[1]], lookup,
                                                   by = c('ID'))
  }
  
  outs_list <- mapply(FUN = function(x, i) {
    writeRaster(x, 
                filename  = file.path(mp_dir, paste0('mostlikely_', i ,'.tif')),
                format    = 'GTiff',
                NAflag    = -9999,
                datatype  = 'INT2S',
                overwrite = TRUE)
    },
    x = map_list,
    i = seq_along(1:nmaps)
    )
  
  assign('most_prob_maps', outs_list, envir = .GlobalEnv)
  message(paste0(nmaps, ' most-likely soils maps produced.'))
      
  ### optional probability surfaces
  if(!is.null(class_probabilities)) {
    suppressWarnings(
      probmap_list <- unstack(class_probabilities[[1:nmaps]])
    )
  }
  
  probsouts_list <- mapply(FUN = function(x, i) {
    writeRaster(x,
                filename = file.path(mp_dir, paste0('mostlikely_prob_', i, '.tif')),
                format = 'GTiff',
                datatype = 'FLT4S',
                NAflag = -9999,
                overwrite = TRUE)
    },
    x = probmap_list, 
    i = seq_along(1:nmaps)
    )
  
  assign('most_prob_prob_maps', probsouts_list, envir = .GlobalEnv)
  message(paste0(nmaps, ' probability maps produced.'))
      
  endCluster()
}
