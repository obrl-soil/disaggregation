####################################################################################################
#                                                                                                  #
#                                 DSMART  - confusion index map                                    #
#                                                                                                  #
#  The confusion index is a way to assess DSMART model performance. The CI is a measure of the     #
#  probability gap between the first and second most probable soils. Ideally, the gap is large,    #
#  implying that the most-probable soil is far more probable than any other.                       #
#                                                                                                  #
#  Usage                                                                                           #
#  eval_confind(class_probabilities = NULL, cpus = 1)                                              #
#                                                                                                  #
#  Arguments                                                                                       #
#  class_probabilities  RasterStack; 'class_probabilities_ranked.tif'; the tallied and sorted soil #
#                       class probabilites from DSMART_tally_reals().                              #
#                                                                                                  #
#  cpus                 Numeric; number of compute nodes to use. Default is 1.                     #
#                                                                                                  #
####################################################################################################

eval_confind <- function(class_probabilities = NULL, cpus = 1) {
  
  conf_index  <- function(x) {1 - (x[[1]] - x[[2]])}

  if (!dir.exists("dsmartOuts/evaluation")) {
      dir.create("dsmartOuts/evaluation", showWarnings = F)
  }
  strs    <- file.path(getwd(), 'dsmartOuts', 'evaluation')
  
  beginCluster(cpus)
  confus_ind  <- clusterR(na.omit(class_probabilities),
                          fun       = conf_index,
                          filename  = file.path(strs, 'confusion_index.tif'),
                          datatype  = 'FLT4S',
                          NAflag    = -9999,
                          overwrite = TRUE)
  assign('confus_ind', confus_ind,   envir = .GlobalEnv)
  endCluster()
  
}
