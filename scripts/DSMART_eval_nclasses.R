
####################################################################################################
#                                                                                                  #
#                        DSMART  - number of classes predicted map(s)                              #
#                                                                                                  #
#  In practice, the most probable soil is often a close race between several possible classes.     #
#  This can happen when a number of soil classes are very similar, or where                        #
#  DSMART is having trouble making effective predictions. A complementary measure to the confusion #
#  index, this function highlights the number of different soils predicted on a given pixel.       #
#  Ideally, this number is low. DSMART outputs also have a certain amount of 'noise' - soils       #
#  predicted only a few times on a given pixel during the `n` model runs, which can perhaps safely #
#  be disregarded.                                                                                 #
#                                                                                                  #
#  Usage                                                                                           #
#  eval_nclasses(class_predictions = NULL, cpus = 1, noise_cutoff = 0.1, nreals = NULL)            #
#                                                                                                  #
#  Arguments                                                                                       #
#  class_predictions  RasterStack; 'class_predictions_ranked.tif'; the tallied and sorted soil     #
#                     classes from DSMART_tally_reals().                                           #
#                                                                                                  #
#  cpus               Integer; number of compute nodes to use. Default is 1.                       #
#                                                                                                  #
#  noise_cutoff       Numeric; noise threshold e.g. 0.1 = soil classes predicted in less than      #
#                     10% of model runs are not counted.                                           #
#                                                                                                  #
#  nreals             Integer; The number of model runs.                                           #
#                                                                                                  #
####################################################################################################

eval_nclasses <- function(class_predictions = NULL, cpus = 1, noise_cutoff = 0.1, nreals = NULL) {
  
  if (is.null(nreals)) {
    stop('total number of DSMART model runs must be supplied.')
  }
  
  if (!dir.exists("dsmartOuts/evaluation")) {
      dir.create("dsmartOuts/evaluation", showWarnings = F)
  }
  strs    <- paste0(getwd(), "/dsmartOuts/evaluation/")
  
  # per-cell functions
  n_classes_predicted <- function(x) { 
    if (is.na(sum(x))) {
      NA
    } else {
      length(x[x > 0])
    }
  }
  
  # as above, with 'noise' removed - only classes predicted on >x% of model runs
  n_classes_predicted_x <- function(x) { 
    if (is.na(sum(x))) {
      NA
    } else {
      length(x[x > (ceiling(nreals * noise_cutoff))])
    }
  }
  
  beginCluster(cpus)
  assign('noise_cutoff', noise_cutoff, envir = .GlobalEnv)
  assign('nreals', nreals, envir = .GlobalEnv)
  
  npred <- clusterR(class_predictions,
                    fun = calc,
                    args = list(fun = n_classes_predicted),
                    filename = paste0(strs, 'n_classes_predicted.tif'),
                    NAflag = -9999,
                    datatype = 'INT2S',
                    overwrite = TRUE)
  
  assign('n_classes_predicted', npred,   envir = .GlobalEnv)
  
  npredx <- clusterR(class_predictions,
                     fun = calc, args = list(fun = n_classes_predicted_x),
                     export = c('nreals', 'noise_cutoff'),
                     filename = paste0(strs, 'n_classes_predicted_over_', (nreals * noise_cutoff), '.tif'),
                     datatype = 'INT2S',
                     NAflag = -9999,
                     overwrite = TRUE)
  
  assign('n_classes_predicted_x',  npredx, envir = .GlobalEnv)
  endCluster()
  
}
