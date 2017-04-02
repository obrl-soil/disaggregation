####################################################################################################
#                                                                                                  #
#                                 DSMART  - tally model realisations                               #
#                                                                                                  #
#  This is the first step in producing DSMART output maps - all model realisations are stacked and #
#  the number of times a given soil is predicted on a given pixel tallied. The tally is then       #
#  sorted descending, producing a stack of most-to-least probable soil classes. The actual         #
#  probabilities are calculated as well.                                                           #
#                                                                                                  #
#  Usage                                                                                           #
#  tally_reals(realstack = NULL, lookup = NULL, cpus = 1, keep_tallies = TRUE)                     #
#                                                                                                  #
#  Arguments                                                                                       #
#  realstack    RasterStack; A stack of model realisations from DSMART_AP().                       #
#                                                                                                  #
#  lookup       Data frame; contains the full range of predicted ID (num) and CLASS (chr) values.  #
#               Can be retrieved from one of the realisation RATs, or (more reliably) imported     #
#               from the .txt file in /dsmartOuts/rasters.                                         #
#                                                                                                  #
#  pid_field    String; the name of the numeric ID field in indata.                                #
#                                                                                                  #
#  cpus         Numeric; number of compute nodes to use. Default is 1.                             #
#                                                                                                  #
#  Keep_tallies Logical; default TRUE. Writes the unordered tally stacks to file - for FALSE they  #
#               will only be retained in R's temp location until the session is ended.             #
#                                                                                                  #
####################################################################################################

tally_reals <- function(realstack = NULL, lookup = NULL, cpus = 1, keep_tallies = TRUE) {
  
  if (!dir.exists("dsmartOuts/summaries")) {
      dir.create("dsmartOuts/summaries", showWarnings = F)
  }
  strs    <- paste0(getwd(), "/dsmartOuts/summaries/")
  classes <- nrow(lookup)
  realn   <- nlayers(realstack)
  
  ### cell by cell calc functions ###
  
  # class_count produces an integer vector counting the number of times a given soil was 
  # predicted on this pixel
  # e.g. 0 5 8 0 == soil 1 was not predicted, soil 2 was predicted 5 times, soil 3 was predicted 
  # 8 times, soil 4 was not predicted.
  class_count <- function(x) { 
    if (is.na(sum(x))) {
      rep(NA, classes)
    } else {
      tabulate(x, nbins = classes)
    }
  }
  
  # class_prop takes the vector produced by class_count and normalises it against the 
  # total number of model runs
  # e.g. for 20 runs, 0 5 8 0 becomes 0 0.25 0.4 0
  class_prop  <- function(x) {round((x / realn), 3)}
  
  # stack_cell_sort sorts the elements of class_order or probs_order from largest to smallest.
  # e.g. 0 0.25 0.4 0 becomes 0.4 0.25 0 0
  stack_cell_sort <- function(x) {
    if (is.na(sum(x))) {
      rep(NA, classes)
    } else { 
      sort(x, decreasing = TRUE, na.last = TRUE)
    }
  }
  
  # stack_cell_order returns a vector of the position of the elements of probs_order
  # from largest to smallest.
  # e.g. 0 0.25 0.4 0 becomes [1] 3 2 1 4
  # these values correspond to the lookup ID column, so they can be linked to soil class name
  stack_cell_order <- function(x) {
    if (is.na(sum(x))) {
      rep(NA, classes)
    } else { 
      order(x, decreasing = TRUE, na.last = TRUE)
    }
  }
  
  # note that its easiest to think about the results of stack_cell_order as the names
  # of stack_cell_sort  e.g. [1] 3 2 1 4 or [1]  3   2   1 4
  #                              8 5 0 0        0.4 0.25 0 0
  
  ###
  beginCluster(cpus)
  assign('classes', classes, envir = .GlobalEnv)
  assign('realn',   realn,   envir = .GlobalEnv)
  
  pb <- txtProgressBar(min = 0, max = 100, initial = 0, style = 3)
  
  # 1. get a stack of predicted class counts (NA cells accounted for)
  counts <- clusterR(x         = realstack,
                     fun       = calc, 
                     args      = list(fun = class_count), 
                     export    = 'classes', 
                     datatype  = 'INT2S')
  
  if (keep_tallies == TRUE) { 
    writeRaster(counts,
                filename  = paste0(strs, 'class_counts.tif'),
                datatype  = 'INT2S',
                format    = 'GTiff',
                NAflag    = -9999,
                overwrite = TRUE)
    assign('counts', counts, envir = .GlobalEnv) 
  }

  setTxtProgressBar(pb, 40)

  # 2. express that stack as a probability
  probs <- clusterR(x         = counts,
                    fun       = calc,
                    args      = list(fun = class_prop),
                    export    = 'realn',
                    datatype  = 'FLT4S')
  
  if (keep_tallies == TRUE) { 
    writeRaster(probs,
                filename  = paste0(strs, 'class_props.tif'),
                datatype  = 'FLT4S',
                format    = 'GTiff',
                NAflag    = -9999,
                overwrite = TRUE)
    assign('probs', probs, envir = .GlobalEnv)
  }

  setTxtProgressBar(pb, 60)
 
  # 3. sort probs by most to least probable
  sorted  <- clusterR(x         = probs,
                      fun       = calc, 
                      args      = list(fun = stack_cell_sort),
                      filename  = paste0(strs, 'class_probabilities_ranked.tif'),
                      datatype  = 'FLT4S',
                      format    = 'GTiff',
                      NAflag    = -9999,
                      overwrite = TRUE)
  assign('sorted', sorted, envir = .GlobalEnv)
  
  setTxtProgressBar(pb, 80)
  
  # 4. order probs by most to least probable
  ordered <- clusterR(x        = na.omit(probs),
                      fun      = calc, 
                      args     = list(fun = stack_cell_order),
                      datatype = 'INT2S')

 # give the raster names for these values
  for (i in 1:nlayers(ordered)) { 
    levels(ordered[[i]]) <- lookup
  }
  writeRaster(ordered,
              filename  = paste0(strs, 'class_predictions_ranked.tif'),
              dataype   = 'INT2S',
              format    = 'GTiff',
              NAflag    = -9999,
              overwrite = TRUE)
  assign('ordered', ordered, envir = .GlobalEnv)
  
  if (keep_tallies == FALSE) { 
    rm(counts)
    rm(probs) 
    }

  setTxtProgressBar(pb, 100)
  close(pb)
  endCluster()

}
