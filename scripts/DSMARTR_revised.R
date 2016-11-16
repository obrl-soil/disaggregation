# this processes all the model realisations into final products. 

DSMARTR_2 <- function(realstack = NULL, n_mpmaps = 3, class_maps = FALSE, lookup = NULL, cpus = 1) {
  
  pb <- txtProgressBar(min = 0, max = 100, style = 3)
  classes <- nrow(lookup)
  realn   <- nlayers(realstack)
  
  # functions
  class_count <- function(x) { 
    if (is.na(sum(x))) {
      rep(NA, classes)
    } else {
      tabulate(x, nbins = classes)
    }
  }
  
  class_prop  <- function(x) { x / realn }
  
  class_order  <- function(x) {
    if (is.na(sum(x))) {
      rep(NA, classes)
    } else { 
      order(x, decreasing = TRUE, na.last = TRUE) 
    }
  }
  
  probs_order  <- function(x) {
    if (is.na(sum(x))) {
      rep(NA, classes)
    } else { 
      sort(x, decreasing = TRUE, na.last = TRUE) 
    }
  }
  
  conf_index  <- function(x) { 1 - (x[[1]] - x[[2]]) }
  
  #debilitatingperfectionism
  nth_fun <- function(x) {
    x <- abs(x)
    if (as.integer(substring(x, seq(nchar(x)), seq(nchar(x))))[nchar(x)] == 1) {
      return(paste0(x, 'st'))
    } else if (as.integer(substring(x, seq(nchar(x)), seq(nchar(x))))[nchar(x)] == 2) {
      return(paste0(x, 'nd'))
    } else if (as.integer(substring(x, seq(nchar(x)), seq(nchar(x))))[nchar(x)] == 3) {
      return(paste0(x, 'rd'))
    } else {
      return(paste0(x, 'th'))
    }
  }
  
  dir.create('dsmartOuts/summaries/rds_files/', showWarnings = F)
  dir.create('dsmartOuts/summaries/most_probable_maps/', showWarnings = F)
  rds_dir   <- paste0(getwd(), '/dsmartOuts/summaries/rds_files/')
  mp_dir    <- paste0(getwd(), '/dsmartOuts/summaries/most_probable_maps/') 
  
  #### process reals ####
  beginCluster(cpus)
  assign('classes', classes, envir = .GlobalEnv)
  assign('realn',   realn,   envir = .GlobalEnv)
  
  # 1. get a stack of predicted class counts (NA cells accounted for)
  # object 'counts' will require (2 * ncells * nreals) bytes of storage space
  counts <- clusterR(reals, calc, args = list(fun = class_count), export = 'classes', datatype = 'INT2S')
  assign('counts', counts, envir = .GlobalEnv)
  setTxtProgressBar(pb, 25)
  
  # 2. express that stack as a probability
  # object 'counts' will require (4 * ncells * nreals) bytes of storage space
  probs <- clusterR(na.omit(counts), calc, args = list(fun = class_prop), export = 'realn', datatype = 'FLT4S')
  assign('probs', probs, envir = .GlobalEnv)
  setTxtProgressBar(pb, 50)
  
  # 3. sort counts and probs by most to least probable, by cell, get confusion index for 1 and 2 mp
  counts_ord  <- clusterR(na.omit(counts), calc, args = list(fun = class_order), datatype = 'INT2S')
  probs_ord   <- clusterR(na.omit(probs),  calc, args = list(fun = probs_order), datatype = 'FLT4S')
  confus_ind  <- clusterR(na.omit(probs_ord), fun = conf_index, datatype = 'FLT4S')
  assign('counts_ord', counts_ord, envir = .GlobalEnv)
  assign('probs_ord',  probs_ord,  envir = .GlobalEnv)
  assign('confus_ind', confus_ind, envir = .GlobalEnv)
  setTxtProgressBar(pb, 75)
  
  # 4. Generate n_most probable soil maps and matching probability maps
  if (!is.null(n_mpmaps)) {
    if (n_mpmaps == 0) { 
      return(print("Warning: No most-probable maps will be generated"))
      } else {
        # retrieve n most probable
        suppressWarnings(nprobs_list <- unstack(counts_ord[[1:n_mpmaps]])) 
        
        # lapply fails here for no apparent reason so we're stuck with the for loop when setting
        # factor levels
          for (i in 1:length(nprobs_list)) {
            levels(nprobs_list[[i]]) <- lookup
          }
        
        nprob_tifs <- mapply(FUN = function(x, i) {
                              writeRaster(x, 
                                          filename  = paste0(mp_dir, nth_fun(i), "_mostprob.tif"),
                                          format    = "GTiff",
                                          NAflag    = -9999,
                                          datatype  = 'INT2S',
                                          overwrite = TRUE)}, 
                              x = nprobs_list,
                              i = seq_along(1:length(nprobs_list)))
        assign('nprob_maps', nprob_tifs, envir = .GlobalEnv)
        rm(nprobs_list)
                       
        # n most probable probabilities
        suppressWarnings(npps_list <- unstack(probs_ord[[1:n_mpmaps]]))
        npp_tifs <- mapply(FUN = function(x, i) {
          writeRaster(x, 
                      filename  = paste0(mp_dir, nth_fun(i), "_probable_probs.tif"),
                      format    = "GTiff",
                      NAflag    = -9999,
                      datatype  = 'FLT4S',
                      overwrite = TRUE)}, 
          x = npps_list,
          i = seq_along(1:length(npps_list)))
        assign('nprob_probs_maps', npp_tifs, envir = .GlobalEnv)
        rm(npps_list)
        
        # confusion index
        writeRaster(confus_ind, 
                    filename  = paste0(mp_dir, "confusion_index.tif"),
                    format    = "GTiff",
                    NAflag    = -9999,
                    datatype  = 'FLT4S',
                    overwrite = TRUE)
        }
    }
  
  # 5. If class_maps = T, make a probability map for each class
  if (class_maps == TRUE) {
    
    dir.create('dsmartOuts/summaries/class_maps/', showWarnings = F)
    class_dir <- paste0(getwd(), '/dsmartOuts/summaries/class_maps/')
    
    suppressWarnings(probs_list <- unstack(probs))
    
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
                                      overwrite = TRUE)}) 
    
    assign('prob_tifs', prob_tifs, envir = .GlobalEnv)
    rm(probs_list)
  }
  
  # Save clusterR outputs as rds objects where possible, GTiff where RAM limits prevent this
  # FLT4S .gri file sizes are about the same as object.size(readAll(object))
  if ((file.size(gsub('grd', 'gri', counts@file@name)) * 2) > (memory.size() * 1048576) {
    message("class counts object too large to fit in memory. Writing GTiff instead of RDS")
    dir.create('dsmartOuts/summaries/bigouts/', showWarnings = F)
    tif_dir   <- paste0(getwd(), '/dsmartOuts/summaries/bigouts/')
    writeRaster(counts, 
                filename = paste0(tif_dir, "class_counts.tif"),
                format = 'GTiff',
                NAflag = -9999,
                datatype = 'INT2S',
                overwrite = TRUE)
  } else {
  saveRDS(readAll(counts), paste0(rds_dir, 'class_counts.rds'))
  rm(counts)
  }
  
  if (file.size(gsub('grd', 'gri', probs@file@name)) > (memory.size() * 1048576)) { 
    message("probabilities object too large to fit in memory. Writing GTiff instead of RDS")
    dir.create('dsmartOuts/summaries/bigouts/', showWarnings = F)
    tif_dir   <- paste0(getwd(), '/dsmartOuts/summaries/bigouts/')
    writeRaster(probs, 
                filename = paste0(tif_dir, "class_probs.tif"),
                format = 'GTiff',
                NAflag = -9999,
                datatype = 'FLT4S',
                overwrite = TRUE)
  } else {
  saveRDS(readAll(probs), paste0(rds_dir, 'class_probs.rds'))
  rm(probs)
  }
  
  if ((file.size(gsub('grd', 'gri', counts_ord@file@name)) * 2) > (memory.size() * 1048576) {
    message("ordered classes object too large to fit in memory. Writing GTiff instead of RDS")
    dir.create('dsmartOuts/summaries/bigouts/', showWarnings = F)
    tif_dir   <- paste0(getwd(), '/dsmartOuts/summaries/bigouts/')
    writeRaster(counts_ordered, 
                filename = paste0(tif_dir, "class_counts_ordered.tif"),
                format = 'GTiff',
                NAflag = -9999,
                datatype = 'INT2S',
                overwrite = TRUE)
    
  } else {
  saveRDS(readAll(counts_ord), paste0(rds_dir, 'class_counts_ordered.rds'))
  rm(counts_ord)
  }
  
  if (file.size(gsub('grd', 'gri', probs_ord@file@name)) > (memory.size() * 1048576)) {
    message("ordered probabilities object too large to fit in memory. Writing GTiff instead of RDS")
    dir.create('dsmartOuts/summaries/bigouts/', showWarnings = F)
    tif_dir   <- paste0(getwd(), '/dsmartOuts/summaries/bigouts/')
    writeRaster(probs_ordered, 
                filename = paste0(tif_dir, "class_probs_ordered.tif"),
                format = 'GTiff',
                NAflag = -9999,
                datatype = 'FLT4S',
                overwrite = TRUE)
  } else {
  saveRDS(readAll(probs_ord),  paste0(rds_dir, 'class_probs_ordered.rds'))
  rm(probs_ord)
  }
  
  if (file.size(gsub('grd', 'gri', confus_ind@file@name)) > (memory.size() * 1048576)) {
    message("confusion index too large to fit in memory. Writing GTiff instead of RDS")
    dir.create('dsmartOuts/summaries/bigouts/', showWarnings = F)
    tif_dir   <- paste0(getwd(), '/dsmartOuts/summaries/bigouts/')
    writeRaster(confus_ind, 
                filename = paste0(tif_dir, "confusion_index_1_2.tif"),
                format = 'GTiff',
                NAflag = -9999,
                datatype = 'FLT4S',
                overwrite = TRUE)
  } else {
  saveRDS(readAll(confus_ind),  paste0(rds_dir, 'confusion_index_1_2.rds'))
  rm(confus_ind)
  }
  endCluster()
  
  setTxtProgressBar(pb, 100)
  close(pb)
  message(paste0("DSMART outputs can be located at: ", getwd()))
}
