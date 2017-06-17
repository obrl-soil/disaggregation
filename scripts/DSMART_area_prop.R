####################################################################################################
#                                                                                                  #
#            DSMART with area-proportional sampling... and a bunch of other stuff                  #
#                                                                                                  #
#  This version of the DSMART function requires inputs to be projected in a CRS that measures      #
#  distance in units of metres. Equal area projections are preferred for maximum accuracy, but     # 
#  UTM can be used across small N/S extents without too much distortion. The function also         #
#  takes in input data in a different format to the official release; see below for details.       #
#                                                                                                  #
#  Usage                                                                                           #
#  DSMART_AP(covariates = NULL, indata = NULL, pid_field = NULL,  obsdat = NULL, reals = NULL,     #
#            t_factor = 10000, cpus = 1, write_files = FALSE)                                      #
#                                                                                                  #
#  Arguments                                                                                       #
#  covariates   RasterStack; A stack of grids or rasters that are surrogates for environmental     #
#               variables to both guide the C5 model fitting and provide the necessary             #
#               environmental data for the spatial mapping of predicted soil classes.              #
#                                                                                                  #
#  indata       SpatialPolygonsDataFrame; with wide-formatted polygon composition data. Format:    #
#               * one row per polygon                                                              #
#               * A numeric unique ID field for polygons                                           #
#               * A character field for map_code                                                   #
#               * 1-n character fields for soil classes named CLASS1 to CLASSn                     #
#               * 1-n numeric fields for soil class percentages named PERC1 to PERCn. PERC1 must   #
#                 relate to CLASS1, etc.                                                           #
#               Other columns may exist in the object; they will be ignored. NA values in CLASS    #
#               and PERC columns are allowed, e.g. in cases where a polygon only contains 1 of n   #
#               allowable classes. NB: best to readOGR with stringsAsFactors = FALSE.              #
#                                                                                                  #
#  pid_field    String; the name of the numeric ID field in indata.                                #
#                                                                                                  #
#  obsdat       Data.frame; Optional extra. Sample points denoting known soil types at particular  #
#               locations. str(obsdat) should be:                                                  #
#                   'data.frame':  n obs of 5 variables:                                           #
#                    $ POLY_NO: int <the polygon in which the site falls - value of pid_field>     #
#                    $ SAMP_NO: int <the site_ID>                                                  #
#                    $ SAMP_X : num or int <X coordinate of site>                                  #
#                    $ SAMP_Y : num or int <Y coordinate of site>                                  #
#                    $ CLASS  : chr <soil name, must occur in indata>                              #
#                                                                                                  #
#  reals        Integer; number of C5 modeling fitting and mapping realisations to implement.      #
#                                                                                                  #
#  t_factor     Integer; a modifier to clamp extreme values from rdirichlet. Default is 10000.     #
#                                                                                                  #
#  cpus         Integer; number of compute nodes to use. Default is 1.                             #
#                                                                                                  #
#  write_files  Logical; default FALSE. Model outputs are saved as rds objects by default. Set     #
#               this parameter to TRUE to write the c5 model as text, the sampled points and       #
#               and extracted covariates as a point shapefile, and the predicted map as GeoTiff.   #
#               Lookup values are embedded in the geotiff metadata.                                #
#                                                                                                  #
####################################################################################################

DSMART_AP <- function (covariates = NULL, indata = NULL, pid_field = NULL,
                       obsdat = NULL, reals = NULL, cpus = 1, write_files = FALSE) 
{
  
  dir.create('dsmartOuts/',          showWarnings = FALSE)
  dir.create('dsmartOuts/samples',   showWarnings = FALSE)
  dir.create('dsmartOuts/rasters',   showWarnings = FALSE)
  dir.create('dsmartOuts/models',    showWarnings = FALSE)
  
  strd   <- file.path(getwd(), 'dsmartOuts', 'samples')
  strr   <- file.path(getwd(), 'dsmartOuts', 'rasters')
  strm   <- file.path(getwd(), 'dsmartOuts', 'models')
  crs    <- indata@proj4string
  
  # used later to set consistent factoring across model runs
  all_classes  <- sort(na.omit(unique(unlist(indata@data[, grep('CLASS', names(indata@data))]))))
  class_levels <- as.factor(all_classes)
  
  # count classes per polygon
  indata$cl_count <- apply(indata@data, MARGIN = 1, FUN = function(x) {
    sum(!is.na(x) & grepl('CLASS', names(indata@data)) == T)
  })
  
  pb <- txtProgressBar(min = 0, max = reals, style = 3)
  
  for (j in 1:reals) {
    
    beginCluster(cpus)
    sample_points <- vector('list', length = length(indata))
    
    # sampling loop (one polygon at a time)
    for (i in 1:length(indata@polygons)) { 
      
      polyid  <- indata@data[i, c(pid_field)] 
      area    <- as.integer(indata@polygons[[i]]@area)
      poly_cl <- as.integer(indata@data[i, c('cl_count')])
      
      # rate - samples per sq km
      minrate <- poly_cl * 10
      nsamp   <- ceiling(area / (1000000 / minrate))
      nsamp   <- max(minrate, nsamp)
      
      # NB spsample docs say non-spherical coordinates are required, but it doesn't seem to
      # matter with type = 'random'
      spoints <- spsample(indata[i, ], n = nsamp, type = 'random', iter = 10)
      
      # get proportions for assigning classes to spoints
      s <- rdirichlet(1, na.omit(unlist(indata@data[i, c(grep('PERC', names(indata@data)))])) * t_factor)
      
      # NB setting size to = length(spoints) below instead of nsamp handles cases where spsample
      # fails to place nsamp points - usually only an issue if you switch from random to 
      # other sampling types, and possibly with very small/irregular polygons
      polyclassnames <- as.vector(na.omit(unlist(indata@data[i, c(grep('CLASS', names(indata@data)))])))
      classdata <- sample(polyclassnames,
                          size     = length(spoints),
                          replace  = TRUE,
                          prob     = s[1, ])
      
      data <- data.frame('POLY_NO' = polyid,
                         'SAMP_NO' = 1:length(spoints),
                         'SAMP_X'  = spoints@coords[, 1],
                         'SAMP_Y'  = spoints@coords[, 2],
                         'CLASS'   = classdata, stringsAsFactors = FALSE)
      
      spointsdf <- SpatialPointsDataFrame(spoints, data, proj4string = crs)
      sample_points[[i]] <- spointsdf
    }
    
    # get all the sampling data for all the polygons for this realisation into one spdf
    all_samplepoints <- do.call('rbind', sample_points)
    rm(sample_points)
    
    # if aux points are in use, add them to all_samplepoints
    # NB aux points now have to be supplied with their intersecting poly ID and their site_id;
    # together these may be non-unique with the sample point IDs generated above. Adding a 99 out
    # front of the aux SAMP_NOs to handle that for now    
    if (!is.null(obsdat)) {
      od         <- obsdat
      names(od)  <- c('POLY_NO', 'SAMP_NO', 'SAMP_X', 'SAMP_Y', 'CLASS')
      od$SAMP_NO <- as.numeric(paste0(99, od$SAMP_NO))
      od         <- SpatialPointsDataFrame(od[, c('SAMP_X', 'SAMP_Y')], 
                                           data = od,
                                           proj4string = crs)
      all_samplepoints <- rbind(all_samplepoints, od)
      rm(od)
    }
    
    # sample covariates, appending the sampled data to the SPDF
    # NB with GSIF::extract, sampling from non-stacked rasters is possible...hmmmm
    all_samplepoints <- raster::extract(covariates, all_samplepoints, sp = TRUE)
    
    # forces all outputs to be on the same scale eg. raster value 1 always equals factor level 1
    all_samplepoints$CLASS         <- as.factor(all_samplepoints$CLASS) 
    levels(all_samplepoints$CLASS) <- class_levels
    
    # ditch any points with NA values in the extracted covariate data, C5.0 does not approve
    all_samplepoints <- all_samplepoints[complete.cases(all_samplepoints@data), ]
    
    # generate decision tree 
    res <- C5.0(all_samplepoints@data[ , c(6:ncol(all_samplepoints@data))],
                y = all_samplepoints@data$CLASS)
    
    # Generate lookup table to match GeoTIFF values
    lookup <- data.frame('ID' = as.integer(as.factor(res$levels)), 'CLASS' = as.factor(res$levels)) 
    
    # make prediction map 
    # The resulting temp file for each realisation will require (2 x ncells) bytes of 
    # storage space in temp folder (C:/Users/username/AppData/Local/Temp by default on Win) 
    r1 <- clusterR(na.omit(covariates), raster::predict, args = list(res), datatype = 'INT2S')
    
    # factorise r1
    levels(r1) <- lookup
    
    # save for later (readALL or your raster rds' will just be pointers to temp files)
    saveRDS(all_samplepoints, file.path(strd, paste0('samples_', j, '.rds')))
    saveRDS(res, file.path(strm, paste0('C5_model_', j, '.rds')))
    suppressWarnings(saveRDS(readAll(r1), file.path(strr, paste0('map_', j, '.rds'))))
    
    # make it optional to write rasters and shapefiles etc
    if (write_files == TRUE) {
      
      # decision tree to plain text tho lets bh the rds is more useful
      out <- capture.output(summary(res))
      f2  <- file.path(strm, paste0('C5_model_', j, '.csv'))
      cat(out, file = f2, sep = '\n', append = TRUE)
      
      # write probability map from this realisation to GeoTIFF (lookup values are embedded)
      nme <- paste0(strr, 'map_', j, '.tif')
      writeRaster(r1, filename = nme, format = 'GTiff', overwrite = TRUE, datatype = 'INT2S', 
                  NAflag = -9999)
      
      # may as well keep the class lookup in txt as can't use it in QGIS
        if (j == 1) {
          write.table(lookup, file.path(strr, 'class_lookup.txt'), 
                      col.names = TRUE, row.names = FALSE,
                      quote = FALSE, sep = ',')
        }
            
      # make a lookup table for all_samplepoints covariate column names, because they're about
      # to get severely abbreviated for writing to shp
      cov_LUT_nm <- file.path(strd, 'covariate_LUT.csv')
      cov_names <- names(all_samplepoints[6:ncol(all_samplepoints)])
      cov_shpnames <- paste0('COV_', 1:length(cov_names))
      if (!file.exists(cov_LUT_nm)) {
        cov_LUT <- data.frame('COV_NAMES' = cov_names, 'SHPCOL_NAMES' = cov_shpnames)
        write.table(cov_LUT, file = cov_LUT_nm, 
                    sep = ', ', quote = FALSE, col.names = TRUE, row.names = FALSE)
      }
      
      # abbreviate colnames and write sample data from this realisation to a shapefile
      names(all_samplepoints)[6:ncol(all_samplepoints)] <- cov_shpnames
      spname <- paste0('samplepoints_', j)
      writeOGR(all_samplepoints,
               dsn = file.path(getwd(), 'dsmartOuts', 'samples'),
               layer = spname,
               driver = 'ESRI Shapefile',
               overwrite = TRUE)
      
    }
    # tidy tidy
    rm(r1)
    endCluster()
    setTxtProgressBar(pb, j)
  }
  
  close(pb)
  message(paste0('DSMART outputs can be located at: ', getwd(), '/dsmartOuts/'))
}
