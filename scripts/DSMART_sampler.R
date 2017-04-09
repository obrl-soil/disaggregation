# this is just the sampling loop from the main dsmart function, its handy to run your input
# polygons through it beforehand to check for errors. If it fails out, its probably one of a few
# things:
#  * wrong data types in columns
#  * wrong column names
#  * a row where n classes doesn't match n percs

DSMART_sampler <- function(indata = NULL, pid_field = NULL, obsdat = NULL) {

  nclass <- length(names(indata@data)[grep('CLASS', names(indata@data))])
  crs    <- indata@proj4string
  sample_points <- vector("list", length = length(indata))

  # count classes per polygon
  indata$cl_count <- apply(indata@data, 1, function(x) {
    sum(!is.na(x) & grepl('CLASS', names(indata@data)) == T)
  })
      
  # sampling loop (one polygon at a time)
  for (i in 1:length(indata@polygons)) { 
    
    polyid  <- indata@data[i, c(pid_field)] 
    area    <- as.integer(indata@polygons[[i]]@area)
    poly_cl <- as.integer(indata@data[i, c('cl_count')])
    
    # samples per sq km
    minrate <- poly_cl * 10
    nsamp   <- ceiling(area / (1000000 / minrate))
    nsamp   <- max(minrate, nsamp)
    
    # sample dat polygon
    spoints <- spsample(indata[i, ], n = nsamp, type = "random", iter = 10)
    
    # get proportions for assigning classes to spoints
    s <- rdirichlet(1, na.omit(unlist(indata@data[i, c(grep('PERC', names(indata@data)))])))
    
    # assign classes to spoints
    polyclassnames <- as.vector(na.omit(unlist(indata@data[i, c(grep('CLASS', names(indata@data)))])))
    classdata <- sample(polyclassnames,
                        size = length(spoints), 
                        replace = TRUE, 
                        prob = s[1, ])
    
    data <- data.frame("POLY_NO" = polyid,
                       "SAMP_NO" = 1:length(spoints),
                       "SAMP_X"  = spoints@coords[, 1],
                       "SAMP_Y"  = spoints@coords[, 2],
                       "CLASS"   = classdata, stringsAsFactors = FALSE)
    
    spointsdf <- SpatialPointsDataFrame(spoints, data, proj4string = crs)
    sample_points[[i]] <- spointsdf
  }

# get all the sampling data for all the polygons for this realisation into one spdf
all_samplepoints <- do.call('rbind', sample_points)

if (!is.null(obsdat)) {
    od         <- obsdat
    names(od)  <- c("POLY_NO", "SAMP_NO", "SAMP_X", "SAMP_Y", "CLASS")
    od$SAMP_NO <- as.numeric(paste0(99, od$SAMP_NO))
    od         <- SpatialPointsDataFrame(od[, c('SAMP_X', 'SAMP_Y')],  
                                         data = od,
                                         proj4string = crs)
    all_samplepoints <- rbind(all_samplepoints, od)
    rm(od)
  }

return(all_samplepoints)

}
