## from https://github.com/NERC-CEH/IDM_comparisons and Suhaimi et al 2021 

# PACKAGES --------------------------------------------------------------

# packages
library(INLA)
library(reshape2)
library(rgeos)
library(fields)
library(deldir)
library(dplyr)
library(RColorBrewer)
library(data.table)
library(raster)
library(tidyverse)

## load original data from FaCeT project? default is FALSE
if (!exists("get_data")) { get_data = FALSE }
## use small spatial subset of the data? default is TRUE
if (!exists("sp_subset")) { sp_subset = TRUE }
## want to simulate some data? this is from Suhaimi and is not tested but retained for posterity
if (!exists("sim_data")) { sim_data = FALSE }
## want to clean up all the random objects that get created? usually this will be true
if (!exists("clean")) { clean = TRUE}
## rebuild bias fields?
if (!exists("build_bias")) { build_bias = FALSE}
## make plots as the data loads and other intermediate outputs are generated?
if (!exists("make_plot")) { make_plot = TRUE}
if (!exists("build_pred")) { build_pred = FALSE}
if (!exists("run_cameletti")) { run_cameletti = FALSE}
if (!exists("run_suhaimi")) { run_suhaimi = FALSE}
if (!exists("coarse_mesh")) { coarse_mesh = TRUE}
if (!exists("scale_rasters")) { scale_rasters = FALSE}


# GENERATE DATA FOR SIMS --------------------------------------------------------------

## we have real data so dont need this but leaving here as a guide for data formats, etc
## this section will run if needed to look at these variables from Suhaimi
if (sim_data){
  ## Set parameters
  
  #genData
  dim = c(100,100)   
  lambda = -3    
  env.beta = 1.2 
  plotdat = FALSE
  sigma2x = 1     # var in rLGCP
  kappa = 0.05    # =1/scale in rLGCP  [0.02, 0.1] good range
  
  #genStrataLam
  strata = 25
  rows = 5
  cols = 5
  plot = TRUE
  
  #addSpatialBias - detection
  # this is where we can control bias
  probs = rep(c(0.5, 0.3, 0.1, 0.05, 0.01),5) #default
  #probs = rep(c(0.2, 0.2, 0.2, 0.2, 0.2),5) # unbiased
  #probs = rep(c(0.5, 0.4, 0.1, 0.01, 0.001),5) # very biased
  
  #sampleStructured
  nsamp = 100  
  #nsamp = 50   #small
  #nsamp = 200  #large
  plotdat = TRUE
  qsize = 1          #~~~ neighborhood; buffer <- (qsize-1)/2
  
  #~~~ mesh
  mesh.edge = c(7, 14)   
  mesh.offset = c(2, 7)   
  
  #~~~ dimension
  resolution = c(5,5)
  
  
  ##read seeds
  
  #seed_all <- read.csv("seed_biased.csv")
  seed_new <- 399
  
  # Generate data ##
  
  # change the seed for each run
  
  i <- as.numeric(commandArgs(trailingOnly = TRUE))[1]
  
  #seed_new <- seed_all[i,2]
  
  print(seed_new)
  
  #start_time_ori <- Sys.time()
  
  start_time <- Sys.time()
  seed = seed_new
  
  source("Functions to generate data and sample.R")
  g1 <- genDataFunctions(dim = dim,
                         lambda = lambda,
                         env.beta = env.beta,
                         seed = seed_new,
                         kappa = kappa,
                         sigma2x = sigma2x,
                         strata = strata,
                         rows = rows,
                         cols = cols,
                         probs = probs,
                         nsamp = nsamp,
                         plot = F,
                         plotdat = F,
                         qsize = 1)
  
  structured_data <- g1$structured_data
  unstructured_data <- g1$unstructured_data
  biasfield <- g1$biasfield
  dat1 <- g1$dat1
  biascov <- g1$biascov
  strata1 <- g1$strata1
  
  print(seed_new)
  print(nrow(unstructured_data))
  
  
  #' Visualise thinned unstructured data
  #+ echo = FALSE
  #par(mfrow=c(1,1))
  image.plot(list(x=dat1$Lam$xcol, y=dat1$Lam$yrow, z=t(dat1$rf.s)),
             main='Thinned unstructured data', asp=1,
             col = hcl.colors(12, "RdBu", rev = TRUE))
  points(unstructured_data$x, unstructured_data$y, pch = 20)#note rescale again - plotting back on original
  
  #' Visualise structured data
  #+ echo = FALSE
  #par(mfrow=c(1,1))
  image.plot(list(x=dat1$Lam$xcol, y=dat1$Lam$yrow, z=t(dat1$rf.s)),
             main='Structured data', asp=1,
             col = hcl.colors(12, "RdBu", rev = TRUE))
  points(structured_data$x,structured_data$y, pch = 21, bg = structured_data$presence, col = "black")
  #par(xpd = TRUE)
  #legend(0,115,c("Absence", "Presence"), pch = 21, col = "black", pt.bg = c(0,1))
  
  ggplot() +
    geom_raster(data = dat1$gridcov , aes(x = x, y = y, fill = sim1))# + 
  
}

# LOAD REAL DATA --------------------------------------------------------------

if (get_data){
  ## observer data "enhanced" with env covariates, incl all data that passed QC
  data_obs <- data.table::fread('../NASA-FaCeT/data/enhance/observer/observer-enhanced.csv')
  data_obs$pres <- ifelse(data_obs$CATCH > 0, 1, 0)
  data_obs %>% group_by(pres) %>% summarise(n=n())
  data.table::fwrite(data_obs, '~/Google Drive/Shared drives/MPG_WHOI/data/bsh_inla/bsh_observer.csv')
  
  ## mark-recapture data "enhanced" with env covariates, incl all data that passed QC
  data_marker <- data.table::fread('../NASA-FaCeT/data/enhance/iccat/iccat-enhanced.csv')
  data_marker <- data_marker %>% filter(SpeciesCode == 'BSH' & pres == 1) # approx 25k combined release/recoveries
  data_marker %>% group_by(ObsType) %>% summarise(n=n())
  data.table::fwrite(data_marker, '~/Google Drive/Shared drives/MPG_WHOI/data/bsh_inla/bsh_marker.csv')
  
  ## electronic data "enhanced" with env covariates, incl all data that passed QC. 
  ## spot data has been fit with movement model (foieGras) and standardized to daily locations (with error from model)
  data_spot <- data.table::fread('../NASA-FaCeT/data/enhance/etag/enhanced/etag-enhanced-bg_sat_extent2.csv')
  data_spot <- data_spot %>% filter(pres == 1) # 
  data_spot$tag_type <- 'spot'
  #data_spot %>% group_by(instrument_name) %>% summarise(n=n())
  
  ## psat data has been fit with movement model (GPE3) and standardized to daily locations (with error from model)
  data_psat <- data.table::fread('../NASA-FaCeT/data/enhance/etag/enhanced/etag-enhanced-bg_psat_extent2.csv')
  data_psat <- data_psat %>% filter(pres == 1) # 
  data_psat$tag_type <- 'psat'
  #data_psat %>% group_by(instrument_name) %>% summarise(n=n())
  ## combine etag data
  data_etag <- rbind(data_spot, data_psat)
  rm(data_psat); rm(data_spot)
  data.table::fwrite(data_etag, '~/Google Drive/Shared drives/MPG_WHOI/data/bsh_inla/bsh_etag.csv')
} else{
  
  data_observer <- data.table::fread('~/Google Drive/Shared drives/MPG_WHOI/data/bsh_inla/bsh_observer.csv')
  data_observer$hookhours <- data_observer$NUMBER_HOOKS_SET / data_observer$SOAK_TIME
  data_marker <- data.table::fread('~/Google Drive/Shared drives/MPG_WHOI/data/bsh_inla/bsh_marker.csv')
  data_etag <- data.table::fread('~/Google Drive/Shared drives/MPG_WHOI/data/bsh_inla/bsh_etag.csv')
  
}

if (sp_subset){
  
  print(paste("**---- Using small subset of study domain to speed things up ----**"))
  
  xl <- c(-75, -65)
  yl <- c(35, 45)
  
  data_observer <- data_observer %>% filter(lon > xl[1] & lon < xl[2] & lat > yl[1] & lat < yl[2])
  data_marker <- data_marker %>% filter(lon > xl[1] & lon < xl[2] & lat > yl[1] & lat < yl[2])
  data_etag <- data_etag %>% filter(lon > xl[1] & lon < xl[2] & lat > yl[1] & lat < yl[2])
  
} else{
  
  print(paste("**---- Using full N Atlantic study domain ----**"))
  
  xl <- c(-100, -5)
  yl <- c(10, 55)
  
}

data_observer <- data_observer %>%
  filter_at(vars(sst:rugosity), all_vars(!is.na(.)))

data_marker <- data_marker %>%
  filter_at(vars(sst:rugosity), all_vars(!is.na(.)))

data_etag <- data_etag %>%
  filter_at(vars(sst:rugosity), all_vars(!is.na(.)))


# BUILD BIAS FIELDS --------------------------------------------------------------

#===
## OBSERVER

if (build_bias){
  
  ## observer will be distance from bathy contour (i.e. representing shelf break as focal point of effort)
  ## no. hooks and soak time -> would make much more difference if worried about no. of individuals captured
  r_bias <- raster(xmn=xl[1], ymn=yl[1], xmx=xl[2], ymx=yl[2], res=.5)
  observer_bias_setcount <- raster::rasterize(cbind(data_observer$lon, data_observer$lat), r_bias, fun='count')
  observer_bias_setcount_normalized <- observer_bias_setcount / raster::cellStats(observer_bias_setcount, 'max')
  
  ## adjust to "hook-hours" which is represented in "EFFORT" column
  observer_bias_hookhours <- raster::rasterize(cbind(data_observer$lon, data_observer$lat), field=data_observer$EFFORT, r_bias, fun='sum')
  observer_bias_hookhours_normalized <- observer_bias_hookhours / raster::cellStats(observer_bias_hookhours, 'max')
  
  par(mfrow=c(2,1))
  plot(observer_bias_setcount_normalized, main='observer set count (normalized)')
  world(add=T)
  plot(observer_bias_hookhours_normalized, main='observer hook-hours (normalized)')
  world(add=T)
  
  observer_bias <- observer_bias_hookhours_normalized
  raster::writeRaster(observer_bias, '~/Google Drive/Shared drives/MPG_WHOI/data/bsh_inla/observer_bias.grd')
  
  ## housekeeping
  if (clean){
    rm(observer_bias_hookhours_normalized); rm(observer_bias_setcount_normalized); gc()
  }
} else{
  observer_bias <- raster::raster('~/Google Drive/Shared drives/MPG_WHOI/data/bsh_inla/observer_bias.grd')
  if (make_plot) plot(observer_bias, main='observer bias'); world(add=T)
}

#===
## MARKER

if (build_bias){
  ## load marker data for all species
  data_marker_allsp <- data.table::fread('../NASA-FaCeT/data/enhance/iccat/iccat-enhanced.csv')
  data_marker_allsp <- data_marker_allsp %>% filter(pres == 1) 
  marker_bias <- raster::rasterize(cbind(data_marker_allsp$lon, data_marker_allsp$lat), r_bias, field=data_marker_allsp$pres, fun=function(x,...) sum(x, na.rm=T) )
  par(mfrow=c(1,1)); plot(marker_bias, main='all species marker tag "effort"'); world(add=T)
  raster::writeRaster(marker_bias, '~/Google Drive/Shared drives/MPG_WHOI/data/bsh_inla/marker_bias.grd')
  
  ## housekeeping
  if (clean){
    rm(data_marker_allsp); gc()
  }
} else{
  marker_bias <- raster::raster('~/Google Drive/Shared drives/MPG_WHOI/data/bsh_inla/marker_bias.grd')
  if (make_plot) plot(marker_bias, main='marker bias'); world(add=T)
}

#===
## ETAG
## etag will be daily distance from previous location? or do we need this?
## sum across all tag UDs to derive bias field

if (build_bias){
  
  calc.errEll_inla <- function(locs, locs.grid){
    ## borrowed and modified from HMMoce  
    # set up a larger grid to base ellipse on and to shift that error, if necessary (GPE only)
    ngrid <- rev(dim(locs.grid$lon))
    lon1 <- seq(min(locs.grid$lon[1,]) - 10, max(locs.grid$lon[1,]) + 10, by = locs.grid$dlo)
    lat1 <- seq(min(locs.grid$lat[,1]) - 10, max(locs.grid$lat[,1]) + 10, by = locs.grid$dla)
    g1 <- meshgrid(lon1, lat1)
    
    
    slon.sd <- locs$longitudeError #semi minor axis
    slat.sd <- locs$latitudeError #semi minor axis
    
    # calc semi minor axis based on longitude error
    #slon.sd <- locs$Error.Semi.minor.axis / 1000 / 111 #semi minor axis
    L.light.lon <- stats::dnorm(t(g1$X), locs$lon, slon.sd) # Longitude data
    #slat.sd <- locs$Error.Semi.major.axis / 1000 / 111 #semi major axis
    L.light.lat <- stats::dnorm(t(g1$Y), locs$lat, slat.sd)
    
    #image.plot(g$lon[1,],g$lat[,1],L.light.lat*L.light.lon)
    
    # create the ellipse by multiplying lat * lon error
    L <- raster::flip(raster::raster(t(L.light.lat * L.light.lon), xmn = min(lon1), 
                                     xmx = max(lon1), ymn = min(lat1), ymx = max(lat1)), direction = 'y')
    
    return(L)
  }
  source('../HMMoce/R/repmat.r'); source('../HMMoce/R/meshgrid.r')
  source('../HMMoce/R/setup.locs.grid.r')
  
  bathy <- raster::raster('~/Google Drive/Shared drives/MPG_WHOI/env_data/bathy/global_bathy_0.01.nc')
  bathy <- raster::rotate(bathy)
  bathy <- raster::crop(bathy, raster::extent(xl[1], xl[2], yl[1], yl[2]))
  
  locs.grid <- setup.locs.grid(list(extent(bathy)@xmin,
                                    extent(bathy)@xmax,
                                    extent(bathy)@ymin,
                                    extent(bathy)@ymax), res = 'hycom')
  
  data_etag_mod <- data_etag #%>% dplyr::select(date, lon, lat, longitudeError, latitudeError)
  
  ncores <- ceiling(parallel::detectCores() * .9)
  cl = parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl, cores = ncores)
  t1 <- Sys.time()
  foreach::foreach(i = 1:nrow(data_etag_mod)) %dopar%{
    if (data_etag_mod$longitudeError[i] < 0.005 | is.na(data_etag_mod$longitudeError[i])) data_etag_mod$longitudeError[i] <- 0.005
    if (data_etag_mod$latitudeError[i] < 0.005 | is.na(data_etag_mod$latitudeError[i])) data_etag_mod$latitudeError[i] <- 0.005
    
    if (nchar(i) == 1) i_paste <- paste0('0000', i)
    if (nchar(i) == 2) i_paste <- paste0('000', i)
    if (nchar(i) == 3) i_paste <- paste0('00', i)
    if (nchar(i) == 4) i_paste <- paste0('0', i)
    if (nchar(i) == 5) i_paste <- i
    
    etag_error <- calc.errEll_inla(data_etag_mod[i,], locs.grid)
    etag_error <- etag_error / raster::cellStats(etag_error, 'sum')
    raster::writeRaster(etag_error, filename = paste0('./etag_error2/etag_error_', i_paste, '.grd'), overwrite=T)
  }
  t2 <- Sys.time()
  
  
  rm_idx <- which(data_etag_mod$longitudeError > 5 & data_etag_mod$latitudeError > 5)
  fList <- list.files('./etag_error2/', full.names=TRUE)
  fList <- fList[grep('.grd', fList)]
  fList <- fList[-rm_idx]
  for (i in 1:length(fList)){
    if ((i - 1) %% 100 == 0 | i - 1 == 0){
      if (i == 1) {
        all_etag_error <- sum(raster::stack(fList[i:(i + 99)]), na.rm=T)
      } else if ((i + 99) > length(fList)){
        all_etag_error <- sum(raster::stack(all_etag_error,
                                            sum(raster::stack(fList[i:length(fList)]), na.rm=T)), na.rm=T)
      } else{
        all_etag_error <- sum(raster::stack(all_etag_error,
                                            sum(raster::stack(fList[i:(i + 99)]), na.rm=T)), na.rm=T)
      }
      gc() 
      print(i)
    }
  }
  plot(all_etag_error)
  raster::writeRaster(all_etag_error, '~/Google Drive/Shared drives/MPG_WHOI/data/bsh_inla/etag_bias_FULLres.grd')
  
  bathy_resamp <- raster::resample(bathy, all_etag_error)
  error_try <- raster::mask(all_etag_error, bathy_resamp)
  error_try <- raster::resample(error_try, r_bias)
  #error_try <- raster::reclassify(error_try, c(c(0, 0, NA)))
  error_try <- error_try / raster::cellStats(error_try, 'max')
  
  if (make_plot){
    #pdf('etag_bias.pdf', height=10, width=12)
    par(mfrow=c(2,1))
    #plot(raster::mask(all_etag_error, bathy_resamp))
    plot(error_try / raster::cellStats(error_try, 'max'), main='sum across normalized UDs')
    world(add=T)
    #points(data_etag_mod$lon, data_etag_mod$lat)
    r_bias <- raster(xmn=xl[1], ymn=yl[1], xmx=xl[2], ymx=yl[2], res=.5)
    etag_bias_pts <- raster::rasterize(cbind(data_etag_mod$lon, data_etag_mod$lat), r_bias, fun='count')
    plot(etag_bias_pts / raster::cellStats(etag_bias_pts, 'max'), main='rasterize summary of locations (count)'); world(add=T)
    #dev.off()
  }
  
  raster::writeRaster(error_try, '~/Google Drive/Shared drives/MPG_WHOI/data/bsh_inla/etag_bias.grd')
  ## housekeeping
  if (clean){
    rm(bathy_resamp); rm(error_try); rm(all_etag_error); rm(fList); rm(etag_error); rm(locs.grid); rm(data_etag_mod); gc()
  }
} else{
  etag_bias <- raster::raster('~/Google Drive/Shared drives/MPG_WHOI/data/bsh_inla/etag_bias.grd')
  if (make_plot) plot(etag_bias, main='etag bias'); world(add=T)
}


# BUILD MESH --------------------------------------------------------------
## some guidance here: https://punama.github.io/BDI_INLA/
## and a book on spatiotemporal models with INLA / R: https://sites.google.com/a/r-inla.org/stbook/

## construct the mesh 
proj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
atl <- rgdal::readOGR('~/Google Drive/Shared drives/MPG_WHOI/env_data/shapefiles/atlantic.shp')
atl <- raster::crop(atl, raster::extent(xl[1], xl[2], yl[1], yl[2]))

bdry <- inla.sp2segment(atl)
bdry$loc <- inla.mesh.map(bdry$loc)
if (coarse_mesh){
  mesh <- inla.mesh.2d(boundary = bdry, max.edge=c(2,4), 
                       offset = c(0.5, 1),
                       cutoff = 0.3)
  
} else{
  stop('Fine mesh has not been defined.')
}

mesh$crs <- proj

if (make_plot){
  plot(mesh)
  world(add=T)
}

# MAKE SPDE & MATRICES --------------------------------------------------------------

# make SPDE 
spde <- inla.spde2.matern(mesh = mesh, alpha = 2)

# Make observation structure for estimation data
# Create A matrices (mapping between the mesh nodes and station locations) #

# make A matrix for structured data
data_observer_A <- inla.spde.make.A(mesh = mesh,
                                    loc = as.matrix(data_observer[,c('lon','lat')]))

# make A matrices for unstructured data
data_marker_A <- inla.spde.make.A(mesh = mesh,
                                  loc = as.matrix(data_marker[,c('lon','lat')]))

data_etag_A <- inla.spde.make.A(mesh = mesh,
                                loc = as.matrix(data_etag[,c('lon','lat')])
                                #group=Piemonte_data$time, ## these temporal groupings are from Cameletti, not sure if we need them?
                                #n.group=n_days
                                )
## ^^^ perhaps grouping argument above is way to account for individual in etag data?


# create integration stack
max_x <- max(mesh$loc[,1])
max_y <- max(mesh$loc[,2])
loc.d <- t(matrix(c(0, 0, max_x, 0, max_x, max_y, 0, max_y, 0, 0), 2))

#make dual mesh
dd <- deldir::deldir(mesh$loc[, 1], mesh$loc[, 2])
tiles <- deldir::tile.list(dd)

#make domain into spatial polygon
domainSP <- SpatialPolygons(list(Polygons(
  list(Polygon(loc.d)), '0')))

#intersection between domain and dual mesh
poly.gpc <- as(domainSP@polygons[[1]]@Polygons[[1]]@coords, "gpc.poly")

# w now contains area of voronoi polygons
w <- sapply(tiles, function(p) rgeos::area.poly(rgeos::intersect(as(cbind(p$x, p$y), "gpc.poly"), poly.gpc)))

nv <- mesh$n
n_marker <- nrow(data_marker)
n_etag <- nrow(data_etag)

# change data to include 0s for nodes and 1s for presences. 
# only necessary for "unstructured" data types (i.e. PO)
y.pp_marker <- rep(0:1, c(nv, n_marker)) ## corresponds to y.pp in Suhaimi "unstructured" example
y.pp_etag <- rep(0:1, c(nv, n_etag)) ## corresponds to y.pp in Suhaimi "unstructured" example

# add expectation vector (area for integration points/nodes and 0 for presences)
e.pp_marker <- c(w, rep(0, n_marker))
e.pp_etag <- c(w, rep(0, n_etag))

#diagonal matrix for integration point A matrix
imat <- Diagonal(nv, rep(1, nv))

A.pp_marker <- rbind(imat, data_marker_A)
A.pp_etag <- rbind(imat, data_etag_A)


# SCALE ENV COVARIATES - MARKER TAGS --------------------------------------------------------------

mean_covariates <- apply(data_marker %>% 
                           dplyr::select(sst:rugosity), 2, function(x) mean(x, na.rm=T))
sd_covariates <- apply(data_marker %>% 
                           dplyr::select(sst:rugosity), 2, function(x) sd(x, na.rm=T))
data_marker <- data_marker %>% as.data.frame()
data_marker[,c(which(names(data_marker) == 'sst'):which(names(data_marker) == 'rugosity'))] <-
  scale(data_marker[,c(which(names(data_marker) == 'sst'):which(names(data_marker) == 'rugosity'))],
        mean_covariates, sd_covariates)


# SCALE ENV COVARIATES - OBSERVER AND ETAG --------------------------------------------------------------

## scale using mean/sd from marker tag data

data_observer <- data_observer %>% as.data.frame()
data_observer[,c(which(names(data_observer) == 'sst'):which(names(data_observer) == 'rugosity'))] <-
  scale(data_observer[,c(which(names(data_observer) == 'sst'):which(names(data_observer) == 'rugosity'))],
        mean_covariates, sd_covariates)

data_etag <- data_etag %>% as.data.frame()
data_etag[,c(which(names(data_etag) == 'sst'):which(names(data_etag) == 'rugosity'))] <-
  scale(data_etag[,c(which(names(data_etag) == 'sst'):which(names(data_etag) == 'rugosity'))],
        mean_covariates, sd_covariates)



if (scale_rasters){
  
  for (i in 1993:2019){
    fList <- list.files(paste0('/Volumes/Elements/glorys_monthly/', i, '/'), full.names = TRUE, recursive = TRUE)
    
    if (i == 1993){
      fList_all <- fList
    } else{
      fList_all <- c(fList_all, fList)
    }
    
  }
  fList_all <- fList_all[grep('.grd', fList_all)]
  
  for (i in 1:length(fList_all)){
    r <- raster::stack(fList_all[i])
    if (i == 1){
      mld_stack <- r[['mld']]
    } else{
      mld_stack <- stack(mld_stack, r[['mld']])
    }
  }
  mld_clim <- mean(mld_stack)
  raster::writeRaster(mld_clim, filename = paste0('/Volumes/Elements/glorys_monthly/climMean/cmems_mod_glo_phy_my_0.083_P1M-m_climMean_mld.grd'))
  mld_clim_scaled <- raster::calc(mld_clim, function(x) {scale(x, mean_covariates['mld'], sd_covariates['mld'])})
  raster::writeRaster(mld_clim_scaled, filename = paste0('/Volumes/Elements/glorys_monthly/climMean_allMonths/mld/cmems_mod_glo_phy_my_0.083_P1M-m_climMean_mld.grd'), overwrite = TRUE)
  
  month_vec <- c('01','02','03','04','05','06','07','08','09','10','11','12')
  env_vars <- c('sst','sss','ssh','log_eke','sst_sd','ssh_sd','sss_sd')
  for (tt in env_vars){
    for (i in month_vec){
      ## get monthly climatology raster
      r <- raster::raster(paste0('/Volumes/Elements/glorys_monthly/climMean/', tt, '/cmems_mod_glo_phy_my_0.083_P1M-m_climMean_', tt, '_', i, '.grd'))
      
      ## stack by month
      if (!exists('r_clim')){
        r_clim <- r
      } else{
        r_clim <- raster::stack(r_clim, r)
      }
      
    }
    
    ## take mean, scale by marker data and save
    r_clim <- mean(r_clim)
    r_clim_scaled <- raster::calc(r_clim, function(x) {scale(x, mean_covariates[tt], sd_covariates[tt])})
    if(!dir.exists(paste0('/Volumes/Elements/glorys_monthly/climMean_allMonths/', tt, '/'))) dir.create(paste0('/Volumes/Elements/glorys_monthly/climMean_allMonths/', tt, '/'), recursive = TRUE)
    raster::writeRaster(r_clim_scaled, filename = paste0('/Volumes/Elements/glorys_monthly/climMean_allMonths/', tt, '/cmems_mod_glo_phy_my_0.083_P1M-m_climMean_', tt, '.grd'), overwrite = TRUE)
    rm(r_clim); rm(r)
  }
  
  fList <- list.files('/Volumes/Elements/glorys_monthly/climMean_allMonths/', full.names = TRUE, recursive = TRUE)
  fList <- fList[grep('.grd', fList)]
  par(mfrow=c(3,3))
  for (i in 1:length(fList)){
    r <- raster::raster(fList[i])
    plot(r, main=fList[i])
  }
  
  ## same for static covars
  bathy <- raster::raster('~/Google Drive/Shared drives/MPG_WHOI/env_data/bathy/global_bathy_0.01.nc')
  bathy <- raster::rotate(bathy)
  bathy <- raster::crop(bathy, raster::extent(xl[1], xl[2], yl[1], yl[2]))
  bathy_scaled <- raster::calc(bathy, function(x) {scale(x, mean_covariates['bathy'], sd_covariates['bathy'])})
  raster::writeRaster(bathy_scaled, filename='/Volumes/Elements/bathy_scaled_BSHmarker.grd', overwrite = TRUE)
  
  rugosity <- raster::raster('~/Google Drive/Shared drives/MPG_WHOI/env_data/bathy/bathy_rugosity.nc')
  rugosity <- raster::rotate(rugosity)
  rugosity <- raster::crop(rugosity, raster::extent(xl[1], xl[2], yl[1], yl[2]))
  rugosity_scaled <- raster::calc(rugosity, function(x) {scale(x, mean_covariates['rugosity'], sd_covariates['rugosity'])})
  raster::writeRaster(rugosity_scaled, filename='/Volumes/Elements/rugosity_scaled_BSHmarker.grd', overwrite = TRUE)
  
} else{
  
  env_vars <- c('sst','sss','ssh','mld','log_eke','sst_sd','ssh_sd','sss_sd')
  for (tt in env_vars){
    if (tt == env_vars[1]) try(rm(clim_stack), TRUE)
    r <- raster::raster(paste0('~/Google Drive/Shared drives/MPG_WHOI/env_data/glorys_monthly/climMean_allMonths/', tt, '/cmems_mod_glo_phy_my_0.083_P1M-m_climMean_', tt, '.grd'))
    names(r) <- tt
    if (!exists('clim_stack')){
      clim_stack <- r
    } else{
      clim_stack <- raster::stack(clim_stack, r)
    }
  }
  
  bathy_scaled <- raster::raster('~/Google Drive/Shared drives/MPG_WHOI/env_data/glorys_monthly/climMean_allMonths/bathy_scaled_BSHmarker.grd')
  rugosity_scaled <- raster::raster('~/Google Drive/Shared drives/MPG_WHOI/env_data/glorys_monthly/climMean_allMonths/rugosity_scaled_BSHmarker.grd')
  
}

# EXTRACT COVARS TO INTEGRATION POINTS --------------------------------------------------------------

# get covariate(s) for integration points
env_covariates <- raster::extract(clim_stack, cbind(mesh$loc[,1], mesh$loc[,2])) %>% as.data.frame()
env_covariates$bathy <- raster::extract(bathy_scaled, cbind(mesh$loc[,1], mesh$loc[,2]))
env_covariates$rugosity <- raster::extract(rugosity_scaled, cbind(mesh$loc[,1], mesh$loc[,2]))

## get bias covar for marker data and extract at integration points
data_marker$bias <- raster::extract(marker_bias, cbind(data_marker$lon, data_marker$lat))
data_marker$bias_scaled <- c(scale(data_marker$bias, mean(data_marker$bias, na.rm=T), sd(data_marker$bias, na.rm=T)))
marker_bias_scaled <- raster::calc(marker_bias, function(x) {scale(x, mean(data_marker$bias, na.rm=T), sd(data_marker$bias, na.rm=T))})
biascovariate_marker_re <- raster::extract(marker_bias_scaled, cbind(mesh$loc[,1], mesh$loc[,2]))

## get bias covar for observer data
data_observer$hookhours_scaled <- c(scale(data_observer$hookhours, mean(data_observer$hookhours, na.rm=T), sd(data_observer$hookhours, na.rm=T)))

# CHECK ENV COVARIATES (COLLINEARITY, ETC) --------------------------------------------------------------

## NA values should have already been removed from all env covariates
# MCA: The issue of collinearity doesn't seem to come up in the lit on iSDMs but given that coefficient estimation
# is done with joint likelihood even in the 'correlative' model structure we are using, I believe we need to 
# check for collinearity in each dataset used in the model. My thought here is that if collinearity is present 
# but not addressed in a given dataset, then the resulting coefficient likelihood from that submodel will skew 
# the overall joint likelihood estimate across submodels. Below, I provide a function I wrote for quickly 
# identifying collinear pairs of covariates in a provided dataset. This is designed to identify those pairs with 
# |r| > 0.7 (as recommended by Dormann et al. 2013 Ecography) without having to stare at a bunch of plots or match
# rows to columns visually in a large X by X matrix while remembering to ignore the diagonals and duplication across
# the upper or lower triangles of the correlation matrix.

collinearity <- function(covariates) {
  # input example:   collinearity(covariates = data_observer[ , c("sst","ssh","chl","bathy")]) # just call covariates by their col names and the function does the rest
  ## function written by MCA
  
  ## if NA values present, can give spurious results
  if (nrow(na.omit(covariates)) < nrow(covariates)){
    warning('NA values present in covariates. Using na.omit() to remove them.')
    print(paste0('Input covariates contained ', nrow(covariates), ' rows of data.'))
    covariates <- na.omit(covariates)
    print(paste0('Output covariates, with NA removed, contained ', nrow(covariates), ' rows of data.'))
  }
  
  correlations <- cor(covariates)
  dimlength <- nrow(covariates)
  diags <- seq(1, dimlength ^ 2, by = (dimlength + 1))
  colinears <- which(abs(correlations) > 0.7 & upper.tri(correlations, diag = FALSE) == TRUE)
  
  if (length(colinears) != 0){
    for (i in 1:length(colinears)){
      ind <- which(correlations == correlations[colinears[i]] & upper.tri(correlations, diag = FALSE) == TRUE, arr.ind = TRUE)
      print(paste(rownames(correlations)[ind[1]], colnames(correlations)[ind[2]], sep = ", "))
    }
  } else {print("No pairwise comparisons with |r| > 0.70")}
}

collinearity(na.omit(data_observer %>% dplyr::select(sst:rugosity)))
collinearity(na.omit(data_marker %>% dplyr::select(sst:rugosity)))
collinearity(na.omit(data_etag %>% dplyr::select(sst:rugosity)))

# INLA STACK & FIT - BSH --------------------------------------------------------------
## build INLA stacks

stack_observer <- inla.stack(
  data = list(pres = cbind(data_observer$pres, NA, NA),
              Ntrials = rep(1, nrow(data_observer))), ## still not clear what exactly Ntrials is but it seems important in "structured" data types (i.e. PA)
  effects = list(
    list( ## element 1 of effects list contains intercept, env covars and bias covar (if any)
      data.frame(
        intercept_observer = rep(1, length(data_observer$pres)), 
        sst = data_observer$sst,
#        sss = data_observer$sss,
#        ssh = data_observer$ssh,
        mld = data_observer$mld,
        log_eke = data_observer$log_eke,
        sst_sd = data_observer$sst_sd,
#        sss_sd = data_observer$sss_sd,
#        ssh_sd = data_observer$ssh_sd,
        bathy = data_observer$bathy,
        rugosity = data_observer$rugosity,
        bias_observer = data_observer$hookhours_scaled
      )),
    list(data.frame(spatial_field = 1:spde$n.spde),
         spatial_field.group = rep(1, spde$n.spde))), ## a group number
  A = list(1, data_observer_A),
  tag = "observer_data")


stack_marker <- inla.stack(
  data = list(y = cbind(NA, y.pp_marker, NA),
              e = e.pp_marker), 
  effects = list(
    list( ## element 1 of effects list contains intercept, env covars and bias covar (if any)
      data.frame(
        intercept_marker = rep(1, nv + n_marker),
      sst = c(env_covariates$sst, data_marker$sst),
#      sss = c(env_covariates$sss, data_marker$sss),
#      ssh = c(env_covariates$ssh, data_marker$ssh),
      mld = c(env_covariates$mld, data_marker$mld),
      log_eke = c(env_covariates$log_eke, data_marker$log_eke),
      sst_sd = c(env_covariates$sst_sd, data_marker$sst_sd),
#      sss_sd = c(env_covariates$sss_sd, data_marker$sss_sd),
#      ssh_sd = c(env_covariates$ssh_sd, data_marker$ssh_sd),
      bathy = c(env_covariates$bathy, data_marker$bathy),
      rugosity = c(env_covariates$rugosity, data_marker$rugosity),
      bias_marker = c(biascovariate_marker_re, data_marker$bias_scaled)
      )),
    list(data.frame(spatial_field = 1:spde$n.spde),
         spatial_field.group = rep(2, spde$n.spde),
         data.frame(bias_field_marker = 1:spde$n.spde))),
  A = list(1, A.pp_marker),
  tag = "marker_data")


stack_etag <- inla.stack(
  data = list(pres = cbind(NA, NA, y.pp_etag),
              e = e.pp_etag), ## still not clear what e is but see "E" argument as input to inla() and "e.pp" in Suhaimi code
  effects = list(
    list( ## element 1 of effects list contains intercept, env covars and bias covar (if any)
      data.frame(
        intercept_etag = rep(1, nv + n_etag),
      sst = c(env_covariates$sst, data_etag$sst),
#      sss = c(env_covariates$sss, data_etag$sss),
#      ssh = c(env_covariates$ssh, data_etag$ssh),
      mld = c(env_covariates$mld, data_etag$mld),
      log_eke = c(env_covariates$log_eke, data_etag$log_eke),
      sst_sd = c(env_covariates$sst_sd, data_etag$sst_sd),
#      sss_sd = c(env_covariates$sss_sd, data_etag$sss_sd),
#      ssh_sd = c(env_covariates$ssh_sd, data_etag$ssh_sd),
      bathy = c(env_covariates$bathy, data_etag$bathy),
      rugosity = c(env_covariates$rugosity, data_etag$rugosity)
      )),
    list(data.frame(spatial_field = 1:spde$n.spde),
         spatial_field.group = rep(3, spde$n.spde), ## a group number
         data.frame(bias_field_etag = 1:spde$n.spde))),
  A = list(1, A.pp_etag),
  tag = "etag_data")



# bias = c(biascovariate, unstructured_data$bias))),

## combine the stacks
stk <- inla.stack(stack_observer, stack_marker, stack_etag)

## create a prediction stack
if (build_pred){
  stop('THIS NEEDS WORK, currently copied from Suhaimi and slightly modified')
  
  ## create an expand.grid() of desired output grid (e.g. coordinates for each cell)
  pred.grid <- expand.grid(x=seq(resolution[1]/2,
                                 max(biasfield$x),
                                 resolution[1]), 
                           y=seq(resolution[2]/2, 
                                 max(biasfield$y),
                                 resolution[2])) # make grid
  
  dim(pred.grid) 
  
  ## extract covariate values at these points
  ## ***** this is the tricky part ****
  
  ## create A.pred matrix based on mesh and pred coordinates
  A.pred <- inla.spde.make.A(mesh, loc=as.matrix(pred.grid[,1:2]))
  
  np = length(pred.grid[,1]) # number of points
  ys <- cbind(rep(NA, nrow(pred.grid)), rep(NA, nrow(pred.grid))) ## fill pred "data" with NAs
  
  stack.pred_observer =
    inla.stack(
      data = list(
        pres = ys
      ),
      effects = list(
        list( ## element 1 of effects list contains intercept, env covars and bias covar (if any)
          data.frame(
            interceptA = rep(1, np), ## assuming these intercepts mirror the corresponding data-based stack
            sst = xxx,
            bathy = xxx
          )
        ),
        list( ## element 2 of effects list contains "groupings" from Suhaimi. still unclear what they do except that field.group groups the data types. we call this one group #1
          data.frame(
            uns_field = 1:spde$n.spde ## 1:n integration points from the mesh?
          ),
          field.group = rep(1, spde$n.spde) ## a group number
        )
      ), ## close effects list
      A = list(1, 1, A.pred), ## not sure why this structure
      tag = 'pred.observer'
    )
  
  stack.pred_marker =
    inla.stack(
      data = list(
        pres = ys
      ),
      effects = list(
        list( ## element 1 of effects list contains intercept, env covars and bias covar (if any)
          data.frame(
            interceptB = rep(1, np), ## assuming these intercepts mirror the corresponding data-based stack
            sst = xxx,
            bathy = xxx
          )
        ),
        list( ## element 2 of effects list contains "groupings" from Suhaimi. still unclear what they do except that field.group groups the data types. we call this one group #1
          data.frame(
            uns_field = 1:spde$n.spde ## 1:n integration points from the mesh?
          ),
          field.group = rep(2, spde$n.spde) ## a group number
        )
      ), ## close effects list
      A = list(1, 1, A.pred), ## not sure why this structure
      tag = 'pred.marker'
    )
  
  stack.pred_etag =
    inla.stack(
      data = list(
        pres = ys
      ),
      effects = list(
        list( ## element 1 of effects list contains intercept, env covars and bias covar (if any)
          data.frame(
            interceptC = rep(1, np), ## assuming these intercepts mirror the corresponding data-based stack
            sst = xxx,
            bathy = xxx
          )
        ),
        list( ## element 2 of effects list contains "groupings" from Suhaimi. still unclear what they do except that field.group groups the data types. we call this one group #1
          data.frame(
            uns_field = 1:spde$n.spde ## 1:n integration points from the mesh?
          ),
          field.group = rep(3, spde$n.spde) ## a group number
        )
      ), ## close effects list
      A = list(1, 1, A.pred), ## not sure why this structure
      tag = 'pred.etag'
    )
  
  stk <- inla.stack(stk, stack.pred_observer, stack.pred_marker, stack.pred_etag)
  
}

#****************
## specify model formulation
## CAMELETTI VERSION
## skeptical of the UTMX + UTMY covars, aren't these projected coordinates in UTM? i believe so
#formula <- (logPM10 ~ -1 + Intercept + A + UTMX + UTMY + WS + TEMP + HMIX + PREC + EMI +
#              f(field, model=spde, group=field.group, control.group=list(model="ar1")))
## SUHAIMI VERSION
#formula = y ~ -1 + interceptA + interceptB + env + bias +
#  f(uns_field, model = spde, group = uns_field.group, control.group = list(model = 'exchangeable'))

## QUESTIONS:
## what intercepts do we need?
## how to incorporate multiple bias fields? maybe in data-specific stacks above?
## check f() and how "uns_field" object needs to be structured

formulaCorrelation = y ~ -1 + 
  intercept_observer + # observer intercept (dataset-specific)
  intercept_marker + # marker intercept (dataset-specific)
  intercept_etag + # etag intercept (dataset-specific)
  #sst +
  #f(inla.group(sst, n=10, method="quantile"), model="rw2", constr=FALSE) +
  f(inla.group(sst, n=25, method="cut"), model="rw2", constr=FALSE) +
  f(inla.group(mld, n=25, method="cut"), model="rw2", constr=FALSE) +
  f(inla.group(log_eke, n=25, method="cut"), model="rw2", constr=FALSE) +
  f(inla.group(sst_sd, n=25, method="cut"), model="rw2", constr=FALSE) +
  #f(inla.group(mld, n=10, method="quantile"), model="rw2", constr=FALSE) +
  #f(inla.group(log_eke, n=10, method="quantile"), model="rw2", constr=FALSE) +
  #f(inla.group(sst_sd, n=10, method="quantile"), model="rw2", constr=FALSE) +
#  f(inla.group(bathy, n=10, method="quantile"), model="rw2", constr=FALSE) +
#  f(inla.group(rugosity, n=10, method="quantile"), model="rw2", constr=FALSE) +
#  sss +
#  ssh +
#  mld +
#  log_eke +
#  sst_sd +
#  ssh_sd +
#  sss_sd + 
  f(inla.group(bathy, n=25, method="cut"), model="rw2", constr=FALSE) +
  f(inla.group(rugosity, n=25, method="cut"), model="rw2", constr=FALSE) +
  #bathy +
  #rugosity +
  #env + # environmental covariate (shared across datasets) estimated via joint likelihood
  bias_observer +
  bias_marker +
  f(spatial_field, model = spde, group = spatial_field.group, control.group = list(model = 'exchangeable')) + # spatial field for each dataset with spatial correlation between them
  f(bias_field_marker, model = spde) + # second spatial field (accounting for unknown bias) specific to marker dataset
  f(bias_field_etag, model = spde) # second spatial field (accounting for unknown bias) specific to etag dataset


#formula = pres ~ -1 + 
#  interceptA + interceptB + interceptC + ## intercepts for each of the three datasets
#  bias + ## how to incorporate multiple bias fields?
#  sst + bathy + ## as example env covariates to start with
#  f(uns_field, model = spde, group = field.group, control.group = list(model = 'exchangeable'))



#****************
## fit INLA
## order in our data stack will be observer (binomial), marker (poisson), etag (poisson)
result <- inla(formulaCorrelation, 
               family = c("binomial", "poisson", "poisson"),
               data = inla.stack.data(stk),
               control.predictor = list(A = inla.stack.A(stk),
                                        compute = TRUE),
               control.family = list(list(link = "cloglog"),
                                     list(link = "log"),
                                     list(link = "log")),
               E = inla.stack.data(stk)$e, ## this applies only to PO data, see ?inla argument "E"
               Ntrials = inla.stack.data(stk)$Ntrials,
               control.compute = list(#dic = TRUE, 
                 #cpo = TRUE,
                 waic = TRUE),
               verbose = TRUE)



res24 <- inla(formula24,
             family="binomial", 
             data=inla.stack.data(stk), 
             keep=FALSE,
             #control.family=list(hyper=list(prec=list(param=c(1,0.5)))),
             control.predictor=list(A=inla.stack.A(stk), compute=TRUE, link=1),#link=1
             control.inla=list(tolerance=1e-5,numint.maxfeval= 10e6),
             #control.fixed = list(expand.factor.strategy='inla'),
             control.compute = list(return.marginals=TRUE,dic=TRUE, cpo=TRUE))
  
  
  
  
  
  

# INLA STACK & FIT V1 - CAMELETTI --------------------------------------------------------------
## we need aspects of the spatiotemporal component(s) from the Cameletti code

if (run_cameletti){
  
  
  field.indices =
    inla.spde.make.index("field",
                         n.spde=spde$n.spde)#,
                         n.group=n_days)
  stack.est =
    inla.stack(data=list(logPM10=Piemonte_data$logPM10),
               A=list(A.est, 1),
               effects=
                 list(c(field.indices,
                        list(Intercept=1)),
                      list(Piemonte_data[,3:10])),
               tag="est")
  
  
  ## this DID run with our observer data as a very small, simple example
  
  field.indices <- inla.spde.make.index("field",
                                        n.spde=spde$n.spde
                                        #n.group=n_days
  )
  
  stack_marker <- inla.stack(data=list(pres = data_marker$pres),
                             A = list(data_marker_A, 1),
                             effects = list(c(field.indices,
                                              list(Intercept=1)),
                                            list(data_marker %>% select(sst))), ## grabs just sst for now
                             #list(data_observer %>% select(sst:bvfreq))), ## grabs all env covariates
                             tag="marker")
  
  stack_observer <- inla.stack(data=list(pres = data_observer$pres),
                               A = list(data_observer_A, 1),
                               effects = list(c(field.indices,
                                                list(Intercept=1)),
                                              list(data_observer %>% select(lon, lat, sst))), ## grabs just sst for now
                               #list(data_observer %>% select(sst:bvfreq))), ## grabs all env covariates
                               tag="observer")
  
  stack_etag <- inla.stack(data=list(pres = data_etag$pres),
                           A = list(data_etag_A, 1),
                           effects = list(c(field.indices,
                                            list(Intercept=1)),
                                          list(data_etag %>% select(sst))), ## grabs just sst for now
                           #list(data_observer %>% select(sst:bvfreq))), ## grabs all env covariates
                           tag="etag")
  
  #stack = inla.stack(stack_marker, stack_observer, stack_etag)
  
  ## skeptical of the UTMX + UTMY covars, aren't these projected coordinates in UTM? i believe so
  #formula <- (logPM10 ~ -1 + Intercept + A + UTMX + UTMY + WS + TEMP + HMIX + PREC + EMI + f(field, model=spde, group=field.group, control.group=list(model="ar1")))
  ## mirroring formula for Cameletti but concerned about how to include spatial component
  formula <- (pres ~ -1 + Intercept + lon + lat + sst + f(field, model=spde, group=field.group, control.group=list(model="ar1")))
  
  
  t1 <- Sys.time()  
  result =
    inla(formula,
         data=inla.stack.data(stack_observer, spde=spde),
         family="binomial",
         control.predictor=list(A=inla.stack.A(stack_observer), compute=TRUE),
         control.compute=list(cpo=FALSE),
         control.inla = list(reordering = "metis"),
         keep=FALSE, verbose=TRUE)
  t2 <- Sys.time()
  t2-t1
  
}

# INLA STACK & FIT V2 - SUHAIMI --------------------------------------------------------------
## we need aspects of the correlation modeling from the Suhaimi code
if (run_suhaimi){
  # integration stack for unstructured data ###
  max_x <- max(coarse_mesh$loc[,1])
  max_y <- max(coarse_mesh$loc[,2])
  
  loc.d <- t(matrix(c(0,0,max_x,0,max_x,max_y,0,max_y,0,0), 2))
  
  #make dual mesh
  dd <- deldir::deldir(coarse_mesh$loc[, 1], coarse_mesh$loc[, 2])
  tiles <- deldir::tile.list(dd)
  
  #make domain into spatial polygon
  domainSP <- SpatialPolygons(list(Polygons(
    list(Polygon(loc.d)), '0')))
  
  #intersection between domain and dual mesh
  poly.gpc <- as(domainSP@polygons[[1]]@Polygons[[1]]@coords, "gpc.poly")
  
  # w now contains area of voronoi polygons
  w <- sapply(tiles, function(p) rgeos::area.poly(rgeos::intersect(as(cbind(p$x, p$y), "gpc.poly"), poly.gpc)))
  
  #check some have 0 weight
  table(w>0)
  
  
  nv <- coarse_mesh$n
  n <- nrow(unstructured_data)
  
  #change data to include 0s for nodes and 1s for presences
  y.pp <- rep(0:1, c(nv, n))
  
  #add expectation vector (area for integration points/nodes and 0 for presences)
  e.pp <- c(w, rep(0, n))
  
  #diagonal matrix for integration point A matrix
  imat <- Diagonal(nv, rep(1, nv))
  
  A.pp <- rbind(imat, data_marker_A)  # new A matrix for unstructured
  
  
  #### CDB STOPPED HERE IN THIS SECTION
  
  #get covariate for integration points
  #refer plot 'covariate'
  
  covariate = dat1$gridcov[Reduce('cbind', 
                                  nearest.pixel(mesh$loc[,1], 
                                                mesh$loc[,2], 
                                                im(dat1$gridcov)))]
  
  biascovariate = biascov[Reduce('cbind',
                                 nearest.pixel(mesh$loc[,1],
                                               mesh$loc[,2],
                                               im(biascov)))]
  
  biascovariate_rasterim = biascov[Reduce('cbind',
                                 nearest.pixel(mesh$loc[,1],
                                               mesh$loc[,2],
                                               im(biascov_r)))]
  
  biascovariate_rastermatim = biascov[Reduce('cbind',
                                          nearest.pixel(mesh$loc[,1],
                                                        mesh$loc[,2],
                                                        im(as.matrix(biascov_r))))]
  
  
  # create data stacks ##
  
  #unstructured data stack with integration points
  stk_unstructured_data <- inla.stack(
    data = list(y = cbind(y.pp, NA),
                e = e.pp),
    effects = list(list(data.frame(interceptB = rep(1, nv+n),
                                   env = c(covariate, unstructured_data$env),
                                   bias = c(biascovariate, unstructured_data$bias))),
                   list(data.frame(uns_field = 1:spde$n.spde),
                        uns_field.group = rep(1, spde$n.spde))),  # named Group 1
    A = list(1, A.pp),
    tag = "unstructured_data"
  )
  
  #structured data stack
  stk_structured_data <- inla.stack(
    data = list(y = cbind(NA, structured_data$presence),
                Ntrials = rep(1, nrow(structured_data))),
    effects = list(list(data.frame(interceptA = rep(1, length(structured_data$x)),
                                   env = structured_data$env)),
                   list(data.frame(uns_field = 1:spde$n.spde),
                        uns_field.group = rep(2, spde$n.spde))), # named Group 2
    A = list(1, structured_data_A),
    tag = "structured_data"
  )
  
  # combine stacks
  stk <- inla.stack(stk_unstructured_data, stk_structured_data)
  
  # prediction stack ###
  ## CDB: this adds a prediction stack to the input stack to INLA. 
  ## I dont understand how to do this or why its needed?
  ## Cameletti has this too...main diff is Suhaimi just wrapped theirs in this custom function
  source("Create prediction stack for correlation model.R")
  join.stack <- create_prediction_stack_corr(data_stack = stk,
                                             resolution = resolution,
                                             biasfield = biasfield,
                                             dat1 = dat1,
                                             mesh = mesh,
                                             spde = spde)
  
  
  # fit model ###
  formulaC = y ~ -1 + interceptA + interceptB + env + bias +
    f(uns_field, model = spde, group = uns_field.group, control.group = list(model = 'exchangeable'))
  
  result <- inla(formulaC, 
                 family = c("poisson", "binomial"),
                 data = inla.stack.data(join.stack),
                 control.predictor = list(A = inla.stack.A(join.stack),
                                          compute = TRUE),
                 control.family = list(list(link = "log"),
                                       list(link = "cloglog")),
                 E = inla.stack.data(join.stack)$e,
                 Ntrials = inla.stack.data(join.stack)$Ntrials,
                 control.compute = list(#dic = TRUE, 
                   #cpo = TRUE,
                   waic = TRUE))
  
  #return(list(join.stack = join.stack, result = result))
  
  
}

# OLD JUNK --------------------------------------------------------------
run_junk <- FALSE
if (run_junk){
  
  ## marker tag will be as above but adding distance from shore to capture recreational fishery spatial bias? or do we need this?
  bathy <- raster::raster('~/Google Drive/Shared drives/MPG_WHOI/env_data/bathy/global_bathy_0.01.nc')
  bathy <- raster::rotate(bathy)
  bathy <- raster::crop(bathy, raster::extent(xl[1], xl[2], yl[1], yl[2]))
  
  gears <- data.table::fread('../NASA-FaCeT/scratch/cdb/CODES_gears.csv')
  data_marker <- merge(data_marker, gears %>% dplyr::select(GearCode, Code, Name), by='GearCode')
  #data.frame(data_marker %>% group_by(ObsType, GearName) %>% summarise(n=n()))
  
  u_types <- expand.grid(unique(data_marker$ObsType), unique(data_marker$Name))
  names(u_types) <- c('ObsType','Name')
  p_list <- list()
  r1 <- raster(raster::extent(xl[1], xl[2], yl[1], yl[2]), resolution=c(1))
  
  for (i in 1:nrow(u_types)){
    data_marker.i <- data_marker %>% filter(ObsType == u_types$ObsType[i] & Name == u_types$Name[i]) %>% as.data.frame()
    if (nrow(data_marker.i) == 0) next
    x <- raster::rasterize(cbind(data_marker.i$lon, data_marker.i$lat), r1, field=data_marker.i$pres, fun=function(x,...) sum(x, na.rm=T) )
    spdf <- as(x, "SpatialPixelsDataFrame")
    meso200 <- as.data.frame(spdf)
    colnames(meso200) <- c("value", "x", "y")
    
    p <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group), fill='grey60') +
      coord_fixed(xlim=xl, ylim=yl, ratio=1.3) + xlab('') + ylab('') 
    p <- p + geom_raster(data=meso200, aes(x = x, y = y, fill=value))#, breaks = list(x = sx, y = sy)) 
    #p <- p + scale_x_continuous(breaks= x.at, labels=x.labels)
    #p <- p + scale_y_continuous(breaks= y.at, labels=y.labels)
    #p <- p + scale_fill_gradientn("presence count")#, colours = jet.colors(10), 
    #                              limits=c(0,50), oob=squish,
    #                              breaks = seq(0, 50, by=10),
    #                              labels = c(seq(0, 40, by=10), '>50'),
    #                              guide = guide_colorbar(barwidth = .5, barheight = 15)) 
    p <- p + theme_bw(base_size = 14) + theme(panel.grid=element_blank()) + ggtitle(paste0(u_types$ObsType[i], ' - ', u_types$Name[i], '(n=', nrow(data_marker.i),')'))#ggtitle('All species - standard deviation')#+ guides(fill=FALSE) #+ ggtitle(sp[i]) + 
    #p <- p + facet_wrap_paginate(. ~ platform, ncol=4, nrow=3) #+ theme_bw(base_size = 10) + theme(panel.grid=element_blank())
    p_list[[i]] <- p
    
  }
  lay <- rbind(c(1,2,3,4),
               c(5,6,7,8),
               c(9,10,11,12),
               c(13,14,15,16))
  g <- gridExtra::arrangeGrob(grobs = p_list, heights = c(4,4,4,4),
                              width = c(6,6,6,6), layout_matrix = lay)
  ggsave(file = paste('./compare_marker_gearcodes.png', sep=''), width=26, height=18, units = 'in', g)
  
  ## convert "sport" to RR
  data_marker$GearCode[which(data_marker$Name == 'sport')] <- 'RR'
  data_marker$Name[which(data_marker$Name == 'sport')] <- 'rood & reel'
  
  
  
  ggplot(data_marker, aes(x=lon, y=lat)) +
    rasterise(geom_point()) +
    facet_wrap(ObsType ~ Name)
  
  g_tmax_map <- ggplot(data = tmax_Jan_09_df) +
    geom_raster(aes(x = x, y = y, fill = `2009-01-01`)) +
    scale_fill_viridis_c() +
    theme_void() +
    theme(
      legend.position = "bottom"
    )
  
  
  sum_n <- data_marker %>% filter(ObsType == 'release') %>% group_by(ObsType, Code) %>% summarise(n=n())
  sum_n$perc <- round(sum_n$n / sum(sum_n$n), 2)
  sum_n
  
  data_marker_rr <- data_marker %>% filter()
  
  ## distance to bathy contour for marker data
  bcontour <- rasterToContour(bathy, levels=round(median(data_marker$bathy, na.rm=T), 0))
  paste0('marker - median bathy=', abs(round(median(data_marker$bathy, na.rm=T),0)), 'm')
  proj4string(bcontour) <- proj4string(bathy)
  r1 <- raster(raster::extent(xl[1], xl[2], yl[1], yl[2]), resolution=c(.05))
  dd <- gDistance(bcontour, as(r1,"SpatialPoints"), byid=T)
  r1[] = apply(dd, 1, min)
  r1_scaled <- (r1 / cellStats(r1, max))*-1 + 1
  r1_scaled <- raster::resample(r1_scaled, bathy)
  r1_scaled <- raster::mask(r1_scaled, bathy)
  r1_binary <- raster::reclassify(r1_scaled, c(c(0, 0.9, 0), c(0.9, 1, 1)))
  r1_scaled_0.9 <- raster::reclassify(r1_scaled, c(0, 0.9, NA))
  r.min = cellStats(r1_scaled_0.9, "min")
  r.max = cellStats(r1_scaled_0.9, "max")
  r1_scaled_0.9 <- ((r1_scaled_0.9 - r.min) / (r.max - r.min))
  r1_scaled_0.9_marker <- raster::mask(raster::reclassify(r1_scaled_0.9, c(NA, NA, 0)), r1_binary)
  #plot(r1_scaled_0.9_marker)
  
  ## distance to bathy contour for observer data
  bcontour <- rasterToContour(bathy, levels=round(median(data_observer$bathy, na.rm=T), 0))
  paste0('observer - median bathy=', abs(round(median(data_observer$bathy, na.rm=T),0)), 'm')
  proj4string(bcontour) <- proj4string(bathy)
  r1 <- raster(raster::extent(xl[1], xl[2], yl[1], yl[2]), resolution=c(.05))
  dd <- gDistance(bcontour, as(r1,"SpatialPoints"), byid=T)
  r1[] = apply(dd, 1, min)
  r1_scaled <- (r1 / cellStats(r1, max))*-1 + 1
  r1_scaled <- raster::resample(r1_scaled, bathy)
  r1_scaled <- raster::mask(r1_scaled, bathy)
  r1_binary <- raster::reclassify(r1_scaled, c(c(0, 0.9, 0), c(0.9, 1, 1)))
  r1_scaled_0.9 <- raster::reclassify(r1_scaled, c(0, 0.9, NA))
  r.min = cellStats(r1_scaled_0.9, "min")
  r.max = cellStats(r1_scaled_0.9, "max")
  r1_scaled_0.9 <- ((r1_scaled_0.9 - r.min) / (r.max - r.min))
  r1_scaled_0.9_observer <- raster::mask(raster::reclassify(r1_scaled_0.9, c(NA, NA, 0)), r1_binary)
  #plot(r1_scaled_0.9_observer)
  
  ## convex hull polygons
  ## trick the hull so doesnt cut off gulf of maine
  data_observer_mod <- rbind(data_observer %>% select(lon, lat), rbind(c(-70, 45),
                                                                       c(-74.68333, 11.71667)), use.names=F) 
  hull_observer <- terra::convHull(terra::vect(data_observer_mod))
  hull_observer <- as(hull_observer, "Spatial")
  ## trick the hull so doesnt cut off gulf of mexico
  data_marker_mod <- rbind(data_marker %>% select(lon, lat), rbind(#c(-100, 25),
    c(-96.36667, 26.13333)), use.names=F) ## trick the hull so doesnt cut off gulf of maine
  hull_marker <- terra::convHull(terra::vect(data_marker_mod %>% select(lon, lat)))
  hull_marker <- as(hull_marker, "Spatial")
  
  ## mask by resulting hull
  r1_scaled_0.9_observer_mask <- raster::mask(r1_scaled_0.9_observer, hull_observer)
  r1_scaled_0.9_marker_mask <- raster::mask(r1_scaled_0.9_marker, hull_marker)
  
  ## combine a rod reel bathy-based field with a longline bathy-based field and weight them according to %s in iccat dataset
  marker_bias <- sum(r1_scaled_0.9_observer_mask * .3, ## roughly 30% of marker tag data comes from longline effort
                     r1_scaled_0.9_marker_mask * .7, na.rm=T) ## roughly 70% of marker tag data comes from rod and reel effort
  marker_bias <- raster::mask(marker_bias, bathy)
  
  ## take a peek
  pdf('bathy_bias_field.pdf', height=12, width=8)
  par(mfcol=c(3,1))
  plot(r1_scaled_0.9_marker, main=paste0('marker - median bathy=', abs(round(median(data_marker$bathy, na.rm=T),0)), 'm'))
  plot(hull_marker, add=T)
  points(data_marker$lon, data_marker$lat)
  plot(r1_scaled_0.9_observer, main=paste0('observer - median bathy=', abs(round(median(data_observer$bathy, na.rm=T),0)), 'm'))
  plot(hull_observer, add=T)
  points(data_observer$lon, data_observer$lat)
  
  plot(marker_bias, main='bathy bias field = 70% marker & 30% observer')
  dev.off()
  
  
  ## example env covariate grid that those variables are extracted from and that we will ultimately predict to
  ## uses native GLORYS grid and resolution
  #mld <- raster::rotate(raster::stack('/Volumes/Elements/roms_nwa/mld.GLORYS.NWAtl.mon.mean.1993-2018.nc'))
  #mld <- mld[[1]]
  #env_coords <- xyFromCell(mld, 1:ncell(mld))
  
}


r <- raster(xm=0, ymn=0, xmx=100, ymx=100, res=1)
