## from https://github.com/NERC-CEH/IDM_comparisons and Suhaimi et al 2021 


t1 <- Sys.time()
sink(paste0('file_', format(t1, '%Y-%m-%dT%H%M%SZ'), '.txt'))


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

sp_subset = FALSE

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
  data_observer$hookhours <- data_observer$NUMBER_HOOKS_SET * data_observer$SOAK_TIME
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
        sss = data_observer$sss,
        ssh = data_observer$ssh,
        mld = data_observer$mld,
        log_eke = data_observer$log_eke,
        sst_sd = data_observer$sst_sd,
        sss_sd = data_observer$sss_sd,
        ssh_sd = data_observer$ssh_sd,
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
        sss = c(env_covariates$sss, data_marker$sss),
        ssh = c(env_covariates$ssh, data_marker$ssh),
        mld = c(env_covariates$mld, data_marker$mld),
        log_eke = c(env_covariates$log_eke, data_marker$log_eke),
        sst_sd = c(env_covariates$sst_sd, data_marker$sst_sd),
        sss_sd = c(env_covariates$sss_sd, data_marker$sss_sd),
        ssh_sd = c(env_covariates$ssh_sd, data_marker$ssh_sd),
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
        sss = c(env_covariates$sss, data_etag$sss),
        ssh = c(env_covariates$ssh, data_etag$ssh),
        mld = c(env_covariates$mld, data_etag$mld),
        log_eke = c(env_covariates$log_eke, data_etag$log_eke),
        sst_sd = c(env_covariates$sst_sd, data_etag$sst_sd),
        sss_sd = c(env_covariates$sss_sd, data_etag$sss_sd),
        ssh_sd = c(env_covariates$ssh_sd, data_etag$ssh_sd),
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
  #bathy +
  f(inla.group(sst, n=25, method="cut"), model="rw2", constr=FALSE) +
  f(inla.group(mld, n=25, method="cut"), model="rw2", constr=FALSE) +
  #f(inla.group(log_eke, n=25, method="cut"), model="rw2", constr=FALSE) +
  f(inla.group(sst_sd, n=25, method="cut"), model="rw2", constr=FALSE) +
  f(inla.group(bathy, n=25, method="cut"), model="rw2", constr=FALSE) +
  f(inla.group(rugosity, n=25, method="cut"), model="rw2", constr=FALSE) +
  ###f(inla.group(sst, n=10, method="quantile"), model="rw2", constr=FALSE) +
  ###f(inla.group(mld, n=10, method="quantile"), model="rw2", constr=FALSE) +
  ###f(inla.group(log_eke, n=10, method="quantile"), model="rw2", constr=FALSE) +
  ###f(inla.group(sst_sd, n=10, method="quantile"), model="rw2", constr=FALSE) +
  ###f(inla.group(bathy, n=10, method="quantile"), model="rw2", constr=FALSE) +
  ###f(inla.group(rugosity, n=10, method="quantile"), model="rw2", constr=FALSE) +
  bias_observer +
  bias_marker +
  f(spatial_field, model = spde, group = spatial_field.group, control.group = list(model = 'exchangeable')) + # spatial field for each dataset with spatial correlation between them
  f(bias_field_marker, model = spde) + # second spatial field (accounting for unknown bias) specific to marker dataset
  f(bias_field_etag, model = spde) # second spatial field (accounting for unknown bias) specific to etag dataset



#****************
## fit INLA
## order in our data stack will be observer (binomial), marker (poisson), etag (poisson)
inla.getOption()
#inla.setOption(inla.mode = 'experimental')
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
               verbose = TRUE, 
               safe = TRUE)

t2 <- Sys.time()

print(str(data_marker))
print(str(data_etag))
print(str(data_observer))
print(formulaCorrelation)
print(paste0('Start time ', t1))
print(paste0('End time ', t2))

save(result, file=paste0('file_', format(t1, '%Y-%m-%dT%H%M%SZ'), '.rda'))
