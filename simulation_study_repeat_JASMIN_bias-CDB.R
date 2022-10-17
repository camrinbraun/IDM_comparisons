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

## load original data from FaCeT project? default is FALSE
if (!exists("get_data")) { get_data = FALSE }

## use small spatial subset of the data? default is TRUE
if (!exists("sp_subset")) { sp_subset = TRUE }

## want to simulate some data? this is from Suhaimi and is not tested but retained for posterity
if (!exists("sim_data")) { sim_data = FALSE }


# PRE-REQS & NOTES --------------------------------------------------------------

# observer binomial
# everything else poisson

#bias covariates:
#  observer is canyon-related/focused
#  mark-recapture is distance to shore + bathy something
#  etag is distance from previous day

# pairwise correlations between each



## something to do with the strata?
#biasfield

# x by y matrix of bias weightings (0 to 1, where 1 is most biased)
#biascov

# mesh parameters
#mesh.edge
#mesh.offset

#dat1
# dat1$gridcov = environmental covariate grid

# inla "stk" is how data is prepped for inla. simply vectors of covariates (e.g. env, bias, etc) and the accompanying grid
# plus some mesh-derived stuff. dont need any info on pres/abs or even location bc its PO
#stk_unstructured_data

# structured data reqs info on pres/abs in addition to the above
#stk_structured_data


# GENERATE DATA FOR SIMS --------------------------------------------------------------
## we have real data so dont need this but leaving here as a guide for data formats, etc
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
  
  seed_all <- read.csv("seed_biased.csv")
  #seed_new <- 399
  
  # Generate data ##
  
  # change the seed for each run
  
  i <- as.numeric(commandArgs(trailingOnly = TRUE))[1]
  
  seed_new <- seed_all[i,2]
  
  print(seed_new)
  
  #start_time_ori <- Sys.time()
  
  start_time <- Sys.time()
  seed = seed_all[i]
  
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

## cols needed are x (lon), y (lat), presence (0 or 1), ...

if (get_data){
  ## observer data "enhanced" with env covariates, incl all data that passed QC
  data_obs <- data.table::fread('../NASA-FaCeT/data/enhance/observer-monthly/observer-enhanced_allAbs.csv')
  data_obs <- data_obs %>% filter(SpeciesCode == 'bsh') # results in approx 22k sets, 7300 have presence
  data_obs %>% group_by(pres) %>% summarise(n=n())
  data.table::fwrite(data_obs, '~/Google Drive/Shared drives/MPG_WHOI/data/bsh_inla/bsh_observer.csv')
  
  ## mark-recapture data "enhanced" with env covariates, incl all data that passed QC
  data_marker <- data.table::fread('../NASA-FaCeT/data/enhance/iccat-monthly/iccat-enhanced.csv')
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
  data_marker <- data.table::fread('~/Google Drive/Shared drives/MPG_WHOI/data/bsh_inla/bsh_marker.csv')
  data_etag <- data.table::fread('~/Google Drive/Shared drives/MPG_WHOI/data/bsh_inla/bsh_etag.csv')
  
}

if (sp_subset){
  
  print(paste("**---- Using small subset of study domain to speed things up ----**"))
  
  xl <- c(-70, -65)
  yl <- c(35, 40)
  
  data_observer <- data_observer %>% filter(lon > xl[1] & lon < xl[2] & lat > yl[1] & lat < yl[2])
  data_marker <- data_marker %>% filter(lon > xl[1] & lon < xl[2] & lat > yl[1] & lat < yl[2])
  data_etag <- data_etag %>% filter(lon > xl[1] & lon < xl[2] & lat > yl[1] & lat < yl[2])
  
} else{
  
  print(paste("**---- Using full N Atlantic study domain ----**"))
  
  xl <- c(-100, -5)
  yl <- c(10, 55)
  
}

## these datasets each need a "bias" col and whatever env covars will be used

#"dat1" needs dat1$gridcov which is matrix of env grid (which i think is just used for getting dimensions?)


# CHECK ENV COVARIATES (COLLINEARITY, ETC) --------------------------------------------------------------
## NA values should have already been removed from all env covariates


## do we need to standardize env? Cameletti did like:
#mean_covariates = apply(data_observer %>% select(lon, lat, sst:bvfreq), 2, mean)
#sd_covariates = apply(data_observer %>% select(lon, lat, sst:bvfreq), 2, sd)

#data_observer_standardized <- scale(data_observer %>% select(lon, lat, sst:bvfreq),
#                                    mean_covariates, sd_covariates)

# BUILD BIAS FIELDS --------------------------------------------------------------

## observer will be distance from bathy contour (i.e. representing shelf break as focal point of effort)


## marker tag will be as above but adding distance from shore to capture recreational fishery spatial bias


## etag will be daily distance from previous location


# BUILD MESH --------------------------------------------------------------
## some guidance here: https://punama.github.io/BDI_INLA/
## and a book on spatiotemporal models with INLA / R: https://sites.google.com/a/r-inla.org/stbook/


## example env covariate grid that those variables are extracted from and that we will ultimately predict to
## uses native GLORYS grid and resolution
#mld <- raster::rotate(raster::stack('/Volumes/Elements/roms_nwa/mld.GLORYS.NWAtl.mon.mean.1993-2018.nc'))
#mld <- mld[[1]]
#env_coords <- xyFromCell(mld, 1:ncell(mld))

# construct the mesh ##
proj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
atl <- rgdal::readOGR('../NASA-FaCeT/data/shapefiles/atlantic.shp')
atl <- raster::crop(atl, raster::extent(xl[1], xl[2], yl[1], yl[2]))

bdry <- inla.sp2segment(atl)
bdry$loc <- inla.mesh.map(bdry$loc)
coarse_mesh <- inla.mesh.2d(boundary = bdry, max.edge=c(2,4), 
                      offset = c(0.5, 1),
                      cutoff = 0.3)
coarse_mesh$crs <- proj

#mesh_plot <- ggplot() +
#  gg(coarse_mesh) +
#  ggtitle('Plot of mesh') +
#  theme_bw() +
#  theme(plot.title = element_text(hjust = 0.5))
#mesh_plot


# MAKE SPDE & MATRICES --------------------------------------------------------------

# make SPDE 
spde <- inla.spde2.matern(mesh = coarse_mesh, alpha = 2)

# make A matrix
data_marker_A <- inla.spde.make.A(mesh = coarse_mesh,
                                  loc = as.matrix(data_marker[,c('lon','lat')])
                                  #group=Piemonte_data$time, ## these temporal groupings are from Cameletti, not sure if we need them?
                                  #n.group=n_days
                                  )

data_observer_A <- inla.spde.make.A(mesh = coarse_mesh,
                                  loc = as.matrix(data_observer[,c('lon','lat')])
                                  #group=Piemonte_data$time, ## these temporal groupings are from Cameletti, not sure if we need them?
                                  #n.group=n_days
                                  )
## ^^^ do we need a temporal grouping as Cameletti did? I dont see why we would.

data_etag_A <- inla.spde.make.A(mesh = coarse_mesh,
                                loc = as.matrix(data_etag[,c('lon','lat')])
                                #group=Piemonte_data$time, ## these temporal groupings are from Cameletti, not sure if we need them?
                                #n.group=n_days
                                )
## ^^^ perhaps grouping above is way to account for individual in etag data?


# INLA STACK & FIT V1 - CAMELETTI --------------------------------------------------------------

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


# INLA STACK & FIT V2 - SUHAIMI --------------------------------------------------------------

# integration stack for unstructured data ###
#max_x <- max(env_coords[,1])
max_x <- max(coarse_mesh$loc[,1])
#max_y <- max(env_coords[,2])
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
  
