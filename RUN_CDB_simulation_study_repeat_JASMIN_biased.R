# simulation - zero bias (env covariate)

# simulation study

# load libraries ####

library(INLA)
#INLA:::inla.dynload.workaround() 
library(reshape2)
library(deldir)
library(rgeos)
library(fields)
library(RColorBrewer)
#display.brewer.all(colorblindFriendly = TRUE)

### Set parameters ####

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
#probs = rep(c(0.5, 0.3, 0.1, 0.05, 0.01),5) #default
#probs = rep(c(0.2, 0.2, 0.2, 0.2, 0.2),5) # unbiased
probs = rep(c(0.5, 0.4, 0.1, 0.01, 0.001),5) # very biased

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

# Generate data #####

# change the seed for each run

i <- as.numeric(commandArgs(trailingOnly = TRUE))[1]

#seed_new <- seed_all[i,2]
seed_new <- 399

print(seed_new)

#start_time_ori <- Sys.time()

  start_time <- Sys.time()
  #seed = seed_all[i]

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
  #  image.plot(list(x=dat1$Lam$xcol, y=dat1$Lam$yrow, z=t(dat1$rf.s)),
  #             main='Thinned unstructured data', asp=1,
  #             col = hcl.colors(12, "RdBu", rev = TRUE))
  #  points(unstructured_data$x, unstructured_data$y, pch = 20)#note rescale again - plotting back on original

  #' Visualise structured data
  #+ echo = FALSE
  #par(mfrow=c(1,1))
   # image.plot(list(x=dat1$Lam$xcol, y=dat1$Lam$yrow, z=t(dat1$rf.s)),
   #            main='Structured data', asp=1,
   #            col = hcl.colors(12, "RdBu", rev = TRUE))
   # points(structured_data$x,structured_data$y, pch = 21, bg = structured_data$presence, col = "black")
    #par(xpd = TRUE)
    #legend(0,115,c("Absence", "Presence"), pch = 21, col = "black", pt.bg = c(0,1))

  ##########################################################
  #' ### Correlation model
  #' set prior for fixed effect
  #'
  #source("Run correlation-bias model - corrected.R")
  #mod_12 <- correlationbias_model(unstructured_data = unstructured_data,
  #                                structured_data = structured_data,
  #                                dat1,
  #                                biasfield,
  #                                biascov,
  #                                dim,
  #                                plotting=F,
  #                                mesh.edge = mesh.edge,
  #                                mesh.offset = mesh.offset,
  #                                resolution = resolution)

  
  # Function to run joint model with spatial correlation with bias covariate.
  #
  # Function will return a list of data stack and the fitted model.
  
  #correlationbias_model <- function(unstructured_data,
  #                              structured_data,
  #                              dat1,
  #                              biasfield,
  #                              biascov,
  #                              dim,
  #                              plotting = FALSE,
  #                              mesh.edge = c(20,40),
  #                              mesh.offset = c(5,20),
  #                              resolution = c(10,10)){
  
  # packages
  library(INLA)
  library(reshape2)
  library(rgeos)
  library(fields)
  library(deldir)
  
  ###############################
  
  # preparation - mesh contruction - use the loc.domain argument ####
  
  mesh <- inla.mesh.2d(loc.domain = biasfield[, c(1,2)],
                       max.edge = mesh.edge,
                       cutoff = 2,
                       offset = mesh.offset)
  
  # make SPDE ####
  spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
  
  
  # make A matrix
  structured_data_A <- inla.spde.make.A(mesh = mesh,
                                        loc = as.matrix(structured_data[, 2:3]))
  
  unstructured_data_A <- inla.spde.make.A(mesh = mesh,
                                          loc = as.matrix(unstructured_data[ , 1:2]))
  
  
  # integration stack for unstructured data ####
  max_x <- max(biasfield$x)
  max_y <- max(biasfield$y)
  
  loc.d <- t(matrix(c(0,0,max_x,0,max_x,max_y,0,max_y,0,0), 2))
  
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
  
  #check some have 0 weight
  table(w>0)
  
  
  nv <- mesh$n
  n <- nrow(unstructured_data)
  
  #change data to include 0s for nodes and 1s for presences
  y.pp <- rep(0:1, c(nv, n))
  
  #add expectation vector (area for integration points/nodes and 0 for presences)
  e.pp <- c(w, rep(0, n))
  
  #diagonal matrix for integration point A matrix
  imat <- Diagonal(nv, rep(1, nv))
  
  A.pp <- rbind(imat, unstructured_data_A)  # new A matrix for unstructured
  
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
  
  
  # create data stacks ####
  
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
    
    tag = "unstructured_data")
  
  #structured data stack
  stk_structured_data <- inla.stack(
    
    data = list(y = cbind(NA, 
                          structured_data$presence),
                
                Ntrials = rep(1, 
                              nrow(structured_data))),
    
    effects = list(list(data.frame(
      interceptA = rep(1, length(structured_data$x)),
      
      env = structured_data$env)),
      
      list(data.frame(uns_field = 1:spde$n.spde),
           uns_field.group = rep(2, spde$n.spde))), # named Group 2
    
    A = list(1, structured_data_A),
    
    tag = "structured_data"
  )
  
  # combine stacks
  stk <- inla.stack(stk_unstructured_data, stk_structured_data)
  
  # prediction stack ####
  source("Create prediction stack for correlation model.R")
  
  join.stack <- create_prediction_stack_corr(data_stack = stk,
                                             resolution = resolution,
                                             biasfield = biasfield,
                                             dat1 = dat1,
                                             mesh = mesh,
                                             spde = spde)
  
  
  ##############################
  
  # fit model ####
  
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
  
  # no plotting for grouped model
  
  #  return(list(join.stack = join.stack, result = result))
  
  #}
  
  mod_12 <- list(join.stack = join.stack, result = result)
  
  source("validation_function for correlation.R")
  validation_12 <- validation_function_str(result=mod_12[[2]],
                                       resolution=resolution,
                                       join.stack=mod_12[[1]],
                                       model_type="correlationbias_str",
                                       unstructured_data = unstructured_data,
                                       structured_data = structured_data,
                                       dat1 = dat1,
                                       summary_results=T,
                                       absolute=TRUE,
                                       dim = dim,
                                       plotting = T)

  validation_13 <- validation_function_uns(result=mod_12[[2]],
                                       resolution=resolution,
                                       join.stack=mod_12[[1]],
                                       model_type="correlationbias_uns",
                                       unstructured_data = unstructured_data,
                                       structured_data = structured_data,
                                       dat1 = dat1,
                                       summary_results=T,
                                       absolute=TRUE,
                                       dim = dim,
                                       plotting = F)

save.image('run_suhaimi_model12.rda')

