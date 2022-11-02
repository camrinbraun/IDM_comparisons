mesh <- inla.mesh.2d()

#make A matrix for structured data
observer_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix())

#make A matrix for unstructured data
marker_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix())
etag_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix())

# create integration stack
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

nv <- mesh$n
n_marker <- nrow(marker_data)
n_etag <- nrow(etag_data)

#change data to include 0s for nodes and 1s for presences
y.pp_marker <- rep(0:1, c(nv, n_marker))
y.pp_etag <- rep(0:1, c(nv, n_etag))

#add expectation vector (area for integration points/nodes and 0 for presences)
e.pp_marker <- c(w, rep(0, n_marker))
e.pp_etag <- c(w, rep(0, n_etag))

#diagonal matrix for integration point A matrix
imat <- Diagonal(nv, rep(1, nv))

A.pp_marker <- rBind(imat, marker_data_A)
A.pp_etag <- rBind(imat, etag_data_A)

#get (bias) covariate(s) for integration points
covariate = dat1$gridcov[Reduce('cbind', nearest.pixel(mesh$loc[,1], mesh$loc[,2], im(dat1$gridcov)))]
biascovariate = biascov[Reduce('cbind', nearest.pixel(mesh$loc[,1], mesh$loc[,2], im(biascov)))]

#unstructured data stack with integration points
stk_marker_data <- inla.stack(data=list(y=cbind(y.pp_marker, NA), e = e.pp_marker),
                              effects=list(list(data.frame(intercept_marker=rep(1,nv+n_marker)), 
                                                env = c(covariate, marker_data$env)), 
                                           list(data.frame(spatial_field=1:spde$n.spde), 
                                                spatial_field.group=rep("marker",spde$n.spde), 
                                                data.frame(bias_field_marker = 1:spde$n.spde))),
                              A=list(1,A.pp_marker),
                              tag="marker_data")	

stk_etag_data <- inla.stack(data=list(y=cbind(y.pp_etag, NA), e = e.pp_etag),
                            effects=list(list(data.frame(intercept_etag=rep(1,nv+n_etag)), 
                                              env = c(covariate, etag_data$env)), 
                                         list(data.frame(spatial_field=1:spde$n.spde), 
                                              spatial_field.group=rep("etag",spde$n.spde), 
                                              data.frame(bias_field_etag = 1:spde$n.spde))),
                            A=list(1,A.pp_etag),
                            tag="etag_data")	

#stack for structured data
stk_observer_data <- inla.stack(data=list(y=cbind(NA, observer_data$presence), Ntrials = rep(1, nrow(observer_data))),
                                effects=list(list(data.frame(intercept_observer=rep(1,length(observer_data$x)), 
                                                             env = observer_data$env, 
                                                             bias_observer_totalhooks = observer_data$totalhooks, 
                                                             bias_observer_soakdur = observer_data$soakdur)), 
                                             list(data.frame(spatial_field=1:spde$n.spde), 
                                                  spatial_field.group=rep("observer",spde$n.spde))),
                                A=list(1,observer_data_A),
                                tag="observer_data")

formulaCorrelation = y ~ -1 + 
  intercept_observer + # observer intercept (dataset-specific)
  intercept_marker + # marker intercept (dataset-specific)
  intercept_etag + # etag intercept (dataset-specific)
  env + # environmental covariate (shared across datasets) estimated via joint likelihood
  bias_observer_totalhooks + # bias covariate specific to observer dataset
  bias_observer_soakdur + # bias covariate specific to observer dataset
  f(spatial_field, model = spde, group = spatial_field.group, control.group = list(model = 'exchangeable')) + # spatial field for each dataset with spatial correlation between them
  f(bias_field_marker, model = spde) + # second spatial field (accounting for unknown bias) specific to marker dataset
  f(bias_field_etag, model = spde) # second spatial field (accounting for unknown bias) specific to etag dataset

