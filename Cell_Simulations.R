
library(doParallel)
library(parallel)
library(foreach)
library(spatstat)
library(sp)
library(rstan)
library(brms)
library(fields)
library(raster)
library(gridExtra)
library(ggplot2)

## Procedure
## 1000 simulations
## 1) Simulate full data
## 2) Fit full model
## 3) Simulate presence-only data set
## 4) Fit presence-only model
## 5) Record full data number of ven and FSW,
##      full model predicted ven and FSW,
##      presence-only model predicted ven and FSW

set.seed(500)
## Create the domain, the unit square
cells = raster(nrows=45, ncols=45, xmn=0, xmx=1, ymn = 0, ymx = 1)
crs(cells) = NA
cells = as(cells, 'SpatialPixels')
plot(cells)
n.cell = ncell(cells)

cells.small = raster(nrow = 71, ncol = 71, xmn = 0, xmx = 1, ymn = 0, ymx = 1)
crs(cells.small) = NA
cells.large = raster(nrow = 22, ncol = 22, xmn = 0, xmx = 1, ymn = 0, ymx = 1)
crs(cells.large) = NA
cells.small = as(cells.small, 'SpatialPixels')
cells.large = as(cells.large, 'SpatialPixels')
n.cell.small = length(cells.small)
n.cell.large = length(cells.large)

## Simulate the covariates at in the square according to MVN
## Simulate at very fine grid in the square
## One is very smooth, one is medium smooth, and one is rough
cov.grid = list(x = seq(0.001, 0.999, length.out = 100), y = seq(0.001, 0.999, length.out = 100))
cov.grid.expand = expand.grid(cov.grid$x, cov.grid$y)

obj.smooth = Exp.image.cov(grid = cov.grid, theta = 0.23, setup = TRUE)
cov.vs.smooth <- sim.rf(obj.smooth)
cov.vc.smooth <- sim.rf(obj.smooth)
obj.medium = Exp.image.cov(grid = cov.grid, theta = 0.07, setup = TRUE)
cov.vs.medium <- sim.rf(obj.medium)
cov.vc.medium <- sim.rf(obj.medium)
obj.rough = Exp.image.cov(grid = cov.grid, theta = 0.01, setup = TRUE)
cov.vs.rough <- sim.rf(obj.rough)
cov.vc.rough <- sim.rf(obj.rough)

par(mfrow = c(2,3))
image.plot(cov.vs.smooth)
image.plot(cov.vs.medium)
image.plot(cov.vs.rough)
image.plot(cov.vc.smooth)
image.plot(cov.vc.medium)
image.plot(cov.vc.rough)

cov.vs = data.frame(x = cov.grid.expand[,1], y = cov.grid.expand[,2],
                    smooth = c(cov.vs.smooth),
                    medium1 = c(cov.vs.medium),
                    rough = c(cov.vs.rough),
                    medium2 = c(cov.vc.medium))
cov.vc = data.frame(x = cov.grid.expand[,1], y = cov.grid.expand[,2],
                    smooth = c(cov.vc.smooth),
                    medium1 = c(cov.vs.medium),
                    rough = c(cov.vc.rough),
                    medium2 = c(cov.vc.medium))
cov.vs.sp = cov.vs
coordinates(cov.vs.sp) = ~ x + y
cov.vc.sp = cov.vc
coordinates(cov.vc.sp) = ~ x + y

ggplot() +
  geom_tile(mapping = aes(x = x, y = y, fill = smooth, col = smooth), data = cov.vs) +
  coord_equal() + scale_fill_gradientn(colours = tim.colors()) +
  scale_color_gradientn(colours = tim.colors())




# Now get covariates at each cell -----------------------------------------

## Get covariates for same cells
cells.vs = matrix(NA, nrow = n.cell, ncol = ncol(cov.vs) - 2)
cells.vc = matrix(NA, nrow = n.cell, ncol = ncol(cov.vc) - 2)
for(i in 1:n.cell){
  tmp = over(cov.vs.sp, cells[i,])
  cells.vs[i,] = colMeans(cov.vs[which(tmp == 1),-c(1,2)])
  cells.vc[i,] = colMeans(cov.vc[which(tmp == 1),-c(1,2)])
}
cells.vs = as.data.frame(cells.vs)
names(cells.vs) = names(cov.vs)[-c(1,2)]
cells.vc = as.data.frame(cells.vc)
names(cells.vc) = names(cov.vc)[-c(1,2)]


## Get covariates for small cells
cells.vs.small = matrix(NA, nrow = n.cell.small, ncol = ncol(cov.vs) - 2)
cells.vc.small = matrix(NA, nrow = n.cell.small, ncol = ncol(cov.vc) - 2)
for(i in 1:n.cell.small){
  tmp = over(cov.vs.sp, cells.small[i,])
  cells.vs.small[i,] = colMeans(cov.vs[which(tmp == 1),-c(1,2)])
  cells.vc.small[i,] = colMeans(cov.vc[which(tmp == 1),-c(1,2)])
}
cells.vs.small = as.data.frame(cells.vs.small)
names(cells.vs.small) = names(cov.vs)[-c(1,2)]
cells.vc.small = as.data.frame(cells.vc.small)
names(cells.vc.small) = names(cov.vc)[-c(1,2)]

## Get covariates for large cells
cells.vs.large = matrix(NA, nrow = n.cell.large, ncol = ncol(cov.vs) - 2)
cells.vc.large = matrix(NA, nrow = n.cell.large, ncol = ncol(cov.vc) - 2)
for(i in 1:n.cell.large){
  tmp = over(cov.vs.sp, cells.large[i,])
  cells.vs.large[i,] = colMeans(cov.vs[which(tmp == 1),-c(1,2)])
  cells.vc.large[i,] = colMeans(cov.vc[which(tmp == 1),-c(1,2)])
}
cells.vs.large = as.data.frame(cells.vs.large)
names(cells.vs.large) = names(cov.vs)[-c(1,2)]
cells.vc.large = as.data.frame(cells.vc.large)
names(cells.vc.large) = names(cov.vc)[-c(1,2)]



## Set betas
betas.vs.bin = rnorm(4, sd = 0.3)
betas.vs.nb = rnorm(4, sd = 0.6)
betas.vc.bin = rnorm(4)
betas.vc.nb = rnorm(4, sd = 0.4)

uniform = FALSE
include.cov = FALSE







# Do the first run --------------------------------------------------------

cells.vs.y = rep(NA, n.cell)
cells.vc.y = rep(NA, n.cell)

# Simulate the actual process ---------------------------------------------

##
#### Simulate venue locations
##

## Simulate data
cells.vc.y = rbinom(n.cell, 1, prob = plogis(-1 + betas.vc.bin %*% t(cells.vc))) *
  rnbinom(n.cell, mu = exp(0 + betas.vc.nb %*% t(cells.vc)), size = 10)

## Now simulate locations inside each cell
## Now also grab covariates associated with it
venues.sim = matrix(NA, nrow = 0, ncol = 2)
for(i in 1:n.cell){
  if(cells.vc.y[i] > 0){
    cells.xy = spsample(cells[i], cells.vc.y[i], type = "random")
    while(length(cells.xy) != cells.vc.y[i]){
      cells.xy = spsample(cells[i], cells.vc.y[i], type = "random")
    }
    venues.sim = rbind(venues.sim, cbind(coordinates(cells.xy), cells.vc[i,], cells.vs[i,]))
  }
}
venues.sim = as.data.frame(venues.sim)
venues.sim = venues.sim[,-c(8,10)]
names(venues.sim) = c("x","y", "smooth.vc", "medium1", "rough.vc", "medium2", "smooth.vs", "rough.vs")

## Simulate venue size
venues.y = rbinom(nrow(venues.sim), 1,
                  prob = plogis(1 + betas.vs.bin %*% t(venues.sim[,c(7, 4, 8, 6)]))) *
  rlnorm(nrow(venues.sim), 2 + betas.vs.nb %*% t(venues.sim[,c(7, 4, 8, 6)]), sdlog = 0.8)

venues.sim$vs.y = venues.y

for(i in 1:ncol(venues.sim)){
  venues.sim[,i] = as.numeric(venues.sim[,i])
}

cells.data = cbind(cells.vc.y, coordinates(cells), cells.vs, cells.vc)
cells.data = cells.data[,-c(9,11)]
cells.data = as.data.frame(cells.data)
names(cells.data) = c("vc.y", "x","y", "smooth.vs", "medium1", "rough.vs", "medium2",
                      "smooth.vc", "rough.vc")


## Count values in each cell for small and large resolution
venues.sp = venues.sim
coordinates(venues.sp) = ~ x + y
venue.cell.over.same = over(venues.sp, cells)
venue.cell.over.small = over(venues.sp, cells.small)
venue.cell.over.large = over(venues.sp, cells.large)
cells.data.small.full = cbind(NA, cells.vc.small)
for(i in 1:length(cells.small)){
  cells.data.small.full[i,1] = sum(venue.cell.over.small == i, na.rm = T)
}
cells.data.small.full = as.data.frame(cells.data.small.full)
names(cells.data.small.full) = c("vc.y", "smooth.vc", "medium1", "rough.vc", "medium2")
cells.data.small.full = cbind(cells.data.small.full, cells.vs.small[,c(1,3)])
names(cells.data.small.full)[c(6,7)] = c("smooth.vs", "rough.vs")

cells.data.large.full = cbind(NA, cells.vc.large)
for(i in 1:length(cells.large)){
  cells.data.large.full[i,1] = sum(venue.cell.over.large == i, na.rm = T)
}
cells.data.large.full = as.data.frame(cells.data.large.full)
names(cells.data.large.full) = c("vc.y", "smooth.vc", "medium1", "rough.vc", "medium2")
cells.data.large.full = cbind(cells.data.large.full, cells.vs.large[,c(1,3)])
names(cells.data.large.full)[c(6,7)] = c("smooth.vs", "rough.vs")




# Simulate sampling -------------------------------------------------------

if(uniform){
  venues.obs = venues.sim[sample(1:nrow(venues.sim), size = round(nrow(venues.sim)/2), replace = F),]
}else{
  ## Use medium1 for the sampling probability
  venues.obs = venues.sim[sample(1:nrow(venues.sim), size = round(nrow(venues.sim)/2), replace = F,
                                 prob = plogis(-2 + 1.5 * venues.sim$medium1)),]
}

scale.fact = nrow(venues.obs) / nrow(venues.sim)

## Count number of observed venues in each cell at all resolutions
venues.obs.sp = venues.obs
coordinates(venues.obs.sp) = ~ x + y

## Get covariates for each venue
venue.cell.over.same = over(venues.obs.sp, cells)
venue.cell.over.small = over(venues.obs.sp, cells.small)
venue.cell.over.large = over(venues.obs.sp, cells.large)

## Count values in each cell
cells.data.same = cbind(NA, cells.vc)
for(i in 1:length(cells)){
  cells.data.same[i,1] = sum(venue.cell.over.same == i, na.rm = T)
}
cells.data.same = as.data.frame(cells.data.same)
names(cells.data.same) = c("vc.y", "smooth.vc", "medium1", "rough.vc", "medium2")
cells.data.same = cbind(cells.data.same, cells.vs[,c(1,3)])
names(cells.data.same)[c(6,7)] = c("smooth.vs", "rough.vs")

cells.data.small = cbind(NA, cells.vc.small)
for(i in 1:length(cells.small)){
  cells.data.small[i,1] = sum(venue.cell.over.small == i, na.rm = T)
}
cells.data.small = as.data.frame(cells.data.small)
names(cells.data.small) = c("vc.y", "smooth.vc", "medium1", "rough.vc", "medium2")
cells.data.small = cbind(cells.data.small, cells.vs.small[,c(1,3)])
names(cells.data.small)[c(6,7)] = c("smooth.vs", "rough.vs")

cells.data.large = cbind(NA, cells.vc.large)
for(i in 1:length(cells.large)){
  cells.data.large[i,1] = sum(venue.cell.over.large == i, na.rm = T)
}
cells.data.large = as.data.frame(cells.data.large)
names(cells.data.large) = c("vc.y", "smooth.vc", "medium1", "rough.vc", "medium2")
cells.data.large = cbind(cells.data.large, cells.vs.large[,c(1,3)])
names(cells.data.large)[c(6,7)] = c("smooth.vs", "rough.vs")

cells.data.same$sampling.prob = scale.fact
cells.data.small$sampling.prob = scale.fact
cells.data.large$sampling.prob = scale.fact


print("Fitting models")
## Use all data
## Fit the venue size model

if(include.cov){
  print("Fitting with all covariates")
  fit.sub.vs = brm(bf(vs.y ~ smooth.vs + medium1 + rough.vs + medium2,
                      hu ~ smooth.vs + medium1 + rough.vs + medium2),
                   data = venues.obs, chains = 1, cores = 1, iter = 1000, refresh = 0, family = "hurdle_lognormal")
}else{
  fit.sub.vs = brm(bf(vs.y ~ smooth.vs + rough.vs + medium2,
                      hu ~ smooth.vs + rough.vs + medium2),
                   data = venues.obs, chains = 1, cores = 1, iter = 1000, refresh = 0, family = "hurdle_lognormal")
}


## Fit the venue count model
if(include.cov){
  fit.sub.vc.same = brm(bf(vc.y ~ smooth.vc + medium1 + rough.vc + medium2 + offset(log(sampling.prob)),
                           zi ~ smooth.vc + medium1 + rough.vc + medium2),
                        data = cells.data.same, chains = 1, cores = 1, iter = 1000, refresh = 0,
                        family = "zero_inflated_negbinomial")
  fit.sub.vc.small = brm(bf(vc.y ~ smooth.vc + medium1 + rough.vc + medium2 + offset(log(sampling.prob)),
                            zi ~ smooth.vc + medium1 + rough.vc + medium2),
                         data = cells.data.small, chains = 1, cores = 1, iter = 1000, refresh = 0,
                         family = "zero_inflated_negbinomial")
  fit.sub.vc.large = brm(bf(vc.y ~ smooth.vc + medium1 + rough.vc + medium2 + offset(log(sampling.prob)),
                            zi ~ smooth.vc + medium1 + rough.vc + medium2),
                         data = cells.data.large, chains = 1, cores = 1, iter = 1000, refresh = 0,
                         family = "zero_inflated_negbinomial")
  
}else{
  fit.sub.vc.same = brm(bf(vc.y ~ smooth.vc + rough.vc + medium2 + offset(log(sampling.prob)),
                           zi ~ smooth.vc + rough.vc + medium2),
                        data = cells.data.same, chains = 1, cores = 1, iter = 1000, refresh = 0,
                        family = "zero_inflated_negbinomial")
  fit.sub.vc.small = brm(bf(vc.y ~ smooth.vc + rough.vc + medium2 + offset(log(sampling.prob)),
                            zi ~ smooth.vc + rough.vc + medium2),
                         data = cells.data.small, chains = 1, cores = 1, iter = 1000, refresh = 0,
                         family = "zero_inflated_negbinomial")
  fit.sub.vc.large = brm(bf(vc.y ~ smooth.vc + rough.vc + medium2 + offset(log(sampling.prob)),
                            zi ~ smooth.vc + rough.vc + medium2),
                         data = cells.data.large, chains = 1, cores = 1, iter = 1000, refresh = 0,
                         family = "zero_inflated_negbinomial")
}








n.cores = 20
print(n.cores)
cl = makeCluster(n.cores)
n.results = 1000
registerDoParallel(cl)
results = foreach(ind=1:n.results, .packages = c("countreg", "spatstat", "sp", "raster", "optimization", "rstan", "brms")) %dopar% {
  library(rstan)
  library(brms)

  cells.vs.y = rep(NA, n.cell)
  cells.vc.y = rep(NA, n.cell)
  
  # Simulate the actual process ---------------------------------------------
  
  ##
  #### Simulate venue locations
  ##
  
  ## Simulate data
  cells.vc.y = rbinom(n.cell, 1, prob = plogis(-1 + betas.vc.bin %*% t(cells.vc))) *
    rnbinom(n.cell, mu = exp(0 + betas.vc.nb %*% t(cells.vc)), size = 10)
  
  ## Now simulate locations inside each cell
  ## Now also grab covariates associated with it
  venues.sim = matrix(NA, nrow = 0, ncol = 2)
  for(i in 1:n.cell){
    if(cells.vc.y[i] > 0){
      cells.xy = spsample(cells[i], cells.vc.y[i], type = "random")
      while(length(cells.xy) != cells.vc.y[i]){
        cells.xy = spsample(cells[i], cells.vc.y[i], type = "random")
      }
      venues.sim = rbind(venues.sim, cbind(coordinates(cells.xy), cells.vc[i,], cells.vs[i,]))
    }
  }
  venues.sim = as.data.frame(venues.sim)
  venues.sim = venues.sim[,-c(8,10)]
  names(venues.sim) = c("x","y", "smooth.vc", "medium1", "rough.vc", "medium2", "smooth.vs", "rough.vs")
  
  ## Simulate venue size
  venues.y = rbinom(nrow(venues.sim), 1,
                    prob = plogis(1 + betas.vs.bin %*% t(venues.sim[,c(7, 4, 8, 6)]))) *
    rlnorm(nrow(venues.sim), 2 + betas.vs.nb %*% t(venues.sim[,c(7, 4, 8, 6)]), sdlog = 0.8)
  ## Don't round these for simulations
  
  venues.sim$vs.y = venues.y
  
  for(i in 1:ncol(venues.sim)){
    venues.sim[,i] = as.numeric(venues.sim[,i])
  }
  
  cells.data = cbind(cells.vc.y, coordinates(cells), cells.vs, cells.vc)
  cells.data = cells.data[,-c(9,11)]
  cells.data = as.data.frame(cells.data)
  names(cells.data) = c("vc.y", "x","y", "smooth.vs", "medium1", "rough.vs", "medium2",
                        "smooth.vc", "rough.vc")
  
  
  ## Count values in each cell for small and large resolution
  venues.sp = venues.sim
  coordinates(venues.sp) = ~ x + y
  venue.cell.over.same = over(venues.sp, cells)
  venue.cell.over.small = over(venues.sp, cells.small)
  venue.cell.over.large = over(venues.sp, cells.large)
  cells.data.small.full = cbind(NA, cells.vc.small)
  for(i in 1:length(cells.small)){
    cells.data.small.full[i,1] = sum(venue.cell.over.small == i, na.rm = T)
  }
  cells.data.small.full = as.data.frame(cells.data.small.full)
  names(cells.data.small.full) = c("vc.y", "smooth.vc", "medium1", "rough.vc", "medium2")
  cells.data.small.full = cbind(cells.data.small.full, cells.vs.small[,c(1,3)])
  names(cells.data.small.full)[c(6,7)] = c("smooth.vs", "rough.vs")
  
  cells.data.large.full = cbind(NA, cells.vc.large)
  for(i in 1:length(cells.large)){
    cells.data.large.full[i,1] = sum(venue.cell.over.large == i, na.rm = T)
  }
  cells.data.large.full = as.data.frame(cells.data.large.full)
  names(cells.data.large.full) = c("vc.y", "smooth.vc", "medium1", "rough.vc", "medium2")
  cells.data.large.full = cbind(cells.data.large.full, cells.vs.large[,c(1,3)])
  names(cells.data.large.full)[c(6,7)] = c("smooth.vs", "rough.vs")
  
  
  
  
  # Simulate sampling -------------------------------------------------------
  
  if(uniform){
    venues.obs = venues.sim[sample(1:nrow(venues.sim), size = round(nrow(venues.sim)/2), replace = F),]
  }else{
    ## Use medium1 for the sampling probability
    venues.obs = venues.sim[sample(1:nrow(venues.sim), size = round(nrow(venues.sim)/2), replace = F,
                                   prob = plogis(-2 + 1.5 * venues.sim$medium1)),]
  }
  
  scale.fact = nrow(venues.obs) / nrow(venues.sim)
  
  ## Count number of observed venues in each cell at all resolutions
  venues.obs.sp = venues.obs
  coordinates(venues.obs.sp) = ~ x + y
  
  ## Get covariates for each venue
  venue.cell.over.same = over(venues.obs.sp, cells)
  venue.cell.over.small = over(venues.obs.sp, cells.small)
  venue.cell.over.large = over(venues.obs.sp, cells.large)
  
  ## Count values in each cell
  cells.data.same = cbind(NA, cells.vc)
  for(i in 1:length(cells)){
    cells.data.same[i,1] = sum(venue.cell.over.same == i, na.rm = T)
  }
  cells.data.same = as.data.frame(cells.data.same)
  names(cells.data.same) = c("vc.y", "smooth.vc", "medium1", "rough.vc", "medium2")
  cells.data.same = cbind(cells.data.same, cells.vs[,c(1,3)])
  names(cells.data.same)[c(6,7)] = c("smooth.vs", "rough.vs")
  
  cells.data.small = cbind(NA, cells.vc.small)
  for(i in 1:length(cells.small)){
    cells.data.small[i,1] = sum(venue.cell.over.small == i, na.rm = T)
  }
  cells.data.small = as.data.frame(cells.data.small)
  names(cells.data.small) = c("vc.y", "smooth.vc", "medium1", "rough.vc", "medium2")
  cells.data.small = cbind(cells.data.small, cells.vs.small[,c(1,3)])
  names(cells.data.small)[c(6,7)] = c("smooth.vs", "rough.vs")
  
  cells.data.large = cbind(NA, cells.vc.large)
  for(i in 1:length(cells.large)){
    cells.data.large[i,1] = sum(venue.cell.over.large == i, na.rm = T)
  }
  cells.data.large = as.data.frame(cells.data.large)
  names(cells.data.large) = c("vc.y", "smooth.vc", "medium1", "rough.vc", "medium2")
  cells.data.large = cbind(cells.data.large, cells.vs.large[,c(1,3)])
  names(cells.data.large)[c(6,7)] = c("smooth.vs", "rough.vs")
  
  cells.data.same$sampling.prob = scale.fact
  cells.data.small$sampling.prob = scale.fact
  cells.data.large$sampling.prob = scale.fact
  
  
  print("Fitting models")
  ## Use all data
  ## Fit the venue size model
  
  fit.sub.vs = update(fit.sub.vs, newdata = venues.obs, recompile = FALSE)
  
  
  ## Fit the venue count model
  fit.sub.vc.same = update(fit.sub.vc.same, newdata = cells.data.same, recompile = FALSE)
  fit.sub.vc.small = update(fit.sub.vc.small, newdata = cells.data.small, recompile = FALSE)
  fit.sub.vc.large = update(fit.sub.vc.large, newdata = cells.data.large, recompile = FALSE)
  
  
  
  
  n.ven = nrow(venues.sim)
  
  cells.data.same$sampling.prob = 1
  cells.data.small$sampling.prob = 1
  cells.data.large$sampling.prob = 1
  
  ven.pred.same.raw = predict(fit.sub.vc.same, type = "response", newdata = cells.data.same, summary = F)
  ven.scale = n.ven / rowSums(ven.pred.same.raw)
  ven.scale = matrix(ven.scale, nrow = nrow(ven.pred.same.raw), ncol = ncol(ven.pred.same.raw))
  ven.pred.same = ven.pred.same.raw * ven.scale
  size.vals.same = predict(fit.sub.vs, type = "response", newdata = cells.data.same, summary = F)
  fsw.pred.total.same = mean(rowSums(ven.pred.same * size.vals.same))
  
  ven.pred.small.raw = predict(fit.sub.vc.small, type = "response", newdata = cells.data.small, summary = F)
  ven.scale = n.ven / rowSums(ven.pred.small.raw)
  ven.scale = matrix(ven.scale, nrow = nrow(ven.pred.small.raw), ncol = ncol(ven.pred.small.raw))
  ven.pred.small = ven.pred.small.raw * ven.scale
  size.vals.small = predict(fit.sub.vs, type = "response", newdata = cells.data.small, summary = F)
  fsw.pred.total.small = mean(rowSums(ven.pred.small * size.vals.small))
  
  ven.pred.large.raw = predict(fit.sub.vc.large, type = "response", newdata = cells.data.large, summary = F)
  ven.scale = n.ven / rowSums(ven.pred.large.raw)
  ven.scale = matrix(ven.scale, nrow = nrow(ven.pred.large.raw), ncol = ncol(ven.pred.large.raw))
  ven.pred.large = ven.pred.large.raw * ven.scale
  size.vals.large = predict(fit.sub.vs, type = "response", newdata = cells.data.large, summary = F)
  fsw.pred.total.large = mean(rowSums(ven.pred.large * size.vals.large))
  
  
  
  return.list = list(fsw.total = sum(venues.sim$vs.y), ven.total = n.ven,
                     total.fsw.same = fsw.pred.total.same, total.fsw.small = fsw.pred.total.small,
                     total.fsw.large = fsw.pred.total.large,
                     ven.same.truth = cells.data$vc.y, ven.small.truth = cells.data.small.full$vc.y,
                     ven.large.truth = cells.data.large.full$vc.y,
                     ven.same.est = colMeans(ven.pred.same), ven.small.est = colMeans(ven.pred.small),
                     ven.large.est = colMeans(ven.pred.large),
                     true.size.coef = c(2, -1, betas.vs.nb, -betas.vs.bin),
                     sub.size.coef = fixef(fit.sub.vs)[,1],
                     true.count.coef = c(0, 1, betas.vc.nb, -betas.vc.bin),
                     sub.count.coef.small = fixef(fit.sub.vc.small)[,1],
                     sub.count.coef.same = fixef(fit.sub.vc.same)[,1],
                     sub.count.coef.large = fixef(fit.sub.vc.large)[,1])
}
stopCluster(cl)

if(uniform == TRUE & include.cov == TRUE){
  save(results, file = "Sim_thin_unif_included.RData")
}else if(uniform == TRUE & include.cov == FALSE){
  save(results, file = "Sim_thin_unif_notincluded.RData")
}else if(uniform == FALSE & include.cov == TRUE){
  save(results, file = "Sim_thin_nonunif_included.RData")
}else if(uniform == FALSE & include.cov == FALSE){
  save(results, file = "Sim_thin_nonunif_notincluded.RData")
}




load("Sim_thin_unif_included.RData")
results.unif.included = results
load("Sim_thin_unif_notincluded.RData")
results.unif.notincluded = results
load("Sim_thin_nonunif_notincluded.RData")
results.nonunif.notincluded = results
load("Sim_thin_nonunif_included.RData")
results.nonunif.included = results





## Look at FSW size against truth
total.fsw.unif.included.truth = total.fsw.unif.notincluded.truth =
  total.fsw.nonunif.included.truth = total.fsw.nonunif.notincluded.truth = rep(NA, 1000)
total.fsw.unif.included.small = total.fsw.unif.included.same = total.fsw.unif.included.large = rep(NA, 1000)
total.fsw.unif.notincluded.small = total.fsw.unif.notincluded.same = total.fsw.unif.notincluded.large = rep(NA, 1000)
total.fsw.nonunif.notincluded.small = total.fsw.nonunif.notincluded.same = total.fsw.nonunif.notincluded.large = rep(NA, 1000)
total.fsw.nonunif.included.small = total.fsw.nonunif.included.same = total.fsw.nonunif.included.large = rep(NA, 1000)

total.fsw.unif.included.small.full = total.fsw.unif.included.same.full = total.fsw.unif.included.large.full = rep(NA, 1000)
total.fsw.unif.notincluded.small.full = total.fsw.unif.notincluded.same.full = total.fsw.unif.notincluded.large.full = rep(NA, 1000)
total.fsw.nonunif.notincluded.small.full = total.fsw.nonunif.notincluded.same.full = total.fsw.nonunif.notincluded.large.full = rep(NA, 1000)
total.fsw.nonunif.included.small.full = total.fsw.nonunif.included.same.full = total.fsw.nonunif.included.large.full = rep(NA, 1000)

for(i in 1:1000){
  ####################### FSW
  ## Get truth
  total.fsw.unif.included.truth[i] = results.unif.included[[i]]$fsw.total
  total.fsw.unif.notincluded.truth[i] = results.unif.notincluded[[i]]$fsw.total
  total.fsw.nonunif.notincluded.truth[i] = results.nonunif.notincluded[[i]]$fsw.total
  total.fsw.nonunif.included.truth[i] = results.nonunif.included[[i]]$fsw.total

  ## Get presence-only models
  ## Get unif included
  total.fsw.unif.included.same[i] = results.unif.included[[i]]$total.fsw.same
  total.fsw.unif.included.small[i] = results.unif.included[[i]]$total.fsw.small
  total.fsw.unif.included.large[i] = results.unif.included[[i]]$total.fsw.large
  
  ## Get unif notincluded
  total.fsw.unif.notincluded.same[i] = results.unif.notincluded[[i]]$total.fsw.same
  total.fsw.unif.notincluded.small[i] = results.unif.notincluded[[i]]$total.fsw.small
  total.fsw.unif.notincluded.large[i] = results.unif.notincluded[[i]]$total.fsw.large
  
  ## Get nonunif notincluded
  total.fsw.nonunif.notincluded.same[i] = results.nonunif.notincluded[[i]]$total.fsw.same
  total.fsw.nonunif.notincluded.small[i] = results.nonunif.notincluded[[i]]$total.fsw.small
  total.fsw.nonunif.notincluded.large[i] = results.nonunif.notincluded[[i]]$total.fsw.large
  
  ## Get nonunif included
  total.fsw.nonunif.included.same[i] = results.nonunif.included[[i]]$total.fsw.same
  total.fsw.nonunif.included.small[i] = results.nonunif.included[[i]]$total.fsw.small
  total.fsw.nonunif.included.large[i] = results.nonunif.included[[i]]$total.fsw.large
}


## Unif included vs truth
df.unif.included.vs.truth = data.frame(resid =
                         c((total.fsw.unif.included.same - total.fsw.unif.included.truth) / total.fsw.unif.included.truth * 100,
                           (total.fsw.unif.included.small - total.fsw.unif.included.truth) / total.fsw.unif.included.truth * 100,
                           (total.fsw.unif.included.large - total.fsw.unif.included.truth) / total.fsw.unif.included.truth * 100),
                           size = c(rep(c("same", "small", "large"), each = 1000)))

gg.truth.unif.included = ggplot(df.unif.included.vs.truth) +
  stat_density(mapping = aes(x = resid, group = size, linetype = size), geom = "line", position = "identity", adjust = 1, size = 1.3) +
  scale_linetype_manual(values = c("solid", "dotted", "twodash"), breaks = c("same", "small", "large")) + 
  ggtitle("Uniform Sampling \nwith All Covariates") + xlab("Percent Error") +
  geom_vline(xintercept = 0, col = "red") + xlim(-40, 40) +
  guides(linetype = guide_legend(title = "Cell Size", keywidth = 5, override.aes = list(size = 1))) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=18),
        legend.text=element_text(size = 13), legend.title=element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold")) +
  ylab("Density")
gg.truth.unif.included

df.same.unif.included = subset(df.unif.included.vs.truth, size == "same")
df.small.unif.included = subset(df.unif.included.vs.truth, size == "small")
df.large.unif.included = subset(df.unif.included.vs.truth, size == "large")
df.same.unif.included$assump = "unif.included"
df.small.unif.included$assump = "unif.included"
df.large.unif.included$assump = "unif.included"


## Unif notincluded vs truth
df.unif.notincluded.vs.truth = data.frame(resid = c((total.fsw.unif.notincluded.same - total.fsw.unif.notincluded.truth) / total.fsw.unif.included.truth * 100,
                                                    (total.fsw.unif.notincluded.small - total.fsw.unif.notincluded.truth) / total.fsw.unif.included.truth * 100,
                                                    (total.fsw.unif.notincluded.large - total.fsw.unif.notincluded.truth) / total.fsw.unif.included.truth * 100),
                                          size = c(rep(c("same", "small", "large"), each = 1000)))

gg.truth.unif.notincluded = ggplot(df.unif.notincluded.vs.truth) +
  stat_density(mapping = aes(x = resid, group = size, linetype = size), geom = "line", position = "identity", adjust = 1, size = 1.3) +
  scale_linetype_manual(values = c("solid", "dotted", "twodash"), breaks = c("same", "small", "large")) + 
  ggtitle("Uniform Sampling \nwith Missing Covariate") + xlab("Percent Error") +
  geom_vline(xintercept = 0, col = "red") + xlim(-40, 40) +
  guides(linetype = guide_legend(keywidth = 3, override.aes = list(size = 1))) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=18),
        legend.text=element_text(size = 13), legend.title=element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold")) +
  ylab("Density")
gg.truth.unif.notincluded


df.same.unif.notincluded = subset(df.unif.notincluded.vs.truth, size == "same")
df.small.unif.notincluded = subset(df.unif.notincluded.vs.truth, size == "small")
df.large.unif.notincluded = subset(df.unif.notincluded.vs.truth, size == "large")
df.same.unif.notincluded$assump = "unif.notincluded"
df.small.unif.notincluded$assump = "unif.notincluded"
df.large.unif.notincluded$assump = "unif.notincluded"



## nonunif notincluded vs truth
df.nonunif.notincluded.vs.truth = data.frame(resid = c((total.fsw.nonunif.notincluded.same - total.fsw.nonunif.notincluded.truth) / total.fsw.unif.included.truth * 100,
                                                       (total.fsw.nonunif.notincluded.small - total.fsw.nonunif.notincluded.truth) / total.fsw.unif.included.truth * 100,
                                                       (total.fsw.nonunif.notincluded.large - total.fsw.nonunif.notincluded.truth) / total.fsw.unif.included.truth * 100),
                                             size = c(rep(c("same", "small", "large"), each = 1000)))

gg.truth.nonunif.notincluded = ggplot(df.nonunif.notincluded.vs.truth) +
  stat_density(mapping = aes(x = resid, group = size, linetype = size), geom = "line", position = "identity", adjust = 1, size = 1.3) +
  scale_linetype_manual(values = c("solid", "dotted", "twodash"), breaks = c("same", "small", "large")) + 
  ggtitle("Nonuniform Sampling \nwith Missing Covariate") + xlab("Percent Error") +
  geom_vline(xintercept = 0, col = "red") + xlim(-40, 40) +
  guides(linetype = guide_legend(keywidth = 3, override.aes = list(size = 1))) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=18),
        legend.text=element_text(size = 13), legend.title=element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold")) +
  ylab("Density")
gg.truth.nonunif.notincluded


df.same.nonunif.notincluded = subset(df.nonunif.notincluded.vs.truth, size == "same")
df.small.nonunif.notincluded = subset(df.nonunif.notincluded.vs.truth, size == "small")
df.large.nonunif.notincluded = subset(df.nonunif.notincluded.vs.truth, size == "large")
df.same.nonunif.notincluded$assump = "nonunif.notincluded"
df.small.nonunif.notincluded$assump = "nonunif.notincluded"
df.large.nonunif.notincluded$assump = "nonunif.notincluded"



## nonunif included vs truth
df.nonunif.included.vs.truth = data.frame(resid = c((total.fsw.nonunif.included.same - total.fsw.nonunif.included.truth) / total.fsw.unif.included.truth * 100,
                                                    (total.fsw.nonunif.included.small - total.fsw.nonunif.included.truth) / total.fsw.unif.included.truth * 100,
                                                    (total.fsw.nonunif.included.large - total.fsw.nonunif.included.truth) / total.fsw.unif.included.truth * 100),
                                          size = c(rep(c("same", "small", "large"), each = 1000)))

gg.truth.nonunif.included = ggplot(df.nonunif.included.vs.truth) +
  stat_density(mapping = aes(x = resid, group = size, linetype = size), geom = "line", position = "identity", adjust = 1, size = 1.3) +
  scale_linetype_manual(values = c("solid", "dotted", "twodash"), breaks = c("same", "small", "large")) + 
  ggtitle("Nonuniform Sampling \nwith All Covariates") + xlab("Percent Error") +
  geom_vline(xintercept = 0, col = "red") + xlim(-40, 40) +
  guides(linetype = guide_legend(keywidth = 3, override.aes = list(size = 1))) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=18),
        legend.text=element_text(size = 13), legend.title=element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold")) +
  ylab("Density")
gg.truth.nonunif.included

df.same.nonunif.included = subset(df.nonunif.included.vs.truth, size == "same")
df.small.nonunif.included = subset(df.nonunif.included.vs.truth, size == "small")
df.large.nonunif.included = subset(df.nonunif.included.vs.truth, size == "large")
df.same.nonunif.included$assump = "nonunif.included"
df.small.nonunif.included$assump = "nonunif.included"
df.large.nonunif.included$assump = "nonunif.included"

## Save at 2200 x 850
grid.arrange(gg.truth.unif.included + theme(legend.position="NA"),
             gg.truth.nonunif.included + theme(legend.position="NA") + ylab(""),
             gg.truth.unif.notincluded + theme(legend.position="NA") + ylab(""),
             gg.truth.nonunif.notincluded + theme(legend.position="NA") + ylab(""),
             widths = c(3, 3, 3, 3))
## Save at 2200 x 850

gg1 = gg.truth.unif.included + theme(legend.position="NA")
gg2 = gg.truth.nonunif.included + theme(legend.position="NA") + ylab("")
gg3 = gg.truth.unif.notincluded + theme(legend.position="NA") + ylab("")
gg4 = gg.truth.nonunif.notincluded + theme(legend.position="NA") + ylab("")
ggsave("cell_sim_a.png", plot = gg1, width = 6, height = 6, units = "in")
ggsave("cell_sim_b.png", plot = gg2, width = 6, height = 6, units = "in")
ggsave("cell_sim_c.png", plot = gg3, width = 6, height = 6, units = "in")
ggsave("cell_sim_d.png", plot = gg4, width = 6, height = 6, units = "in")





dat.same.1 = subset(df.unif.included.vs.truth, size == "same")
dat.same.2 = subset(df.unif.notincluded.vs.truth, size == "same")
dat.same.3 = subset(df.nonunif.notincluded.vs.truth, size == "same")
dat.same.4 = subset(df.nonunif.included.vs.truth, size == "same")

quantile(dat.same.1$resid, probs = c(0.025, 0.975))
quantile(dat.same.2$resid, probs = c(0.025, 0.975))
quantile(dat.same.3$resid, probs = c(0.025, 0.975))
quantile(dat.same.4$resid, probs = c(0.025, 0.975))

mean(dat.same.1$resid)
mean(dat.same.2$resid)
mean(dat.same.3$resid)
mean(dat.same.4$resid)

dat.small.1 = subset(df.unif.included.vs.truth, size == "small")
dat.large.1 = subset(df.unif.included.vs.truth, size == "large")
dat.small.2 = subset(df.unif.notincluded.vs.truth, size == "small")
dat.large.2 = subset(df.unif.notincluded.vs.truth, size == "large")
dat.small.3 = subset(df.nonunif.notincluded.vs.truth, size == "small")
dat.large.3 = subset(df.nonunif.notincluded.vs.truth, size == "large")
dat.small.4 = subset(df.nonunif.included.vs.truth, size == "small")
dat.large.4 = subset(df.nonunif.included.vs.truth, size == "large")

mean(dat.same.1$resid)
mean(dat.small.1$resid)
mean(dat.large.1$resid)

mean(dat.same.2$resid)
mean(dat.small.2$resid)
mean(dat.large.2$resid)

mean(dat.same.3$resid)
mean(dat.small.3$resid)
mean(dat.large.3$resid)

mean(dat.same.4$resid)
mean(dat.small.4$resid)
mean(dat.large.4$resid)

quantile(dat.same.1$resid, probs = c(0.025, 0.975))
quantile(dat.small.1$resid, probs = c(0.025, 0.975))
quantile(dat.large.1$resid, probs = c(0.025, 0.975))


