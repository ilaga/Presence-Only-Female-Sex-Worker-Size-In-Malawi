

library(rstan)
library(brms)
library(bayesplot) ## Only for plotting
library(ggplot2) ## Only for plotting
library(fields) ## Only for plotting
library(raster) ## Only for plotting


# Initial Data ------------------------------------------------------------

myscale<-function(x){(x-mean(x, na.rm = T))/sd(x, na.rm = T)}

dat.venue = readRDS("Malawi_Venues_nightlight.rds")
dat.cell = readRDS("Malawi_Cells_worldPop.rds")
dat.venue$District = as.character(dat.venue$b5a.y)
dat.venue$District[dat.venue$District == "Chikhwawa"] = "Chikwawa"

dat.venue = dat.venue[,c(3,4,103,111:155)]
dat.venue$fswestimate[dat.venue$fswestimate > 100] = 100


## Now transform the worldPop and nightlight values
dat.venue$worldPop = log(dat.venue$worldPop)
dat.venue$nightlight = log(dat.venue$nightlight)
dat.cell$worldPop = log(dat.cell$worldPop)
dat.cell$nightlight = log(dat.cell$nightlight)

which.inf = which(dat.cell$worldPop == -Inf)
dat.cell = dat.cell[-which.inf,]

########## Now scale everything
dat.cell.or = dat.cell
df = dat.cell[,-c(1,2,3,26,27)]
cov.names = names(df)
meancov<-apply(df, 2, mean)
sdcov<-apply(df, 2, sd)
df<-cbind(dat.cell[,c("fswave","District","numplace","x","y")],
          apply(df, 2, myscale))
dat.cell = df
dat.cell = cbind(dat.cell)

## Scale place.dat
dat.venue$mv012<-(dat.venue$mv012-meancov["mv012"])/sdcov["mv012"]
dat.venue$mv035<-(dat.venue$mv035-meancov["mv035"])/sdcov["mv035"]
dat.venue$mv167<-(dat.venue$mv167-meancov["mv167"])/sdcov["mv167"]
dat.venue$mv191<-(dat.venue$mv191-meancov["mv191"])/sdcov["mv191"]
dat.venue$mv201<-(dat.venue$mv201-meancov["mv201"])/sdcov["mv201"]
dat.venue$mv245<-(dat.venue$mv245-meancov["mv245"])/sdcov["mv245"]
dat.venue$mv531<-(dat.venue$mv531-meancov["mv531"])/sdcov["mv531"]
dat.venue$mv766a<-(dat.venue$mv766a-meancov["mv766a"])/sdcov["mv766a"]
dat.venue$hv204<-(dat.venue$hv204-meancov["hv204"])/sdcov["hv204"]
dat.venue$mv168<-(dat.venue$mv168-meancov["mv168"])/sdcov["mv168"]
dat.venue$mv791<-(dat.venue$mv791-meancov["mv791"])/sdcov["mv791"]
dat.venue$mv793<-(dat.venue$mv793-meancov["mv793"])/sdcov["mv793"]
dat.venue$v155<-(dat.venue$v155-meancov["v155"])/sdcov["v155"]
dat.venue$v168<-(dat.venue$v168-meancov["v168"])/sdcov["v168"]
dat.venue$hivpos<-(dat.venue$hivpos-meancov["hivpos"])/sdcov["hivpos"]
dat.venue$worldPop<-(dat.venue$worldPop-meancov["worldPop"])/sdcov["worldPop"]
dat.venue$nightlight<-(dat.venue$nightlight-meancov["nightlight"])/sdcov["nightlight"]
dat.venue$aridity<-(dat.venue$aridity-meancov["aridity"])/sdcov["aridity"]
dat.venue$built<-(dat.venue$built-meancov["built"])/sdcov["built"]
dat.venue$borders<-(dat.venue$borders-meancov["borders"])/sdcov["borders"]
dat.venue$water<-(dat.venue$water-meancov["water"])/sdcov["water"]
dat.venue$travel<-(dat.venue$travel-meancov["travel"])/sdcov["travel"]


######################################################
## Make plots for WorldPop and Nightlight
worldpop_map = ggplot() +
  geom_tile(data = dat.cell, mapping = aes(x = x, y = y, fill = worldPop), show.legend = T) +
  labs(x = "Longitude", y = "Latitude", fill = "WorldPop") +
  theme(text = element_text(size = 15), legend.title = element_text(size = 15)) +
  scale_fill_gradientn(colours = tim.colors()) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14)) +
  coord_equal()

nightlight_map = ggplot() +
  geom_tile(data = dat.cell, mapping = aes(x = x, y = y, fill = nightlight), show.legend = T) +
  labs(x = "Longitude", y = "Latitude", fill = "Nightlight") +
  theme(text = element_text(size = 15), legend.title = element_text(size = 15)) +
  scale_fill_gradientn(colours = tim.colors()) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))+
  coord_equal()

## Save these at 500 x 700
ggsave("worldpop_map.png", plot = worldpop_map, width = 5, height = 7, units = "in")
ggsave("nightlight_map.png", plot = nightlight_map, width = 5, height = 7, units = "in")
######################################################



## Percent visited in each district
dist.name = c("Lilongwe", "Blantyre","Mangochi", "Machinga", "Mzuzu", "Zomba",
              "Balaka", "Chikwawa", "Dedza", "Dowa", "Karonga",
              "Kasungu", "Mchinji", "Mwanza", "Mzimba", "Neno", "Nkhata Bay",
              "Nkhotakota", "Ntcheu", "Rumphi", "Salima")
dist.excluded = c(14,21,3,11,19,14,37,21,45,8,48,0,18,0,5,9,11,0,1,18,24)
dist.attempt.yes = c(1338,879,353,145,327,319,211,259,219,213,193,176,117,102,246,81,160,243,196,157,211)
dist.attempt.no = c(350,810,83,130,1,109,4,0,6,13,1,25,0,0,55,0,1,0,0,0,18)
dist.oper = c(793,565,179,131,162,188,153,191,175,167,150,140,97,77,207,61,122,164,150,112,178)
dist.est.op = c(993,1104,221,248,162,252,160,197,180,178,160,153,102,78,251,65,124,179,148,113,192)
have.visited = rep(NA, length(dist.name))
j = 1
for(i in dist.name){
  have.visited[j] = sum(dat.cell$numplace[dat.cell$District == i], na.rm = T)
  j = j + 1
}
dat.venues.info = data.frame(Name = dist.name, Excluded = dist.excluded, Attempted.Yes = dist.attempt.yes,
                             Attempted.No = dist.attempt.no, Operational = dist.oper,
                             Est.Op = dist.est.op, Have.visit = have.visited)
dat.venues.info$sampling.prob = dat.venues.info$Have.visit / dat.venues.info$Est.Op
dat.venues.info
dat.cell$sampling.prob = NA
for(dist in dist.name){
  rel.dist = which(dat.cell$District == dist)
  sampling.prob = dat.venues.info$sampling.prob[which(dat.venues.info$Name == dist)]
  if(length(sampling.prob) == 1){
    dat.cell$sampling.prob[rel.dist] = sampling.prob
  }
}



dist.place = c("Balaka","Chikwawa","Dedza","Dowa","Karonga","Kasungu","Mchinji","Mwanza","Mzimba","Neno","Nkhata Bay",
               "Nkhotakota","Ntcheu","Rumphi","Salima")
dat.cell.all = dat.cell
dat.cell = subset(dat.cell, District %in% dist.place)

## Remove Mzimba entirely from count modeling
## Only need to remove it from dat.cell
## dat.venue still uses Mzimba
dat.cell = subset(dat.cell, District != "Mzimba")


dat.venue$District = factor(dat.venue$District)
dat.cell$District = factor(dat.cell$District)
n.dist = length(levels(dat.cell$District))

# Now start the model fitting ---------------------------------------------
count.model = brm(bf(numplace ~ v155 + built +
                       worldPop + (1|District) + offset(log(sampling.prob)),
                     zi ~ mv167 +
                       v155 + built +
                       worldPop + nightlight + (1|District) + gp(x, y, k = 5, c = 5/4)),
                  data = dat.cell, chains = 4, cores = 4, iter = 5000,
                  control = list(adapt_delta = 0.97, max_treedepth = 12),
                  family = "zero_inflated_negbinomial")
summary(count.model)

## Look at posterior pairwise scatterplot
pairs(count.model, par = c("b", "b_zi"))
## Look at graphical posterior checks
pp_check(count.model, type = "xyz")
pp_check(count.model, type = "rootogram", style = "hanging", transform = "round") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14)) +
  xlim(-1, 25) + theme_grey()
pp_check(count.model, type = "scatter_avg") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))
pp_check(count.model, type = "stat_2d") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))
pp_check(count.model, type = "error_scatter_avg_vs_x", x = "worldPop") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))
pp_check(count.model, type = "error_scatter_avg_vs_x", x = "built") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))
pp_check(count.model, type = "error_scatter_avg_vs_x", x = "nightlight") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))
pp_check(count.model, type = "error_scatter_avg_vs_x", x = "v155") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))
pp_check(count.model, type = "error_scatter_avg_vs_x", x = "mv791") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))
pp_check(count.model, type = "error_scatter_avg_vs_x", x = "mv167") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))


## Posterior predictive p-values
## mean, sd, maximum, proportion of zeros

## Calculate posterior predictive p-value
tmp = posterior_predict(count.model, newdata = dat.cell, nsamples = 1000, summary = F)
mean.vec = rowMeans(tmp)
pos.mean.vec = pos.sd.vec = rep(NA, nrow(tmp))
for(i in 1:nrow(tmp)){
  pos.ind = which(tmp[i,] > 0)
  pos.mean.vec[i] = mean(tmp[i,pos.ind])
  pos.sd.vec[i] = sd(tmp[i,pos.ind])
}
sd.vec = apply(tmp, 1, sd)
max.vec = apply(tmp, 1, max)
prob.zero.vec = apply(tmp, 1, function(x){sum(x == 0)})
overdisp.vec = sd.vec / mean.vec
mean.obs = mean(dat.cell$numplace)
sd.obs = sd(dat.cell$numplace)
pos.mean.obs = mean(dat.cell$numplace[which(dat.cell$numplace > 0)])
pos.sd.obs = sd(dat.cell$numplace[which(dat.cell$numplace > 0)])
max.obs = max(dat.cell$numplace)
prob.zero.obs = sum(dat.cell$numplace == 0)
overdisp.obs = sd.obs / mean.obs

## Mean
mean(mean.vec > mean.obs)
mean(pos.mean.vec > pos.mean.obs)
## SD
mean(sd.vec > sd.obs)
mean(pos.sd.vec > pos.sd.obs)
## Max
mean(max.vec > max.obs)
## prob.zero
mean(prob.zero.vec > prob.zero.obs)
## Overdisp
mean(overdisp.vec > overdisp.obs)



## Diagnostic plots for paper
color_scheme_set(scheme = "darkgray")
count_root = pp_check(count.model, type = "rootogram", style = "hanging", nsamples = 1000) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14)) +
  xlim(-1, 30)

count_root

count_scatter = pp_check(count.model, type = "scatter_avg", nsamples = 1000) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))

count_scatter




# Venue size ----------------------------------------------------------
dat.venue.exist = subset(dat.venue, !is.na(fswestimate))
dat.venue.exist = subset(dat.venue.exist, fswestimate < 100)

size.model = brm(bf(fswestimate ~ mv167 +
                      hivpos + worldPop + nightlight + (1|District),
                    hu ~ mv201 +
                      built + worldPop +
                      nightlight + (1|District)),
                 data = dat.venue.exist, chains = 7, cores = 7, iter = 10000, family = "hurdle_lognormal")


size_root = pp_check(size.model, type = "rootogram", style = "hanging", transform = "round", nsamples = 1000) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14)) +
  xlim(-1, 81)

size_root

size_scatter = pp_check(size.model, type = "scatter_avg", nsamples = 1000) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))

size_scatter

grid.diag = grid.arrange(size_scatter, size_root, count_scatter, count_root, ncol = 4)

## Save all diagnostic plots
ggsave("count_rootogram.png", plot = count_root, width = 6, height = 6, units = "in")
ggsave("count_scatter.png", plot = count_scatter, width = 6, height = 6, units = "in")
ggsave("size_rootogram.png", plot = size_root, width = 6, height = 6, units = "in")
ggsave("size_scatter.png", plot = size_scatter, width = 6, height = 6, units = "in")


## Calculate posterior predictive p-values
tmp = posterior_predict(size.model, newdata = dat.venue.exist, nsamples = 4000, summary = F)
mean.vec = rowMeans(tmp)
pos.mean.vec = pos.sd.vec = rep(NA, nrow(tmp))
for(i in 1:nrow(tmp)){
  pos.ind = which(tmp[i,] > 0)
  pos.mean.vec[i] = mean(tmp[i,pos.ind])
  pos.sd.vec[i] = sd(tmp[i,pos.ind])
}
sd.vec = apply(tmp, 1, sd)
max.vec = apply(tmp, 1, max)
prob.zero.vec = apply(tmp, 1, function(x){sum(x == 0)})
overdisp.vec = sd.vec / mean.vec
mean.obs = mean(dat.venue.exist$fswestimate)
sd.obs = sd(dat.venue.exist$fswestimate)
pos.mean.obs = mean(dat.venue.exist$fswestimate[which(dat.venue.exist$fswestimate > 0)])
pos.sd.obs = sd(dat.venue.exist$fswestimate[which(dat.venue.exist$fswestimate > 0)])
max.obs = max(dat.venue.exist$fswestimate)
prob.zero.obs = sum(dat.venue.exist$fswestimate == 0)
overdisp.obs = sd.obs / mean.obs

## Mean
mean(mean.vec > mean.obs)
mean(pos.mean.vec > pos.mean.obs)
## SD
mean(sd.vec > sd.obs)
mean(pos.sd.vec > pos.sd.obs)
## Max
mean(max.vec > max.obs)
## prob.zero
mean(prob.zero.vec > prob.zero.obs)
## Overdisp
mean(overdisp.vec > overdisp.obs)


# Estimate for all cells and districts ------------------------------------
## (1) Estimate for all cells
## (2) For each district, scale values in each cell

## Create dist.pred matrix with columns:
## Name - Average scaling factor - PLACE Venue Size Estimate - Original Venue Size Estimate - Original Venue lower bound -
## Original Venue upper bound - Original FSW Size Estimate - Original FSW Size LB - Original FSW Size UB
## Scaled FSW Size Estimate - Scaled LB - Scaled UB - PLACE FSW

dist.pred = data.frame(Name = c("Dedza", "Dowa", "Kasungu", "Lilongwe", "Mchinji", "Nkhotakota",
                                "Ntcheu", "Salima", "Ntchisi", "Nkhata Bay", "Karonga", "Rumphi",
                                "Mzimba", "Chitipa", "Balaka", "Chikwawa", "Mwanza", "Neno",
                                "Blantyre", "Mangochi", "Machinga", "Zomba", "Chiradzulu", "Mulanje",
                                "Nsanje", "Thyolo", "Phalombe", "Likoma"),
                       Scaling.Factor = rep(NA, 28),
                       PLACE.Venue = rep(NA, 28),
                       OG.Venue.est = rep(NA, 28), OG.Venue.lower = rep(NA, 28), OG.Venue.upper = rep(NA, 28),
                       OG.FSW.est = rep(NA, 28), OG.FSW.lower = rep(NA, 28), OG.FSW.upper = rep(NA, 28),
                       Scaled.FSW.est = rep(NA, 28), Scaled.FSW.lower = rep(NA, 28), Scaled.FSW.upper = rep(NA, 28),
                       PLACE.FSW = c(500,600,1500,7000,600,1300,1100,2000,400,700,800,400,1700,400,
                                     300,900,800,300,6200,900,1000,1800,200,1500,200,1600,300,NA))
dist.pred$PLACE.Venue = dat.venues.info$Est.Op[match(dist.pred$Name, dat.venues.info$Name)]
dist.pred$PLACE = 3
dist.pred$PLACE[dist.pred$Name %in% levels(dat.cell$District)] = 2
dist.pred$PLACE[!(dist.pred$Name %in% levels(dat.cell$District)) & !is.na(dist.pred$PLACE.Venue)] = 1


## Predict to all cells
dat.cell.all$sampling.prob = 1
count.model.pred = fitted(count.model, newdata = dat.cell.all, allow_new_levels = TRUE, summary = F, nsamples = 5000)
size.model.pred = fitted(size.model, newdata = dat.cell.all, allow_new_levels = TRUE, summary = F, nsamples = 5000)
fsw.pred = count.model.pred * size.model.pred

## Estimate 
scaling.factors.df = matrix(NA, nrow = 28, ncol = nrow(fsw.pred))
for(i in 1:28){
  which.cells = which(dat.cell.all$District == dist.pred$Name[i])
  dist.count.pred = count.model.pred[,which.cells]
  dist.size.pred = size.model.pred[,which.cells]
  dist.venue.total.pred = sum(colMeans(dist.count.pred))
  if(!is.na(dist.pred$PLACE.Venue[i])){
    ## Scale each posterior sample
    scale.fact = dist.pred$PLACE.Venue[i] / rowSums(dist.count.pred)
    scaling.factors.df[i,] = scale.fact
    scale.fact = matrix(scale.fact, nrow = nrow(count.model.pred), ncol = ncol(dist.count.pred))
    
    dist.pred$Scaling.Factor[i] = mean(scale.fact)
    interval = quantile(rowSums(dist.count.pred))
    dist.pred$OG.Venue.est[i] = dist.venue.total.pred
    dist.pred$OG.Venue.lower[i] = interval[1]
    dist.pred$OG.Venue.upper[i] = interval[2]
    
    interval = quantile(rowSums(dist.count.pred * dist.size.pred), probs = c(0.05, 0.95))
    dist.pred$OG.FSW.est[i] = sum(colMeans(dist.count.pred * dist.size.pred))
    dist.pred$OG.FSW.lower[i] = interval[1]
    dist.pred$OG.FSW.upper[i] = interval[2]
    
    interval = quantile(rowSums(dist.count.pred * scale.fact * dist.size.pred), probs = c(0.05, 0.95))
    dist.pred$Scaled.FSW.est[i] = sum(colMeans(dist.count.pred * scale.fact * dist.size.pred))
    dist.pred$Scaled.FSW.lower[i] = interval[1]
    dist.pred$Scaled.FSW.upper[i] = interval[2]
    
    fsw.pred[,which.cells] = fsw.pred[,which.cells] * scale.fact
  }else{
    interval = quantile(rowSums(dist.count.pred))
    dist.pred$OG.Venue.est[i] = dist.venue.total.pred
    dist.pred$OG.Venue.lower[i] = interval[1]
    dist.pred$OG.Venue.upper[i] = interval[2]
    
    interval = quantile(rowSums(dist.count.pred * dist.size.pred), probs = c(0.05, 0.95))
    dist.pred$OG.FSW.est[i] = sum(colMeans(dist.count.pred * dist.size.pred))
    dist.pred$OG.FSW.lower[i] = interval[1]
    dist.pred$OG.FSW.upper[i] = interval[2]
  }
  print(c(i,"/",28))
}
dist.pred

scaling.plot.df = data.frame(Name = dist.pred$Name, mean = rowMeans(scaling.factors.df),
                             lower = apply(scaling.factors.df, 1, quantile, probs = c(0.025), na.rm = T),
                             upper = apply(scaling.factors.df, 1, quantile, probs = c(0.975), na.rm = T),
                             PLACE = dist.pred$PLACE)
scaling.plot.df = scaling.plot.df[-28,]
scaling.plot.df = data.frame(scaling.plot.df)
scaling.plot.df$Name = factor(scaling.plot.df$Name, levels = scaling.plot.df$Name[order(colMeans(dist.mat), decreasing = FALSE)])
scaling.plot.df = subset(scaling.plot.df, !is.na(mean))
scaling.plot.df$PLACE = factor(scaling.plot.df$PLACE)

gg.scaling = ggplot(scaling.plot.df,aes(x = Name, y = mean, shape = PLACE)) + geom_point(size = 3) +
  scale_shape_manual(values = shapes, name = "Place Study", labels = c("1", "2", "None")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) +
  geom_abline(slope = 0, intercept = 1) +
  ylab("Scaling Factor") + xlab("District") +
  theme(axis.title = element_text(size = 15)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)
gg.scaling

ggsave("scaling_factors.png", plot = gg.scaling, width = 9, height = 6, units = "in")



cell.est = data.frame(est = colMeans(fsw.pred), var = apply(fsw.pred, 2, var))

ggplot() + geom_point(aes(x = dat.cell.all$x, y= dat.cell.all$y, col = log(colMeans(fsw.pred)))) +
  coord_equal() +
  scale_color_gradientn(colours = tim.colors())


## bayesplot plots
## Need each column to be a new parameter (District)
dist.mat = matrix(NA, nrow = nrow(fsw.pred), ncol = 28)
place.points = matrix(NA, nrow = nrow(fsw.pred), ncol = 28)
for(i in 1:28){
  which.cell = which(dat.cell.all$District == dist.pred$Name[i])
  dist.mat[, i] = rowSums(fsw.pred[,which.cell])
  place.points[, i] = dist.pred$PLACE.FSW[i]
}
colnames(dist.mat) = dist.pred$Name
colnames(place.points) = dist.pred$Name


dist.pred$PLACE = as.character(dist.pred$PLACE)
shapes = c(8, 17, 1)
names(shapes) = c("1", "2", "3")

dist.mat = dist.mat[,-ncol(dist.mat)] ## Remove Likoma
color_scheme_set(scheme = "darkgray")
dist_pred_plot = mcmc_areas(dist.mat[,order(colMeans(dist.mat), decreasing = TRUE)],
                            prob = 0.5, prob_outer = 0.95, point_est = "mean", area_method = "equal height") +
  xlim(0, 7100) +
  geom_point(mapping = aes(y = as.factor(Name), x = PLACE.FSW, shape = PLACE), data = dist.pred, size = 2) +
  theme_bw() + 
  scale_shape_manual(values = shapes, name = "Place Study", labels = c("1", "2", "None")) +
  coord_flip() + theme(axis.text.x = element_text(angle = 90)) +
  xlab("FSW Estimate") + ylab("District") +
  theme(axis.title = element_text(size = 15))

ggsave("fsw_dist_post.png", plot = dist_pred_plot, width = 9, height = 6, units = "in")



## The following are non-aggregated plots, not included in the paper
ggplot() + geom_tile(mapping = aes(x = dat.cell.all$x, y = dat.cell.all$y, col = cell.est$est, fill = cell.est$est), size = 0.01) +
  coord_equal() + 
  scale_fill_gradientn(colours = tim.colors()) +
  scale_color_gradientn(colours = tim.colors()) +
  labs(fill = "FSW", colour = "FSW", x = "Longitude", y = "Latitude") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15))

ggplot() + geom_tile(mapping = aes(x = dat.cell.all$x, y = dat.cell.all$y, col = log(cell.est$est), fill = log(cell.est$est)), size = 0.01) +
  coord_equal() + 
  scale_fill_gradientn(colours = tim.colors()) +
  scale_color_gradientn(colours = tim.colors()) +
  labs(fill = "Log-FSW", colour = "Log-FSW", x = "Longitude", y = "Latitude") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15))

ggplot() + geom_tile(mapping = aes(x = dat.cell.all$x, y = dat.cell.all$y, col = log(cell.est$var), fill = log(cell.est$var)), size = 0.01) +
  coord_equal() + 
  scale_fill_gradientn(colours = tim.colors()) +
  scale_color_gradientn(colours = tim.colors()) +
  labs(fill = "Log-Variance", colour = "Log-Variance", x = "Longitude", y = "Latitude") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15))









## Aggregate to larger cells for plotting

malawi.regs = raster::getData("GADM", country = "MWI", level = 1)
cell.est.sp = as.data.frame(cell.est)
names(cell.est.sp) = c("Estimate", "Variance")
cell.est.sp$x = dat.cell.all$x
cell.est.sp$y = dat.cell.all$y
coordinates(cell.est.sp) = ~x + y
crs(cell.est.sp) = crs(malawi.regs)

ptsreg = spsample(malawi.regs, 5000, type = "regular", offset = c(0.5, 0.5))
cells = SpatialPixels(ptsreg)

over.cells = over(cell.est.sp, cells)

new.vals = rep(NA, length(cells))
new.var = rep(NA, length(cells))
for(i in 1:length(cells)){
  new.vals[i] = sum(cell.est[which(over.cells == i), 1], na.rm = T)
  new.var[i] = sum(cell.est[which(over.cells == i), 2], na.rm = T)
}

ggplot() + geom_tile(mapping = aes(x = coordinates(cells)[,1], y = coordinates(cells)[,2],
                                   color = new.vals, fill = new.vals)) +
  coord_equal() + 
  scale_color_gradientn(colours = tim.colors()) +
  scale_fill_gradientn(colours = tim.colors(), guide = F) +
  labs(colour = "FSW", x = "Longitude", y = "Latitude") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15))

new.vals.zero = new.vals
new.vals.zero[new.vals.zero == 0] = NA
fsw_plot = ggplot() + geom_tile(mapping = aes(x = coordinates(cells)[,1], y = coordinates(cells)[,2],
                                              color = new.vals.zero, fill = new.vals.zero)) +
  coord_equal() + 
  scale_color_gradientn(colours = tim.colors()) +
  scale_fill_gradientn(colours = tim.colors(), guide = F) +
  labs(colour = "FSW", x = "Longitude", y = "Latitude") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15))

logfsw_plot = ggplot() + geom_tile(mapping = aes(x = coordinates(cells)[,1], y = coordinates(cells)[,2],
                                                 color = log(new.vals), fill = log(new.vals))) +
  coord_equal() + 
  scale_color_gradientn(colours = tim.colors()) +
  scale_fill_gradientn(colours = tim.colors(), guide = F) +
  labs(colour = "Log-FSW", x = "Longitude", y = "Latitude") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15))

logvar_plot = ggplot() + geom_tile(mapping = aes(x = coordinates(cells)[,1], y = coordinates(cells)[,2],
                                                 color = log(new.var), fill = log(new.var))) +
  coord_equal() + 
  scale_color_gradientn(colours = tim.colors()) +
  scale_fill_gradientn(colours = tim.colors(), guide = F) +
  labs(colour = "Log-Variance", x = "Longitude", y = "Latitude") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15))


ggsave("malawi_fsw_map.png", plot = fsw_plot, width = 6, height = 8, units = "in")
ggsave("malawi_fsw_map_log.png", plot = logfsw_plot, width = 6, height = 8, units = "in")
ggsave("malawi_fsw_map_var.png", plot = logvar_plot, width = 6, height = 8, units = "in")













# Do Three District Diagnostics -------------------------------------------

## (1) Get data for just 3 districts
## (2) Fit size and count model (full model, no offset)
## (3) Perform venue scaling
## (4) Check ZINB assumption
## (5) Check cell-level ven and cell-level and total FSW predictions


#################
## Step (1)
#################
three.dist = c("Mchinji","Mwanza","Neno")
dat.cell.three = subset(dat.cell.all, District %in% three.dist)
dat.venue.three = subset(dat.venue, District %in% three.dist)
dat.venue.three = subset(dat.venue.three, !is.na(fswestimate))
dat.venue.three = subset(dat.venue.three, fswestimate < 100)


#################
## Step (2)
#################
## No offset since fitting the full model
count.model.three = brm(bf(numplace ~ mv791 + v155 + built +
                             nightlight + (1|District),
                           zi ~ mv167 +
                             mv791 + v155 + built +
                             worldPop + nightlight + (1|District)),
                        data = dat.cell.three, chains = 1, cores = 1, iter = 1000,
                        family = "zero_inflated_negbinomial")


tmp = posterior_predict(count.model.three, newdata = dat.cell.three, nsamples = 5000, summary = F)
mean.vec = rowMeans(tmp)
pos.mean.vec = pos.sd.vec = rep(NA, nrow(tmp))
for(i in 1:nrow(tmp)){
  pos.ind = which(tmp[i,] > 0)
  pos.mean.vec[i] = mean(tmp[i,pos.ind])
  pos.sd.vec[i] = sd(tmp[i,pos.ind])
}
sd.vec = apply(tmp, 1, sd)
max.vec = apply(tmp, 1, max)
prob.zero.vec = apply(tmp, 1, function(x){sum(x == 0)})
overdisp.vec = apply(tmp, 1, var) / mean.vec
mean.obs = mean(dat.cell.three$numplace)
sd.obs = sd(dat.cell.three$numplace)
pos.mean.obs = mean(dat.cell.three$numplace[which(dat.cell.three$numplace > 0)])
pos.sd.obs = sd(dat.cell.three$numplace[which(dat.cell.three$numplace > 0)])
max.obs = max(dat.cell.three$numplace)
prob.zero.obs = sum(dat.cell.three$numplace == 0)
overdisp.obs = var(dat.cell.three$numplace) / mean.obs

## Mean
mean(mean.vec > mean.obs)
mean(pos.mean.vec > pos.mean.obs)
## SD
mean(sd.vec > sd.obs)
mean(pos.sd.vec > pos.sd.obs)
## Max
mean(max.vec > max.obs)
## prob.zero
mean(prob.zero.vec > prob.zero.obs)
## Overdisp
mean(overdisp.vec > overdisp.obs)



pp_check(count.model.three, nsamples = 5, type = "bars")

size.model.three = brm(bf(fswestimate ~ mv167 +
                            hivpos + worldPop + nightlight + (1|District),
                          hu ~ mv201 +
                            built + worldPop +
                            nightlight + (1|District)),
                       data = dat.venue.three, chains = 2, cores = 2, iter = 1000, family = "hurdle_lognormal")





tmp = posterior_predict(size.model.three, newdata = dat.venue.three, nsamples = 5000, summary = F)
mean.vec = rowMeans(tmp)
pos.mean.vec = pos.sd.vec = rep(NA, nrow(tmp))
for(i in 1:nrow(tmp)){
  pos.ind = which(tmp[i,] > 0)
  pos.mean.vec[i] = mean(tmp[i,pos.ind])
  pos.sd.vec[i] = sd(tmp[i,pos.ind])
}
sd.vec = apply(tmp, 1, sd)
max.vec = apply(tmp, 1, max)
prob.zero.vec = apply(tmp, 1, function(x){sum(x == 0)})
overdisp.vec = apply(tmp, 1, var) / mean.vec
mean.obs = mean(dat.venue.three$fswestimate)
sd.obs = sd(dat.venue.three$fswestimate)
pos.mean.obs = mean(dat.venue.three$fswestimate[which(dat.venue.three$fswestimate > 0)])
pos.sd.obs = sd(dat.venue.three$fswestimate[which(dat.venue.three$fswestimate > 0)])
max.obs = max(dat.venue.three$fswestimate)
prob.zero.obs = sum(dat.venue.three$fswestimate == 0)
overdisp.obs = var(dat.venue.three$fswestimate) / mean.obs

## Mean
mean(mean.vec > mean.obs)
mean(pos.mean.vec > pos.mean.obs)
## SD
mean(sd.vec > sd.obs)
mean(pos.sd.vec > pos.sd.obs)
## Max
mean(max.vec > max.obs)
## prob.zero
mean(prob.zero.vec > prob.zero.obs)
## Overdisp
mean(overdisp.vec > overdisp.obs)



#################
## Step (3)
#################


dist.pred.three = subset(dist.pred, Name %in% three.dist)

## Predict to all cells
count.model.pred.three = fitted(count.model.three, newdata = dat.cell.three, summary = F, nsamples = 5000)
size.model.pred.three = fitted(size.model.three, newdata = dat.cell.three, summary = F, nsamples = 5000)
fsw.pred.three = count.model.pred.three * size.model.pred.three
ven.pred.three = matrix(NA, nrow = nrow(fsw.pred.three), ncol = ncol(fsw.pred.three))

## Estimate 
for(i in 1:3){
  which.cells = which(dat.cell.three$District == dist.pred.three$Name[i])
  dist.count.pred = count.model.pred.three[,which.cells]
  dist.size.pred = size.model.pred.three[,which.cells]
  ## so that posterior mean is equal to expected number
  dist.venue.total.pred = sum(colMeans(dist.count.pred))
  if(!is.na(dist.pred.three$PLACE.Venue[i])){
    ## Scale each posterior sample
    scale.fact = dist.pred.three$PLACE.Venue[i] / rowSums(dist.count.pred)
    scale.fact = matrix(scale.fact, nrow = nrow(count.model.pred.three), ncol = ncol(dist.count.pred))
    
    
    dist.pred.three$Scaling.Factor[i] = mean(scale.fact)
    interval = quantile(rowSums(dist.count.pred))
    dist.pred.three$OG.Venue.est[i] = dist.venue.total.pred
    dist.pred.three$OG.Venue.lower[i] = interval[1]
    dist.pred.three$OG.Venue.upper[i] = interval[2]
    
    interval = quantile(rowSums(dist.count.pred * dist.size.pred), probs = c(0.05, 0.95))
    dist.pred.three$OG.FSW.est[i] = sum(colMeans(dist.count.pred * dist.size.pred))
    dist.pred.three$OG.FSW.lower[i] = interval[1]
    dist.pred.three$OG.FSW.upper[i] = interval[2]
    
    interval = quantile(rowSums(dist.count.pred * scale.fact * dist.size.pred), probs = c(0.05, 0.95))
    dist.pred.three$Scaled.FSW.est[i] = sum(colMeans(dist.count.pred * scale.fact * dist.size.pred))
    dist.pred.three$Scaled.FSW.lower[i] = interval[1]
    dist.pred.three$Scaled.FSW.upper[i] = interval[2]
    
    ven.pred.three[,which.cells] = dist.count.pred * scale.fact
    fsw.pred.three[,which.cells] = fsw.pred.three[,which.cells] * scale.fact
  }else{
    interval = quantile(rowSums(dist.count.pred))
    dist.pred.three$OG.Venue.est[i] = dist.venue.total.pred
    dist.pred.three$OG.Venue.lower[i] = interval[1]
    dist.pred.three$OG.Venue.upper[i] = interval[2]
    
    interval = quantile(rowSums(dist.count.pred * dist.size.pred), probs = c(0.05, 0.95))
    dist.pred.three$OG.FSW.est[i] = sum(colMeans(dist.count.pred * dist.size.pred))
    dist.pred.three$OG.FSW.lower[i] = interval[1]
    dist.pred.three$OG.FSW.upper[i] = interval[2]
  }
  print(c(i,"/",3))
}
dist.pred.three




#################
## Step (4)
#################
pp_check(count.model.three, type = "rootogram", style = "hanging", transform = "round", nsamples = 4000) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14)) +
  xlim(-1, 100)
pp_check(count.model.three, type = "scatter_avg", nsamples = 4000) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))
pp_check(count.model.three, type = "stat_2d", nsamples = 4000) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))





# Look at conditional distribution ----------------------------------------

which.pos = which(dat.cell.three$numplace > 0)
cell.three.pos = subset(dat.cell.three, numplace > 0)
fsw.pred.three.ave = predict(size.model.three, newdata = cell.three.pos, nsamples = 5000)
ven.resid = cell.three.pos$numplace - colMeans(ven.pred.three)[which.pos]
fsw.resid = cell.three.pos$fswave - fsw.pred.three.ave[,1]


three_resid = ggplot() + geom_point(mapping = aes(x = ven.resid, y = fsw.resid)) +
  xlab("Venue Count Residual") + ylab("Venue Size Residual") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))


ggsave("venue_resid_corr.png", plot = three_resid, width = 8, height = 6, units = "in")

