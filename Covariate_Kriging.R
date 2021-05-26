
library(sp)
library(raster)
library(ggplot2)
library(ggmap) ## register_google
library(foreign) ## Read dbf
library(fields)
library(geosphere) ## For distm



## Read in the WorldPop data
worldPop.raw = raster("MWI_ppp_2015_adj_v2.tif")
worldPop = cbind(coordinates(worldPop.raw), getValues(worldPop.raw))
worldPop = as.data.frame(worldPop)
names(worldPop) = c("x", "y", "wp")
## Want to keep the NA's for averaging
## NA's represent inhabitable zones, which means zero population for the purpose of this analysis
worldPop$wp[is.na(worldPop$wp)] = 0
worldPop.sp = worldPop
coordinates(worldPop.sp) = ~ x + y


## Look at other covariates

dat.place.1 = readRDS("dat_place1.rds")
dat.place = place.dat = readRDS("dat_place2.rds")
dat.place = dat.place[!is.na(dat.place$Latitude),] ## Can't use venues without Latitude and Longitude
dat.place$b5a.y = factor(dat.place$b5a.y)
dat.place$fswestimate[dat.place$fswestimate > 100] = 100
dat.place.1$fswestimate[dat.place$fswestimate > 100] = 100
dhs.cov.tmp = read.csv("MWGC7AFL.csv")
dhs.dat.tmp = read.dbf("MWGE7AFL.dbf")
dat.dhs = cbind(dhs.cov.tmp, dhs.dat.tmp)


nightlight.dat = read.csv("malawi2016.csv")
nightlight.dat = nightlight.dat[!is.na(nightlight.dat$FID_),]
nightlight.ave = rowMeans(nightlight.dat[,-c(1,2,3)], na.rm = T)


register_google('AIzaSyB3w8SKisb8Fe5cvkG2cO6snR7L-EmbqlU')
malawi.map = get_map(location = "Malawi", zoom=6, maptype= "hybrid")


nightlight.sp = data.frame(x = nightlight.dat$Xlongitude,
                           y = nightlight.dat$Ylatitude,
                           nightlight = nightlight.ave)
coordinates(nightlight.sp) = ~x + y

malawi = raster::getData("GADM", country="MWI", level=0)
malawi.2 = raster::getData("GADM", country = "MWI", level = 1)
plot(malawi.2, xlim = c(33, 35.01), ylim = c(-14.5, -13.5))

ptsreg = spsample(malawi, 50000, type = "regular", offset = c(0.5, 0.5))
cells = SpatialPixels(ptsreg)

crs(nightlight.sp) = crs(malawi)
crs(worldPop.sp) = crs(malawi)


# Get just cells with positive WorldPop -----------------------------------
cells.all = cells
worldpop.over = over(worldPop.sp, cells, fn = mean, na.rm = T)
worldPop.sp@data$cell = worldpop.over
worldpop.ave = (tapply(worldPop.sp@data$wp, worldPop.sp@data$cell, FUN=mean, na.rm = T))
which.in = as.numeric(names(worldpop.ave))

cells.out = cells[-which.in,]
cells = cells[which.in,]

ggplot(as.data.frame(cells), aes(x = x1, y = x2, fill = worldpop.ave)) +
  geom_tile() +
  coord_equal()


tmp = cbind(coordinates(cells), worldpop.ave)
tmp = as.data.frame(subset(tmp, worldpop.ave == 0))

ggmap(malawi.map, extent = "panel") +
  scale_x_continuous( limits = c( 32.7 , 36 ) , expand = c( 0 , 0 ) ) +
  scale_y_continuous( limits = c( -17.2 , -9.3 ) , expand = c( 0 , 0 ) ) +
  geom_point(data = tmp, mapping = aes(x = x1, y = x2, fill = worldpop.ave),
             shape = 3, color = "red") +
  coord_equal()


# Now get nightlight and district name ------------------------------------

nightlight.over = over(nightlight.sp, cells)
nightlight.sp@data$cell = nightlight.over
nightlight.ave = tapply(nightlight.sp@data$nightlight, nightlight.sp@data$cell, FUN = mean)

ggplot(as.data.frame(cells), aes(x = x1, y = x2, fill = nightlight.ave)) +
  geom_tile() +
  coord_equal()


dist.name = over(cells, malawi.2)$NAME_1
dist.name.old = dist.name

## Also get the number of venues in each cell
place.sp = dat.place
coordinates(place.sp) = ~ Longitude + Latitude
crs(place.sp) = crs(cells)
ven.over = over(place.sp, cells)
numplace = table(factor(ven.over, levels = 1:length(cells)))
numplace = as.data.frame(numplace)$Freq
fsw.ave = tapply(place.sp@data$fswestimate, ven.over, FUN = mean, na.rm = T)




### Assign cells over borders to the district with the most venues (fix dist.name)
ven.over.unique = unique(ven.over)
ven.over.unique = ven.over.unique[complete.cases(ven.over.unique)]
for(i in ven.over.unique){
  tmp = over(place.sp, cells[i])
  tmp = which(!is.na(tmp))
  dist.tab = as.data.frame(table(dat.place$b5a.y[tmp]))
  most.dist = as.character(dist.tab[order(dist.tab$Freq, decreasing = T),"Var1"][1])
  dist.name[i] = most.dist
}

dist.name[dist.name == "Chikhwawa"] = "Chikwawa"
cbind(dist.name.old, dist.name)




cells.new = cells
cells.new$wp = worldpop.ave
cells.new$nightlight = nightlight.ave
cells.new$district = dist.name
cells.new$numplace = numplace
cells.new$fswave = NA
cells.new$fswave[as.numeric(names(fsw.ave))] = fsw.ave


regress.ind = which(numplace > 0)
regress.cells = cells.new[regress.ind,]
cell.fsw.tmp = cells.new$fswave
cell.fsw = rep(NA, length(cells.new$fswave))
cell.fsw[regress.ind] = cell.fsw.tmp[regress.ind]

plot(malawi)
plot(regress.cells)
plot(cells[regress.cells,], add = T)
points(dat.place$Longitude, dat.place$Latitude, pch = 20, cex = 0.02, col="red")



## Extract info for Mwanza, Mchinji, and Neno
cells.new$fswsum = rep(NA, nrow(cells.new))
for(i in 1:length(cells)){
  which.i = which(ven.over == i) ## This pulls the venue ids that are in cell i
  cells.new$fswsum[i] = sum(place.sp$fswestimate[which.i], na.rm = T)
}
cells.new.three = subset(cells.new, district %in% c("Mchinji", "Mwanza", "Neno"))

saveRDS(object = cells.new.three, file = "cells_new_three.rds")







# Now do kriging ----------------------------------------------------------




## Convert all relevent covariates
dhs.birth = read.csv("birth_agg.csv")
dhs.child = read.csv("child_agg.csv")
dhs.couple = read.csv("couple_agg.csv")
dhs.hiv = read.csv("hiv_agg.csv")
dhs.household = read.csv("household_agg.csv")
dhs.housemember = read.csv("housemember_agg.csv")
dhs.individual = read.csv("individual_agg.csv")
dhs.men = read.csv("men_agg.csv")
dhs.perc = read.csv("percent_agg.csv")
dhs.geospat = read.csv("MWGC7AFL.csv")
geospat.rem = which(dhs.geospat$All_Population_Count_2005 < -9000) ## Removes 10 columns with bad data


dhs.men = dhs.men[,c("mv012","mv035", "mv167","mv191","mv201","mv245","mv531", "mv766a")]
dhs.household = dhs.household[,c("hv204")]
dhs.individual = dhs.individual[,c("v104","v115","v167")]


mv168.na = is.na(dhs.perc$mv168)
mv791.na = is.na(dhs.perc$mv793)
mv793.na = is.na(dhs.perc$mv793)
v168.na = is.na(dhs.perc$v168)


tmp.all = data.frame("x" = dat.dhs$LONGNUM, "y"= dat.dhs$LATNUM, 
                     dhs.men,
                     dhs.household,
                     dhs.perc)
coordinates(tmp.all) = ~ x + y
proj4string(tmp.all) = proj4string(cells)

dhs.geospat.cov = data.frame("x" = dat.dhs$LONGNUM[-geospat.rem], "y" = dat.dhs$LATNUM[-geospat.rem],
                             dhs.geospat[-geospat.rem,])
coordinates(dhs.geospat.cov) = ~x + y
proj4string(dhs.geospat.cov) = proj4string(cells)


## Plot some DHS clusters
dat.dhs.tmp = dat.dhs[,c(87,88)]
dat.dhs.tmp$rural = as.factor(dat.dhs$URBAN_RURA)
dat.dhs.tmp$wealth = dhs.men$mv191

malawi.terrain = get_map(location = "Malawi", zoom=6, maptype= "terrain")
rural_map = ggmap(malawi.terrain, extent = "panel") +
  scale_x_continuous( limits = c( 32.7 , 36 ) , expand = c( 0 , 0 ) ) +
  scale_y_continuous( limits = c( -17.2 , -9.3 ) , expand = c( 0 , 0 ) ) +
  geom_point(mapping = aes(x = LONGNUM,
                             y = LATNUM,
                             col = rural), data = dat.dhs.tmp, show.legend = T) +
  labs(x = "Longitude", y = "Latitude", col = "Urban-Rural") +
  scale_color_discrete(labels = c("Rural", "Urban")) + 
  theme(text = element_text(size = 15), legend.title = element_text(size = 15))

wealth_map = ggmap(malawi.terrain, extent = "panel") +
  scale_x_continuous( limits = c( 32.7 , 36 ) , expand = c( 0 , 0 ) ) +
  scale_y_continuous( limits = c( -17.2 , -9.3 ) , expand = c( 0 , 0 ) ) +
  geom_point(mapping = aes(x = LONGNUM,
                           y = LATNUM,
                           col = wealth), data = dat.dhs.tmp, show.legend = T) +
  labs(x = "Longitude", y = "Latitude", col = "Wealth index\nfactor score") +
  theme(text = element_text(size = 15), legend.title = element_text(size = 15)) +
  scale_colour_gradientn(colours = tim.colors())

ggsave("rural_map.png", plot = rural_map, width = 5, height = 7, units = "in")
ggsave("wealth_map.png", plot = wealth_map, width = 5, height = 7, units = "in")



# Kriging for cells and venues -------------------------------------------------
## This uses the fields package, and the spatialProcess function
dat.place.1$place = 1
dat.place$place = 2
dat.place.both = rbind(dat.place.1, dat.place)
venue.sp = dat.place.both
coordinates(venue.sp) = ~ Longitude + Latitude
proj4string(venue.sp) = proj4string(tmp.all)
bin = tmp.all
poly.coords = coordinates(cells.new)

fit.mv012 = spatialProcess(x = coordinates(bin), y = bin$mv012)
pred.cell.mv012 = predict.mKrig(fit.mv012, xnew = coordinates(cells.new))
se.cell.mv012 = predictSE.mKrig(fit.mv012, xnew = coordinates(cells.new))
pred.ven.mv012 = predict.mKrig(fit.mv012, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.mv012 = predict.mKrig(fit.mv012, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))


## mv035
fit.mv035 = spatialProcess(x = coordinates(bin), y = bin$mv035)
pred.cell.mv035 = predict.mKrig(fit.mv035, xnew = coordinates(cells.new))
se.cell.mv035 = predictSE.mKrig(fit.mv035, xnew = coordinates(cells.new))
pred.ven.mv035 = predict.mKrig(fit.mv035, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.mv035 = predictSE.mKrig(fit.mv035, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))


## mv167
fit.mv167 = spatialProcess(x = coordinates(bin), y = bin$mv167)
pred.cell.mv167 = predict.mKrig(fit.mv167, xnew = coordinates(cells.new))
se.cell.mv167 = predictSE.mKrig(fit.mv167, xnew = coordinates(cells.new))
pred.ven.mv167 = predict.mKrig(fit.mv167, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.mv167 = predictSE.mKrig(fit.mv167, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))


## mv191
fit.mv191 = spatialProcess(x = coordinates(bin), y = bin$mv191)
pred.cell.mv191 = predict.mKrig(fit.mv191, xnew = coordinates(cells.new))
se.cell.mv191 = predictSE.mKrig(fit.mv191, xnew = coordinates(cells.new))
pred.ven.mv191 = predict.mKrig(fit.mv191, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.mv191 = predictSE.mKrig(fit.mv191, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))

## mv201
fit.mv201 = spatialProcess(x = coordinates(bin), y = bin$mv201)
pred.cell.mv201 = predict.mKrig(fit.mv201, xnew = coordinates(cells.new))
se.cell.mv201 = predictSE.mKrig(fit.mv201, xnew = coordinates(cells.new))
pred.ven.mv201 = predict.mKrig(fit.mv201, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.mv201 = predictSE.mKrig(fit.mv201, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))



## mv245
fit.mv245 = spatialProcess(x = coordinates(bin), y = bin$mv245)
pred.cell.mv245 = predict.mKrig(fit.mv245, xnew = coordinates(cells.new))
se.cell.mv245 = predictSE.mKrig(fit.mv245, xnew = coordinates(cells.new))
pred.ven.mv245 = predict.mKrig(fit.mv245, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.mv245 = predictSE.mKrig(fit.mv245, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))


## mv531
fit.mv531 = spatialProcess(x = coordinates(bin), y = bin$mv531)
pred.cell.mv531 = predict.mKrig(fit.mv531, xnew = coordinates(cells.new))
se.cell.mv531 = predictSE.mKrig(fit.mv531, xnew = coordinates(cells.new))
pred.ven.mv531 = predict.mKrig(fit.mv531, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.mv531 = predictSE.mKrig(fit.mv531, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))



## mv766a
fit.mv766a = spatialProcess(x = coordinates(bin), y = bin$mv766a)
pred.cell.mv766a = predict.mKrig(fit.mv766a, xnew = coordinates(cells.new))
se.cell.mv766a = predictSE.mKrig(fit.mv766a, xnew = coordinates(cells.new))
pred.ven.mv766a = predict.mKrig(fit.mv766a, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.mv766a = predictSE.mKrig(fit.mv766a, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))


## dhs.household
fit.dhs.household = spatialProcess(x = coordinates(bin), y = bin$dhs.household)
pred.cell.dhs.household = predict.mKrig(fit.dhs.household, xnew = coordinates(cells.new))
se.cell.dhs.household = predictSE.mKrig(fit.dhs.household, xnew = coordinates(cells.new))
pred.ven.dhs.household = predict.mKrig(fit.dhs.household, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.dhs.household = predictSE.mKrig(fit.dhs.household, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))

## mv168
fit.mv168 = spatialProcess(x = coordinates(bin), y = bin$mv168)
pred.cell.mv168 = predict.mKrig(fit.mv168, xnew = coordinates(cells.new))
se.cell.mv168 = predictSE.mKrig(fit.mv168, xnew = coordinates(cells.new))
pred.ven.mv168 = predict.mKrig(fit.mv168, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.mv168 = predictSE.mKrig(fit.mv168, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))


## mv791
fit.mv791 = spatialProcess(x = coordinates(bin), y = bin$mv791)
pred.cell.mv791 = predict.mKrig(fit.mv791, xnew = coordinates(cells.new))
se.cell.mv791 = predictSE.mKrig(fit.mv791, xnew = coordinates(cells.new))
pred.ven.mv791 = predict.mKrig(fit.mv791, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.mv791 = predictSE.mKrig(fit.mv791, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))

## mv793
fit.mv793 = spatialProcess(x = coordinates(bin), y = bin$mv793)
pred.cell.mv793 = predict.mKrig(fit.mv793, xnew = coordinates(cells.new))
se.cell.mv793 = predictSE.mKrig(fit.mv793, xnew = coordinates(cells.new))
pred.ven.mv793 = predict.mKrig(fit.mv793, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.mv793 = predictSE.mKrig(fit.mv793, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))


## v155
fit.v155 = spatialProcess(x = coordinates(bin), y = bin$v155)
pred.cell.v155 = predict.mKrig(fit.v155, xnew = coordinates(cells.new))
se.cell.v155 = predictSE.mKrig(fit.v155, xnew = coordinates(cells.new))
pred.ven.v155 = predict.mKrig(fit.v155, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.v155 = predictSE.mKrig(fit.v155, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))


## v168
fit.v168 = spatialProcess(x = coordinates(bin), y = bin$v168)
pred.cell.v168 = predict.mKrig(fit.v168, xnew = coordinates(cells.new))
se.cell.v168 = predictSE.mKrig(fit.v168, xnew = coordinates(cells.new))
pred.ven.v168 = predict.mKrig(fit.v168, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.v168 = predictSE.mKrig(fit.v168, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))

## hiv.test
fit.hiv.test = spatialProcess(x = coordinates(bin), y = bin$hiv.test)
pred.cell.hiv.test = predict.mKrig(fit.hiv.test, xnew = coordinates(cells.new))
se.cell.hiv.test = predictSE.mKrig(fit.hiv.test, xnew = coordinates(cells.new))
pred.ven.hiv.test = predict.mKrig(fit.hiv.test, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.hiv.test = predictSE.mKrig(fit.hiv.test, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))

## Also need to krige geospatial covariates


## Work with the geospatial covariates
## Keep Aridity, BUILT_Population_2014, Proximity_to_National_Borders, Proximity_to_Water, Travel_Times

## Aridity
dhs.cov.tmp$Aridity[dhs.cov.tmp$Aridity < 0] = NA
fit.Aridity = spatialProcess(x = coordinates(bin), y = dhs.cov.tmp$Aridity)
pred.cell.Aridity = predict.mKrig(fit.Aridity, xnew = coordinates(cells.new))
se.cell.Aridity = predictSE.mKrig(fit.Aridity, xnew = coordinates(cells.new))
pred.ven.Aridity = predict.mKrig(fit.Aridity, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.Aridity = predictSE.mKrig(fit.Aridity, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))


## built
dhs.cov.tmp$BUILT_Population_2014[dhs.cov.tmp$BUILT_Population_2014 < 0] = NA
fit.built = spatialProcess(x = coordinates(bin), y = dhs.cov.tmp$BUILT_Population_2014)
pred.cell.built = predict.mKrig(fit.built, xnew = coordinates(cells.new))
se.cell.built = predictSE.mKrig(fit.built, xnew = coordinates(cells.new))
pred.ven.built = predict.mKrig(fit.built, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.built = predictSE.mKrig(fit.built, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))

## borders
fit.borders = spatialProcess(x = coordinates(bin), y = dhs.cov.tmp$Proximity_to_National_Borders)
pred.cell.borders = predict.mKrig(fit.borders, xnew = coordinates(cells.new))
se.cell.borders = predictSE.mKrig(fit.borders, xnew = coordinates(cells.new))
pred.ven.borders = predict.mKrig(fit.borders, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.borders = predictSE.mKrig(fit.borders, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))

## water
fit.water = spatialProcess(x = coordinates(bin), y = dhs.cov.tmp$Proximity_to_Water)
pred.cell.water = predict.mKrig(fit.water, xnew = coordinates(cells.new))
se.cell.water = predictSE.mKrig(fit.water, xnew = coordinates(cells.new))
pred.ven.water = predict.mKrig(fit.water, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.water = predictSE.mKrig(fit.water, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))

## travel
dhs.cov.tmp$Travel_Times[dhs.cov.tmp$Travel_Times < 0] = NA
fit.travel = spatialProcess(x = coordinates(bin), y = dhs.cov.tmp$Travel_Times)
pred.cell.travel = predict.mKrig(fit.travel, xnew = coordinates(cells.new))
se.cell.travel = predictSE.mKrig(fit.travel, xnew = coordinates(cells.new))
pred.ven.travel = predict.mKrig(fit.travel, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))
se.ven.travel = predictSE.mKrig(fit.travel, xnew = cbind(dat.place.both$Longitude, dat.place.both$Latitude))



input.venue = dat.place.both
input.venue$mv012 = pred.ven.mv012
input.venue$mv035 = pred.ven.mv035
input.venue$mv167 = pred.ven.mv167
input.venue$mv191 = pred.ven.mv191
input.venue$mv201 = pred.ven.mv201
input.venue$mv245 = pred.ven.mv245
input.venue$mv531 = pred.ven.mv531
input.venue$mv766a = pred.ven.mv766a
input.venue$hv204 = pred.ven.dhs.household
input.venue$mv168 = pred.ven.mv168
input.venue$mv791 = pred.ven.mv791
input.venue$mv793 = pred.ven.mv793
input.venue$v155 = pred.ven.v155
input.venue$v168 = pred.ven.v168
input.venue$hivpos = pred.ven.hiv.test
input.venue$aridity = pred.ven.Aridity
input.venue$built = pred.ven.built
input.venue$borders = pred.ven.borders
input.venue$water = pred.ven.water
input.venue$travel = pred.ven.travel

input.venue$mv012.se = se.ven.mv012
input.venue$mv035.se = se.ven.mv035
input.venue$mv167.se = se.ven.mv167
input.venue$mv191.se = se.ven.mv191
input.venue$mv201.se = se.ven.mv201
input.venue$mv245.se = se.ven.mv245
input.venue$mv531.se = se.ven.mv531
input.venue$mv766a.se = se.ven.mv766a
input.venue$hv204.se = se.ven.dhs.household
input.venue$mv168.se = se.ven.mv168
input.venue$mv791.se = se.ven.mv791
input.venue$mv793.se = se.ven.mv793
input.venue$v155.se = se.ven.v155
input.venue$v168.se = se.ven.v168
input.venue$hivpos.se = se.ven.hiv.test
input.venue$aridity.se = se.ven.Aridity
input.venue$built.se = se.ven.built
input.venue$borders.se = se.ven.borders
input.venue$water.se = se.ven.water
input.venue$travel.se = se.ven.travel

## Get closest worldPop and nightlight for the venues

venue.poly.ind = over(venue.sp, cells.new)
input.venue$newdistrict = venue.poly.ind$district

pos.cells = cells.new[unique(which(cells.new$nightlight %in% venue.poly.ind$nightlight)),]
nightlight.pos.ind = over(nightlight.sp, pos.cells)
nightlight.pos.ind = which(!is.na(nightlight.pos.ind$wp))
nightlight.small = nightlight.sp[nightlight.pos.ind,]
nightlight.venue.dist = distm(coordinates(venue.sp), coordinates(nightlight.small))
min.dist = apply(nightlight.venue.dist, 1, function(x) order(x, decreasing=F)[2])
input.venue$nightlight = nightlight.small$nightlight[min.dist]

ggplot() + geom_point(mapping = aes(x = coordinates(nightlight.small)[,1], y = coordinates(nightlight.small)[,2])) +
  coord_equal()

ggplot() + geom_point(mapping = aes(x = coordinates(venue.sp)[,1], y = coordinates(venue.sp)[,2])) +
  coord_equal()


pos.cells = cells.new[unique(which(cells.new$wp %in% venue.poly.ind$wp)),]
worldPop.pos.ind = over(worldPop.sp, pos.cells)
worldPop.pos.ind = which(!is.na(worldPop.pos.ind$wp))
worldPop.small = worldPop.sp[worldPop.pos.ind,]
worldPop.small = subset(worldPop.small, wp > 0)
worldPop.venue.dist = distm(coordinates(venue.sp), coordinates(worldPop.small))
min.dist = apply(worldPop.venue.dist, 1, function(x) order(x, decreasing=F)[2])
input.venue$worldPop = worldPop.small$wp[min.dist]

ggplot() + geom_point(mapping = aes(x = coordinates(worldPop.small)[,1], y = coordinates(worldPop.small)[,2])) +
  coord_equal()

ggplot() + geom_point(mapping = aes(x = coordinates(venue.sp)[,1], y = coordinates(venue.sp)[,2])) +
  coord_equal()





## Also get the number of places
poly.sp = data.frame(poly.coords)
coordinates(poly.sp) = ~ x1 + x2
crs(poly.sp) = crs(malawi.2)


all.cell.cov = data.frame(cells.new$fswave, cells.new$district, cells.new$numplace, pred.cell.mv012, pred.cell.mv035, pred.cell.mv167, pred.cell.mv191,
                          pred.cell.mv201, pred.cell.mv245, pred.cell.mv531, pred.cell.mv766a, pred.cell.dhs.household,
                          pred.cell.mv168, pred.cell.mv791, pred.cell.mv793, pred.cell.v155, pred.cell.v168,
                          pred.cell.hiv.test, pred.cell.Aridity, pred.cell.built, pred.cell.borders,
                          pred.cell.water, pred.cell.travel)

all.cell.var = data.frame(se.cell.mv012, se.cell.mv035, se.cell.mv167, se.cell.mv191, se.cell.mv201,
                          se.cell.mv245,
                          se.cell.mv531, se.cell.mv766a, se.cell.dhs.household, se.cell.mv168,
                          se.cell.mv791, se.cell.mv793, se.cell.v155, se.cell.v168, se.cell.hiv.test,
                          se.cell.Aridity, se.cell.built, se.cell.borders, se.cell.water, se.cell.travel)


names(all.cell.cov) = c("fswave", "District", "numplace", "mv012","mv035","mv167","mv191","mv201", "mv245","mv531","mv766a","hv204","mv168","mv791",
                        "mv793","v155","v168","hivpos","aridity","built","borders","water","travel")
names(all.cell.var) = c("mv012","mv035","mv167","mv191","mv201", "mv245","mv531","mv766a","hv204","mv168","mv791",
                        "mv793","v155","v168","hivpos","aridity","built","borders","water","travel")


output = data.frame(all.cell.cov, worldPop = cells.new@data$wp, nightlight = cells.new@data$nightlight)
output.cells = output
output.cells$x = poly.coords[,1]
output.cells$y = poly.coords[,2]

all.cell.var$x = poly.coords[,1]
all.cell.var$y = poly.coords[,2]

## Save all cells
saveRDS(object = output.cells, file = "Malawi_Cells_worldPop.rds")
## and the variances
saveRDS(object = all.cell.var, file = "Malawi_Cells_worldPop_se.rds")

## Save all venues
saveRDS(object = input.venue, file = "Malawi_Venues_worldPop.rds")


## Save average fsw per cell
saveRDS(object = cells.new, file = "Malawi_Cells_fsw.rds")




















