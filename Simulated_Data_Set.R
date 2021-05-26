######################################################################
##
## This file is used to create a pseudo-dataset to run Malawi_Fit_With_Sim.R
##
######################################################################
dat.cell = readRDS("Malawi_Cells_worldPop_sim.rds")
cell.coords = readRDS("Malawi_coords.rds")
venue.names = readRDS("Venue_column_names.rds")

set.seed(737)

districts = c("Balaka", "Blantyre", "Chikwawa", "Chiradzulu", "Chitipa", "Dedza",
              "Dowa", "Karonga", "Kasungu", "Likoma", "Lilongwe", "Machinga",
              "Mangochi", "Mchinji", "Mulanje", "Mwanza", "Mzimba", "Neno", 
              "Nkhata Bay", "Nkhotakota", "Nsanje", "Ntcheu", "Ntchisi", "Phalombe",
              "Rumphi", "Salima", "Thyolo", "Zomba"
)

venue.dist = c("Balaka", "Blantyre", "Chikwawa", "Dedza", "Dowa", "Karonga",
               "Kasungu", "Lilongwe", "Machinga", "Mangochi", "Mchinji", "Mwanza",
               "Mzimba", "Mzuzu", "Neno", "Nkhata Bay", "Nkhotakota", "Ntcheu",
               "Rumphi", "Salima", "Zomba")

pos.dist = c("Balaka", "Chikwawa", "Dedza", "Dowa", "Karonga", "Kasungu", "Mchinji",
             "Mwanza", "Mzimba", "Neno", "Nkhata Bay", "Nkhotakota", "Ntcheu", 
             "Rumphi", "Salima")

## First create dat.venue
dat.venue = matrix(rnorm(2541 * 154), nrow = 2541, ncol = 154)
dat.venue = as.data.frame(dat.venue)
names(dat.venue) = venue.names
dat.venue$b5a.y = sample(venue.dist, size = 2541, replace = T)
dat.venue$fswestimate = rbinom(2541, 1, prob = plogis(dat.venue$mv201 + dat.venue$built +
                                  dat.venue$worldPop +
                                    dat.venue$nightlight)) *
                        round(rlnorm(2541, meanlog = exp(-1 + 0.1 * dat.venue$mv167 + 0.1 * dat.venue$hivpos +
                                  0.1 * dat.venue$worldPop + 0.1 * dat.venue$nightlight)))
venue.loc.sample = sample(1:nrow(cell.coords), nrow(dat.venue))
dat.venue$Longitude = cell.coords$x[venue.loc.sample]
dat.venue$Latitude = cell.coords$y[venue.loc.sample]


dat.cell = matrix(runif(39312 * 27), nrow = 39312)
dat.cell = as.data.frame(dat.cell)
names(dat.cell) = c("fswave", "District", "numplace", "mv012", "mv035", "mv167",
                    "mv191", "mv201", "mv245", "mv531",
                    "mv766a", "hv204", "mv168", "mv791", "mv793", "v155", "v168",
                    "hivpos", "aridity", "built", "borders", "water", "travel",
                    "worldPop", "nightlight", "x", "y")
dat.cell$District = sample(districts, size = 39312, replace = T)
dat.cell$worldPop[4:5] = 0
dat.cell$nightlight[8:9] = 0

dat.cell$numplace = rbinom(39312, 1, prob = plogis(-6.5 + dat.cell$mv167 + 
                 dat.cell$v155 + dat.cell$built + dat.cell$worldPop + dat.cell$nightlight)) *
  round(rnbinom(39312, mu = exp(-1 + 0.1 * dat.cell$v155 + 1 * dat.cell$built +
                 0.5 * dat.cell$worldPop), size = 1))

dat.cell$numplace[!(dat.cell$District %in% pos.dist)] = 0

dat.cell$x = cell.coords$x
dat.cell$y = cell.coords$y


saveRDS(dat.venue, file = "Malawi_Venues_nightlight_sim.rds")
saveRDS(dat.cell, file = "Malawi_Cells_worldPop_sim.rds")

