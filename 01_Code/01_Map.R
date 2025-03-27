#Author: Marion Chevrinais
#Date: July 2024

# Load required packages ####
library(terra)
library(geodata)
library(CoordinateCleaner)  # species occurrence data cleaning
library(sdmpredictors)  # environmental (including marine) data download
library(fuzzySim)  # extract presence-absence and environmental values for modelling
library(raster)  # needed for converting SpatVector to Spatial (e.g. for 'cartogram')
library(cartogram)  # cartograms
library(mapview) # interactive maps
library(car)
library(flextable)
library(magrittr)
library(dplyr)
library(ggplot2)

# load csv file with gps coordinates of DNA samples
eDNA_phoca <- read.csv2(file.path(here::here(),"00_Data/01_Map/GPS coordinate_phoca.csv"))
str(eDNA_phoca) #verify that latitude and longitude are numerical

plot(eDNA_phoca[ , c("long", "lat")])
eDNA_phoca_points <- vect(eDNA_phoca, geom = c("long", "lat"), keepgeom = TRUE, crs = "+proj=longlat") #vectorize latitude and longitude

# load vector files for waterbodies 
countries <- vect(file.path(here::here(),"00_Data/01_Map/countries.gpkg"))
hydro <- vect(file.path(here::here(),"00_Data/01_Map/waterbody_2.shp"))
riv <- hydro[perim(hydro) > 15000,]
crit<-vect(file.path(here::here(),"00_Data/01_Map/Habitat potentiel.shp"))
habitat.essentiel <- vect(file.path(here::here(),"00_Data/01_Map/HabitatEssentiel_PhoqueCommunLLM.shp"))
BV<-vect(file.path(here::here(),"00_Data/01_Map/BV_N1M_S.shp"))

# Plot characteristics
e<- ext(8e+05,8.7e+05,6770000,6850000)
Ca <- crop(countries,ext(c(-85,-50,43,70)))

# the right projection for the area
Ca2<-project(Ca, "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45")
crs(Ca2)
plot(Ca2)

#apply the right projection and extent to all layers
Countries2<-project(countries, "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45")

Nunavik3 <- crop(Countries2,ext(c(8e+05,1300000,6450000,6850000)))
Nunavik3 <- project(Nunavik3, "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45")

crit2<-project(crit, "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45")
crs(crit2)

riv2<-project(riv, "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45")
riv3 <- crop(riv2,ext(c(8e+05,1300000,6450000,6850000)))
crs(riv3)

BV2<-project(BV, "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45") 
BV3<- crop(BV2,ext(c(8e+05,1300000,6450000,6850000)))

habitat.essentiel2 <-project(habitat.essentiel, "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45") 
habitat.essentiel3<- crop(habitat.essentiel2,ext(c(8e+05,1300000,6450000,6850000)))

eDNA_points2 <-project(eDNA_phoca_points, "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45")

# subset by sampling year
A2019 <- eDNA_points2[eDNA_points2$Year=='2019',]
A2019_det<- A2019[A2019$max_copies_detected !=0,]
A2020<- eDNA_points2[eDNA_points2$Year=='2020',]
A2020_det<- A2020[A2020$max_copies_detected !=0,]
A2021<- eDNA_points2[eDNA_points2$Year=='2021',]
A2021_det<- A2021[A2021$max_copies_detected !=0,]
A2022<- eDNA_points2[eDNA_points2$Year=='2022',]
A2022_det<- A2022[A2022$max_copies_detected !=0,]
A2023<- eDNA_points2[eDNA_points2$Year=='2023',]
A2023_det<- A2023[A2023$max_copies_detected !=0,]

# plot - figure_Phoca
plot(Nunavik3, col='grey95', pax = list(side=1:2), border='grey95' )
#plot(countries2, ext=window, col='grey10', add=T)
plot(riv3, add=TRUE, col='lightskyblue2',border='lightskyblue2', lty = 1
) #add rivers
plot(BV3,add=TRUE, border='dodgerblue4') #add watersheds
plot(crit2,add=T, col='grey60',border='grey60') #add potential habitat of LLMHS
plot(habitat.essentiel2, add=T, col='grey30', border="grey30") #add critical habitat of LLMHS
#select the year of interest below to add sampling and detection points
plot(A2019, add=TRUE, col='black', pch=21, cex=1, bg='white')
plot(A2019_det, add=TRUE, col='black', pch=21, cex=1, bg='black')
plot(A2020, add=TRUE, col='black', pch=21, cex=1, bg='white')
plot(A2021, add=TRUE, col='black', pch=21, cex=1, bg='white')
plot(A2021_det, add=TRUE, col='black', pch=21, cex=1, bg='black')
plot(A2022, add=TRUE, col='black', pch=21, cex=1, bg='white')
plot(A2022_det, add=TRUE, col='black', pch=21, cex=1, bg='black')
plot(A2023, add=TRUE, col='black', pch=21, cex=1, bg='white')
plot(A2023_det, add=TRUE, col='black', pch=21, cex=1, bg='black')