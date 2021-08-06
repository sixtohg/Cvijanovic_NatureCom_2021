rm(list = ls())
graphics.off()
gc()

library(raster)
library(rgdal)        # para leer shapefiles
library(maptools)
library(SPEI)

domain='baileys_cal'

filename = '/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/shapefiles/baileys/Baileys_ecoregions_cal.shp';
shp <- readOGR(filename)
shps = spTransform(shp,
                  CRS(
                    "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
                  ))


cal <- readOGR("/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/shapefiles/california/california.shp")
cal = spTransform(cal,
                   CRS(
                     "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
                   ))


plot(cal)
plot(shps,add=TRUE)


cc

area(shps)/1000000000 #133.2879 1000km²
area(cal)/1000000000  #423.967 1000km²

area(shps)/area(cal)

setwd('/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/shapefiles/baileys/');

writeOGR(shps, dsn = "./", layer = "Baileys_ecoregions_cal_wgs84",
         driver = "ESRI Shapefile" )
