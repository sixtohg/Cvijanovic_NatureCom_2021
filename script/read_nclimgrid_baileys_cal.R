#####################################################################
##    code to prepare the climate data
#####################################################################
# Se limpa el espacio de trabajo
rm(list = ls())
graphics.off()
gc()
##########################
## PAQUETES Y FUNCIONES ##
##########################
library(SPEI)         # to calculate PET
library(R.matlab)     # to write matlab files
library(ncdf4)
library(fields)
library(maptools)
library(rgdal)        # para leer shapefiles
library(raster)


dir_oss <- "/Users/marco/Documents/dati/obs/nclimgrid"
dir_out = '/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/'
data(wrld_simpl)
domain='baileys_cal'

## load data
fname <- file.path(dir_oss, 'nclimgrid_tmax_cal.nc')
obs.nc <- nc_open(fname)
tmax <- ncvar_get(obs.nc, "tmax")
tmax[tmax == "NaN"] <- NA
obs.nc$dim$lon$vals -> lon
obs.nc$dim$lat$vals -> lat
lat = rev(lat)
tmax = tmax[, ncol(tmax):1, ]
image.plot(lon, lat, apply(tmax, c(1, 2), mean, na.rm = TRUE))
data(wrld_simpl)
plot(wrld_simpl, add = TRUE)

fname <- file.path(dir_oss, 'nclimgrid_tmin_cal.nc')
obs.nc <- nc_open(fname)
tmin <- ncvar_get(obs.nc, "tmin")
tmin[tmin == "NaN"] <- NA
tmin = tmin[, ncol(tmin):1, ]

fname <- file.path(dir_oss, 'nclimgrid_tavg_cal.nc')
obs.nc <- nc_open(fname)
tavg <- ncvar_get(obs.nc, "tavg")
tavg[tavg == "NaN"] <- NA
tavg = tavg[, ncol(tavg):1, ]

fname <- file.path(dir_oss, 'nclimgrid_prcp_cal.nc')
obs.nc <- nc_open(fname)
prcp <- ncvar_get(obs.nc, "prcp")
prcp[prcp == "NaN"] <- NA
prcp = prcp[, ncol(prcp):1, ]

filename = '/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/shapefiles/baileys/Baileys_ecoregions_cal.shp';
shp <- readOGR(filename)
shp = spTransform(shp,
                  CRS(
                    "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
                  ))
# plot(shp)

points <- expand.grid(lon, lat)
pts = SpatialPoints(points, proj4string = CRS(proj4string(wrld_simpl)))

# Se combina el fichero .nc con el .shp y se realiza la media
# para cada region
tx = matrix(data = NA,
            nrow = dim(tmax)[3],
            ncol = dim(shp)[1])

tn = matrix(data = NA,
            nrow = dim(tmin)[3],
            ncol = dim(shp)[1])

ta = matrix(data = NA,
            nrow = dim(tavg)[3],
            ncol = dim(shp)[1])

pr = matrix(data = NA,
            nrow = dim(prcp)[3],
            ncol = dim(shp)[1])



for (ireg in 1:dim(shp)[1]) {
  aux1 = shp[ireg,]@polygons[[1]]@Polygons[[1]]@coords
  inout = point.in.polygon(points[, 1], points[, 2], aux1[, 1], aux1[, 2])
  dim(inout) <- c(nrow(tmax), ncol(tmax))
  inout[inout == 0] = NA
  if (sum(inout, na.rm = TRUE) > 0) {
    image.plot(lon, lat, inout)
    tmax2 = tmax
    tmin2 = tmin
    tavg2 = tavg
    prcp2 = prcp
    for (i in 1:dim(tmax)[3]) {
      tmax2[, , i] = tmax[, , i] * inout
      tx[i, ireg] = mean(tmax2[, , i], na.rm = TRUE)
      
      tmin2[, , i] = tmin[, , i] * inout
      tn[i, ireg] = mean(tmin2[, , i], na.rm = TRUE)
      
      tavg2[, , i] = tavg[, , i] * inout
      ta[i, ireg] = mean(tavg2[, , i], na.rm = TRUE)
      
      
    }
    
    for (i in 1:dim(prcp)[3]) {
      prcp2[, , i] = prcp[, , i] * inout
      pr[i, ireg] = mean(prcp2[, , i], na.rm = TRUE)
      
    }
    
  }
}

nclimgrid = array(data = NA, dim = c(12, dim(tx)[1])) #spi1,2,3,4,5,6,9,12,24,tasmax.tasmin,tas
nclimgrid[10, ] = tx
nclimgrid[11, ] = tn
nclimgrid[12, ] = ta

## calculate SPI
ts = c(1, 2, 3, 4, 5, 6, 9, 12, 24)
prets = ts(pr, start = c(1895, 01), frequency = 12)
for (i in 1:length(ts)) {
  print(ts[i])
  # dum <- spi(prets, ts[i], ref.start = c(1980, 1), ref.end = c(2000, 1), na.rm = TRUE)
  dum <- spi(prets, ts[i], na.rm = TRUE)
  nclimgrid[i, ] = dum$fitted[(length(dum$fitted)-dim(nclimgrid)[2]+1):length(dum$fitted)]
  dum <- NULL
}

save(nclimgrid, file = file.path(dir_out, paste0("nclimgrid_", domain, ".RData")))
writeMat(file.path(dir_out, paste0("nclimgrid_", domain, ".mat")), nclimgrid = nclimgrid)

