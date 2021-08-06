rm(list = ls())
graphics.off()
gc()

library(ncdf4)
library(fields)
library(maps)
library(Hmisc)
library(sp)
library(R.matlab)
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
# shp = spTransform(shp,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# shps <- gSimplify(shp, tol=0.01, topologyPreserve=TRUE)
# plot(shps)



setwd("/Users/marco/Documents/dati/ivana_ice_fires/")
files = list.files(pattern = "cam.2040y_ptt.nc")

#PREC_ctrl.cam.2040y_ptt_gpcpgrid0_CALIFORNIA.RData

data(wrld_simpl)



for (ifile in 1:6) {
  if (ifile == 1) {
    namefile_ctrl = paste0('ctrl.cam.2040y_ptt.nc')
    namefile_sim1 = paste0('extlow001_iv3p.cam.2040y_ptt.nc')
    namefile_sim2 = paste0('extlow002_iv3p.cam.2040y_ptt.nc')
  } else {
    namefile_ctrl = paste0('ctrl_0',ifile-1,'.cam.2040y_ptt.nc')
    namefile_sim1 = paste0('extlow001_iv3p_0',ifile-1,'.cam.2040y_ptt.nc')
    namefile_sim2 = paste0('extlow002_iv3p_0',ifile-1,'.cam.2040y_ptt.nc')
  }
  
  
  
  ## load ctrl
  obs.nc <- nc_open(namefile_ctrl)
  obs <- ncvar_get(obs.nc, "TSMX")
  # obs[obs == "1.e+20f"] <- NA
  obs.nc$dim$lon$vals -> lon
  obs.nc$dim$lat$vals -> lat
  # lat <- rev(lat)
  # obs = obs[, ncol(obs):1, ]
  lon=lon-180
  dum=obs
  dum[1:(nrow(obs)/2),,]=obs[(nrow(obs)/2+1):nrow(obs),,]
  dum[(nrow(obs)/2+1):nrow(obs),,]=obs[1:(nrow(obs)/2),,]
  rm(obs)
  ctrl=dum[,,13:dim(dum)[3]]
  rm(dum)
  
  ## load sim1
  obs.nc <- nc_open(namefile_sim1)
  obs <- ncvar_get(obs.nc, "TSMX")
  # obs[obs == "1.e+20f"] <- NA
  obs.nc$dim$lon$vals -> lon
  obs.nc$dim$lat$vals -> lat
  # lat <- rev(lat)
  # obs = obs[, ncol(obs):1, ]
  lon=lon-180
  dum=obs
  dum[1:(nrow(obs)/2),,]=obs[(nrow(obs)/2+1):nrow(obs),,]
  dum[(nrow(obs)/2+1):nrow(obs),,]=obs[1:(nrow(obs)/2),,]
  rm(obs)
  sim1=dum[,,13:dim(dum)[3]]
  rm(dum)
  
  ## load sim2
  obs.nc <- nc_open(namefile_sim2)
  obs <- ncvar_get(obs.nc, "TSMX")
  # obs[obs == "1.e+20f"] <- NA
  obs.nc$dim$lon$vals -> lon
  obs.nc$dim$lat$vals -> lat
  # lat <- rev(lat)
  # obs = obs[, ncol(obs):1, ]
  lon=lon-180
  dum=obs
  dum[1:(nrow(obs)/2),,]=obs[(nrow(obs)/2+1):nrow(obs),,]
  dum[(nrow(obs)/2+1):nrow(obs),,]=obs[1:(nrow(obs)/2),,]
  rm(obs)
  sim2=dum[,,13:dim(dum)[3]]
  rm(dum)
  
  ni = dim(ctrl)[1]
  nj = dim(ctrl)[2]
  points <- expand.grid(lon, lat)
  pts = SpatialPoints(points, proj4string = CRS(proj4string(wrld_simpl)))
  
  image.plot(lon, lat, apply(ctrl, c(1, 2), mean, na.rm = TRUE)-273.15)
  plot(wrld_simpl,add=TRUE)
  #map("world2", add = TRUE)
  plot(shps,add=TRUE)
  
  
  
  
  # Se combina el fichero de datos con el .shp y se realiza la media
  # para cada region
  # ii <- !is.na(over(pts, shps))
  # # inout = ii #ii[, 1]
  # inout = ii[, 1]
  # dim(inout) <- c(nrow(ctrl), ncol(ctrl))
  # inout[inout == 0] = NA
  # image.plot(lon, lat, inout)
  
  aux1 = shps[1,]@polygons[[1]]@Polygons[[1]]@coords
  inout = point.in.polygon(points[, 1], points[, 2], aux1[, 1], aux1[, 2])
  dim(inout) <- c(nrow(ctrl), ncol(ctrl))
  inout[inout == 0] = NA
  image.plot(lon, lat, inout)
  
  ctrl2=ctrl
  sim12=sim1
  sim22=sim2
  tx_ctrl=vector()
  tx_sim1=vector()
  tx_sim2=vector()
  for (i in 1:dim(ctrl)[3]) {
    ctrl2[,,i]=inout*ctrl[,,i]
    tx_ctrl[i]=mean(ctrl2[, , i], na.rm = TRUE)-273.15
    
    sim12[,,i]=inout*sim1[,,i]
    tx_sim1[i]=mean(sim12[, , i], na.rm = TRUE)-273.15
    
    sim22[,,i]=inout*sim2[,,i]
    tx_sim2[i]=mean(sim22[, , i], na.rm = TRUE)-273.15
  }
    
  
  ilon=which(lon> -140 & lon < -100)
  ilat=which(lat> 10 & lat < 50)
  
  
  image.plot(lon[ilon], lat[ilat], apply(ctrl2[ilon,ilat,], c(1, 2), mean, na.rm = TRUE)-273.15)
  plot(wrld_simpl,add=TRUE)
  #map("world2", add = TRUE)
  plot(shps,add=TRUE)
  
  
  tx=c(tx_ctrl,tx_sim1)
  
  
  save(tx, file = paste0('TX_low001_0',ifile-1,'.cam.2040y_ptt_',domain,'.RData'))
  
  tx=c(tx_ctrl,tx_sim2)
  
  save(tx, file = paste0('TX_low002_0',ifile-1,'.cam.2040y_ptt_',domain,'.RData'))
  
}

