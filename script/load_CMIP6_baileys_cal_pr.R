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
library(rgeos)


filename = '/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/shapefiles/baileys/Baileys_ecoregions_cal.shp'

shp <- readOGR(filename)
shps = spTransform(shp,
                   CRS(
                     "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
                   ))

dir_gcm="/Users/marco/Documents/dati/CMIP5_CMIP6/pr/ssp585/"
listGCMs=file.ls <- list.files(path=dir_gcm,
                               pattern=glob2rx("*CMIP6_*_r1i*1980-2050_hist.nc4"))

# listGCMs = list.files(pattern = "hist")

# 001 ACCESS1-3
# 002 bcc-csm1-1
#
# 037 NorESM1-M
# 038 NorESM1-ME

# pr_day_BCSD_historical_r1i1p1_ACCESS1-0_1950.nc
data(wrld_simpl)

# listGCMs=c("ACCESS1-0","bcc-csm1-1","BNU-ESM","CanESM2","CCSM4","CESM1-BGC","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR","MIROC5","MIROC-ESM","MIROC-ESM-CHEM","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M")
# listGCMs=c("ACCESS1-0","bcc-csm1-1","BNU-ESM","CanESM2","CCSM4","CESM1-BGC","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR","MIROC5","MIROC-ESM","MIROC-ESM-CHEM","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M")
# GFDL-ESM2G
anni_hist = 1981:2000
anni_sim = 2031:2050

for (ifile in 1:length(listGCMs)) {
  ## load ctrl
  namefile = paste0(dir_gcm, listGCMs[ifile])
  obs.nc <- nc_open(namefile)
  obs <- ncvar_get(obs.nc, "pr")
  obs[obs == "1.00000002004088e+20"] <- NA
  obs.nc$dim$lon$vals -> lon
  obs.nc$dim$lat$vals -> lat
  obs = obs[, , 13:dim(obs)[3]]
  ctrl = obs
  rm(obs)
  
  image.plot(lon, lat, apply(ctrl, c(1, 2), mean, na.rm = TRUE))
  plot(wrld_simpl, add = TRUE)
  plot(shps, add = TRUE)
  
  ## load fut
  namefile=paste0(dir_gcm,gsub("_hist.nc4","",listGCMs[ifile]),'_fut.nc4')
  obs.nc <- nc_open(namefile)
  obs <- ncvar_get(obs.nc, "pr")
  obs[obs == "1.00000002004088e+20"] <- NA
  obs = obs[, , 13:dim(obs)[3]]
  sim = obs
  rm(obs)
  
  image.plot(lon, lat, apply(sim, c(1, 2), mean, na.rm = TRUE))
  plot(wrld_simpl, add = TRUE)
  plot(shps, add = TRUE)
  
  
  
  ni = dim(ctrl)[1]
  nj = dim(ctrl)[2]
  points <- expand.grid(lon, lat)
  pts = SpatialPoints(points, proj4string = CRS(proj4string(wrld_simpl)))
  # Se combina el fichero de datos con el .shp y se realiza la media
  # para cada region
  ii <- !is.na(over(pts, shps))
  inout = ii[, 1]
  dim(inout) <- c(nrow(ctrl), ncol(ctrl))
  inout[inout == 0] = NA
  image.plot(lon, lat, inout)
  plot(shps, add = TRUE)
  
  ctrl2 = ctrl
  sim2 = sim
  prec_ctrl = vector()
  prec_sim = vector()
  k = 0
  for (i in 1:length(anni_hist)) {
    for (j in 1:12) {
      k = k + 1
      paste(anni_hist[i], "-", j, "-01", sep = "")
      monthDays(as.Date(paste(anni_hist[i], "-", j, "-01", sep = "")))
      
      ctrl2[, , k] = inout * ctrl[, , k]
      prec_ctrl[k] = mean(ctrl2[, , k], na.rm = TRUE) #* monthDays(as.Date(paste(anni_hist[i], "-", j, "-01", sep = "")))
      
      sim2[, , k] = inout * sim[, , k]
      prec_sim[k] = mean(sim2[, , k], na.rm = TRUE) #* monthDays(as.Date(paste(anni_hist[i], "-", j, "-01", sep = "")))
      
    }
  }
  
  length(prec_sim)
  plot(prec_ctrl)
  
  prec = c(prec_ctrl, prec_sim)
  prets = ts(prec, start = c(1981, 01), frequency = 12)
  
  spi4 = vector()
  dum <-
    spi(
      prets,
      4,
      ref.start = c(1981, 1),
      ref.end = c(2000, 1),
      na.rm = TRUE
    )
  spi4 = dum$fitted
  
  
  #plot(prec*86400000)
  # save(spi6, file = paste0('SPI6_low001_gpcpgrid',ifile-1,'.RData'))
  #namefile_ctrl=paste0(dir_gcm,"tasmax_day_BCSD_historical_r1i1p1_",listGCMs[ifile],"-CAL_1980-2000.nc")
  
  save(spi4, file = paste0(dir_gcm,gsub("_hist.nc4","",listGCMs[ifile]),'_spi4.RData'))
  
}
