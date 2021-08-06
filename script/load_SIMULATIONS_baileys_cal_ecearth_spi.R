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
anni = c(1981:2000) #this is only used to calculate the number of days in a month

dir_oss='/Users/marco/Documents/dati/ivana_ice_fires/Sixto_ECEarth3/'
setwd(dir_oss)
files = list.files(pattern = "pr_allmon_y1030.mat")





for (ifile in 1:length(files)) {
  namefile_ctrl = files[ifile]
  namefile_sim1 = paste0('low_ice/',files[ifile])
  
  prec_ctrl=readMat(namefile_ctrl)
  prec_sim1=readMat(namefile_sim1)
  
  k = 0
  for (i in 1:length(anni)) {
    for (j in 1:12) {
      k = k + 1
      paste(anni[i], "-", j, "-01", sep = "")
      monthDays(as.Date(paste(anni[i], "-", j, "-01", sep = "")))
      prec_ctrl$pr.eco[k] = prec_ctrl$pr.eco[k] * monthDays(as.Date(paste(anni[i], "-", j, "-01", sep = "")))
      prec_sim1$pr.eco[k] = prec_sim1$pr.eco[k] * monthDays(as.Date(paste(anni[i], "-", j, "-01", sep = "")))
    }
  }
  
  prec=c(prec_ctrl$pr.eco,prec_sim1$pr.eco)
  
  plot(prec)
  
  prets=ts(prec, start=c(1981,01), frequency=12)
  
  spi4=vector()
  dum <- spi(prets, 4, ref.start=c(1981,1), ref.end=c(2000,1), na.rm = TRUE)
  spi4 = dum$fitted
  
  save(spi4, file = paste0(sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(namefile_ctrl)),'_spi4.RData'))
  
}

