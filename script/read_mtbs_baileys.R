rm(list = ls())
graphics.off()
gc()

require(rgdal)
library(R.matlab)
library(raster)
library(rgeos)

dir_fire<-"/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/"

# The input file geodatabase
mtbs <- readOGR("/Users/marco/Documents/dati/fire_us/MTBS_1984_2019/mtbs_perims_DD/mtbs_perims_DD_baileys.shp", verbose = FALSE)
summary(mtbs)
mtbs = spTransform(mtbs,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))




#plot(rr)

BA = array(data = 0, dim = c(1, length(1984:2019) * 12))

anni=as.numeric(format(as.Date(mtbs@data$Ig_Date),"%Y"))
mesi=as.numeric(format(as.Date(mtbs@data$Ig_Date),"%m"))

k = 0
for (iyear in 1984:2019) {
  for (imonth in 1:12) {
    k = k + 1
    idx2 = which(anni == iyear &
                   mesi == imonth)
    
    
    
    
    if (length(idx2) >= 1) {
      
      
      
      for (iidx2 in 1:length(idx2)) {
        
        a=unlist(mtbs[idx2[iidx2],]@polygons)
        
        
        BA[1, k] = sum(BA[1, k],a[[1]]@Polygons[[1]]@area,na.rm = TRUE)
        
        
        
        #BA[1, k] = sum(as.numeric(mtbs[idx2,]$BurnBndAc))
        
      } 
    }
  }
}




save(BA, file = file.path(dir_fire,"fire_1984_2019_monthly_baileys.RData"))
writeMat(file.path(dir_fire, paste("fire_1984_2019_monthly_baileys.mat",sep="")), BA=BA)

