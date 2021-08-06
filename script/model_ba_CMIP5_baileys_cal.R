rm(list = ls())
graphics.off()
gc()


library(R.matlab)

# years = 1950:2018

## fix parameters
# dir_sim='/Users/marco/Documents/dati/ivana_ice_fires/Sixto_ECEarth3/'
dir_oss = '/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/'
dir_out = '/Users/marco/Documents/output/ivana_ice_fires/'
source('/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/script/scale_1981_2000.R')

years = 1950:2020
years_model = 1951:2020
iok_model=match(years_model,years)

model_name = 'baileys_cal_tasmax'
domain = 'baileys'


## load data
dum = readMat(paste0(
  dir_oss,
  paste0('fire_1950_2020_', domain, '_mjjas_gt300.mat')
))
BAS = as.numeric(dum$FIRE)

BAS=BAS[iok_model]/1000

load(paste0(dir_oss, paste0("nclimgrid_", domain, "_cal.RData")))
#delete months of the years 1948 and 1949 (used to have complete SPI24 series)
nclimgrid = nclimgrid[,-(1:24)]

# delete years until 1981
nclimgrid = nclimgrid[,-(1:((iok_model[1]-1)*12))]


TX = rep(NA, length(years_model))
for (iyear in 1:length(years_model)) {
  i1 = (iyear - 1) * 12 + 4
  i2 = (iyear - 1) * 12 + 9
  TX[iyear] = mean(nclimgrid[10, i1:i2], na.rm = TRUE)
}
SPI = nclimgrid[4, seq(3, dim(nclimgrid)[2], 12)]




## the regression model
mydata = data.frame(
  "y" = log(BAS),
  "x1" = scale_1981_2000(SPI,years_model),
  "x2" = scale_1981_2000(TX,years_model),
  "x3" = scale_1981_2000(years_model,years_model)
)
fit <- lm(y ~ x1 + x2 + x3 , data = mydata)

# fit <- lm(y ~ x1  + x3 , data = mydata)
# confint.default(fit)
pred = fit$coefficients[1] + fit$coefficients[2] *
  mydata$x1 + fit$coefficients[3] *
  mydata$x2 + fit$coefficients[4] *
  mydata$x3


plot.ts(log(BAS))
lines((pred), col = "red")

print(cor.test(log(BAS),pred))


summary(fit)




dir_gcm_pr="/Users/marco/Documents/dati/CMIP5_CMIP6/pr/rcp85/"
dir_gcm_tasmax="/Users/marco/Documents/dati/CMIP5_CMIP6/tasmax/rcp85/"

listGCMs=c('ACCESS1-3','ACCESS1-3','BCC-CSM1-1','CNRM-CM5','CNRM-CM5','CanESM2','EC-EARTH','EC-EARTH','FGOALS-g2','GFDL-ESM2G',
           'GFDL-ESM2M','inmcm4','inmcm4','IPSL-CM5A-LR','MIROC-ESM','MIROC5','MPI-ESM-LR','MRI-CGCM3','NorESM1-M','NorESM1-M','HadGEM2-ES')


changes=vector()
changes_tasmax=vector()
changes_spi=vector()
years_gcm=1981:2020
pred_gcm=matrix(NA, nrow = length(listGCMs), ncol =length(years_gcm) )

for (ifile in 1:length(listGCMs)) {
  load(paste0(dir_gcm_tasmax,"tasmax_CMIP5_",listGCMs[ifile],"_r1i1p1_rcp85_mon_1980-2050.RData"))
  load(paste0(dir_gcm_pr,"spi4_Amon_mod_rcp85_",sprintf("%03d", ifile),".RData"))
  # namefile = paste0(dir_gcm,
  #                   "spi4_Amon_mod_rcp85_",
  #                   sprintf("%03d", ifile),
  #                   ".RData")
  # 
  # "/Users/marco/Documents/dati/CMIP5_CMIP6/tasmax/rcp85/tasmax_CMIP5_NorESM1-M_r1i1p1_rcp85_mon_1980-2050.RData"
  aux_spi4=vector()
  aux_spi4=spi4[seq(3,length(spi4),12)]
  
  aux_tx=vector()
  for (iyear in 1:40) {
    i1 = (iyear - 1) * 12 + 4
    i2 = (iyear - 1) * 12 + 9
    aux_tx[iyear] = mean(tx[i1:i2], na.rm = TRUE)
  }
  
  
  spi4_3=(aux_spi4)
  tx_4_9=(aux_tx)
  
  writeMat(file.path(dir_gcm_pr, paste0('spi4_3_',ifile-1,'_nostd_',domain,'.mat')), spi4_3=spi4_3)
  writeMat(file.path(dir_gcm_tasmax, paste0('tx_4_9_',ifile-1,'_nostd_',domain,'.mat')), tx_4_9=tx_4_9)
  writeMat(file.path(dir_oss, paste0('CMIP5/spi4_3_',ifile-1,'_nostd_',domain,'.mat')), spi4_3=spi4_3)
  writeMat(file.path(dir_oss, paste0('CMIP5/tx_4_9_',ifile-1,'_nostd_',domain,'.mat')), tx_4_9=tx_4_9)
  
  
  # spi4_3=scale(aux_spi6)
  # tx_4_10=scale(aux_tx)
  # 
  # writeMat(file.path(dir_gcm_pr, paste0('spi3_4_',ifile-1,'_',domain,'.mat')), spi4_3=spi4_3)
  # writeMat(file.path(dir_gcm_tasmax, paste0('tx_4_10_',ifile-1,'_',domain,'.mat')), tx_4_10=tx_4_10)
  # 
  
  
  # pred = fit$coefficients[1] + fit$coefficients[2] * scale(aux)
  pred = fit$coefficients[1] + fit$coefficients[2] *
    scale_1981_2000(aux_spi4,1981:2020) + fit$coefficients[3] *
    scale_1981_2000(aux_tx,1981:2020) 
  
  pred_gcm[ifile,]=pred
  
  changes_spi[ifile]=mean((aux_spi4[21:40]),na.rm=TRUE)-mean((aux_spi4[1:20]),na.rm=TRUE)
  changes_tasmax[ifile]=mean((aux_tx[21:40]),na.rm=TRUE)-mean((aux_tx[1:20]),na.rm=TRUE)
  
  plot.ts(aux_spi4)
  plot.ts(aux_tx)
  
  ctrl_ens=mean(exp(pred[1:20]),na.rm=TRUE)
  ice_ens=mean(exp(pred[21:40]),na.rm=TRUE)
  
  #mean_change.ba1(:,ieco)=nanmean(change_ba1(:,ieco),1);
  #change_ba1(ircm,:)=100*(nanmean(exp(BA_FUT1),1)-nanmean(exp(BA_PAST),1))./nanmean(exp(BA_PAST),1);
  changes[ifile]=100*(ice_ens-ctrl_ens)/ctrl_ens
  
  
  
}



quantile(changes_tasmax,c(0.1,0.5,0.9))
boxplot(changes_tasmax)

quantile(changes_spi,c(0.1,0.5,0.9))
boxplot(changes_spi)

quantile(changes,c(0.1,0.5,0.9))
boxplot(changes)
  

