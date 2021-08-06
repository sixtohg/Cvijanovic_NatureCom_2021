rm(list = ls())
graphics.off()
gc()

library(R.matlab)

## fix parameters
dir_oss = '/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/'
dir_out = '/Users/marco/Documents/output/ivana_ice_fires/'
years = 1950:2020
years_model = 1951:2020
iok_model=match(years_model,years)
model_name = 'baileys_cal_tasmax'
domain = 'baileys_cal'


## load data
load(paste0(dir_oss, paste0("nclimgrid_", domain, ".RData")))
writeMat(file.path(
  dir_oss,
  paste0("nclimgrid_", domain, ".mat")
), nclimgrid = nclimgrid)


#delete months of the years 1948 and 1949
nclimgrid = nclimgrid[, -(1:24)]
# delete years until 1951
nclimgrid = nclimgrid[,-(1:((iok_model[1]-1)*12))]


tx_4_9 = rep(NA, length(years_model))
for (iyear in 1:length(years_model)) {
  i1 = (iyear - 1) * 12 + 4
  i2 = (iyear - 1) * 12 + 9
  tx_4_9[iyear] = mean(nclimgrid[10, i1:i2], na.rm = TRUE)
}
spi4_3 = nclimgrid[4, seq(3, dim(nclimgrid)[2], 12)]





writeMat(file.path(
  dir_oss,
  paste0('spi4_3_nclimgrid.mat')
), spi4_3 = spi4_3)
writeMat(file.path(
  dir_oss,
  paste0('tx_4_9_nclimgrid.mat')
), tx_4_9 = tx_4_9)


