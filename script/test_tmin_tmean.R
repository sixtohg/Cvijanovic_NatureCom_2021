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
# savename = [dir_fire,'fire_1950_2020_',domain,'_mjjas_gt300.mat'];
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
  TX[iyear] = mean(nclimgrid[12, i1:i2], na.rm = TRUE)
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