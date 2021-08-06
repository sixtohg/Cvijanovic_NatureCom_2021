rm(list = ls())
graphics.off()
gc()

library(R.matlab)

## fix parameters
source('/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/script/scale_1981_2000.R')
dir_oss = '/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/'
dir_out = '/Users/marco/Documents/output/ivana_ice_fires/'
years = 1950:2020
# years_model = 1981:2020
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


#BAS=BAS+1

plot.ts(BAS)
plot.ts(log(BAS))
# BAS[1]=NA
BAS=BAS[iok_model]


load(paste0(dir_oss, paste0("nclimgrid_", domain, "_cal.RData")))
#delete months of the years 1948 and 1949 (used to have complete SPI24 series)
nclimgrid = nclimgrid[,-(1:24)]

# delete years until 1951
nclimgrid = nclimgrid[,-(1:((iok_model[1]-1)*12))]


# nclimgrid[10, ] = tx
# nclimgrid[11, ] = tn
# nclimgrid[12, ] = ta
# ts = c(1, 2, 3, 4, 5, 6, 9, 12, 24)






tx = rep(NA, length(years_model))
for (iyear in 1:length(years_model)) {
  i1 = (iyear - 1) * 12 + 5
  i2 = (iyear - 1) * 12 + 8
  tx[iyear] = mean(nclimgrid[12, i1:i2], na.rm = TRUE)
}
SPI_aux = nclimgrid[7, seq(7, dim(nclimgrid)[2], 12)]

## all the series
mydata_reg = data.frame(
  "y" = log(BAS/1000),
  "x1" = years_model,
  "x2" = years_model^2,
  "x3" = scale_1981_2000(tx,years_model),
  "x4" = scale_1981_2000(SPI_aux,years_model)
)
fit <- lm(y ~ x1 + x2 + x3 + x4, data = mydata_reg)
# confint.default(fit)
summary(fit)
pred = fit$coefficients[1] + fit$coefficients[2] *
  mydata_reg$x1 + fit$coefficients[3] *
  mydata_reg$x2 + fit$coefficients[4] *
  mydata_reg$x3 + fit$coefficients[5] *
  mydata_reg$x4

# dev.off()
# par(mfrow=c(1,3))
# plot(log(BAS))
# plot.ts(log(BAS))
plot.ts(log(BAS))
# lines(log(BAS))
lines((pred), col = "red")
print(cor.test(log(BAS), pred))



cc


## find the best predictor
num_iter = 3645
corr = array(NA, dim = c(num_iter))
sig = array(NA, dim = c(num_iter))
best_m_SPI = array(NA, dim = c(num_iter))
best_t_SPI = array(NA, dim = c(num_iter))
best_start_tx = array(NA, dim = c(num_iter))
best_stop_tx = array(NA, dim = c(num_iter))
best_x1 = array(NA, dim = c(num_iter, length(years)))
best_x2 = array(NA, dim = c(num_iter, length(years)))
best_x3 = array(NA, dim = c(num_iter, length(years)))

best_x1_1 = array(NA, dim = c(num_iter, length(years)))
best_x1_12 = array(NA, dim = c(num_iter, length(years)))
best_x1_13 = array(NA, dim = c(num_iter, length(years)))
best_x1_123 = array(NA, dim = c(num_iter, length(years)))

best_x2_12 = array(NA, dim = c(num_iter, length(years)))
best_x2_123 = array(NA, dim = c(num_iter, length(years)))

best_x3_13 = array(NA, dim = c(num_iter, length(years)))
best_x3_123 = array(NA, dim = c(num_iter, length(years)))


years = scale(years_model)

# climdiv = array(data = NA,dim = c(12,dim(MyData)[1])) #spi1,2,3,4,5,6,9,12,24,tasmax.tasmin,tas
#iscC = c(3, 6, 8)
iscC = c(1, 2, 3, 4, 5, 6, 7, 8, 9)


start_time <- Sys.time()

k = 0
for (im in 1:9) {
  tmp_stop = im:9
  for (istop in 1:length(tmp_stop)) {
    tx = rep(NA, length(years))
    for (iyear in 1:length(years)) {
      i1 = (iyear - 1) * 12 + im
      i2 = (iyear - 1) * 12 + tmp_stop[istop]
      tx[iyear] = mean(nclimgrid[12, i1:i2], na.rm = TRUE)
    }
    
    tx = scale(tx)
    for (im_spi in 1:9) {
      for (isc in 1:length(iscC)) {
        k = k + 1
        print(k)
        
       
        SPI_aux = scale(nclimgrid[iscC[isc], seq(im_spi, dim(nclimgrid)[2], 12)])
        
        pre1 = vector()
        pre12 = vector()
        pre13 = vector()
        pre123 = vector()
        
        for (iy in 1:length(years)) {
          BAS_train = BAS[-iy]
          SPI_train = SPI_aux[-iy]
          tx_train = tx[-iy]
          years_train = years[-iy]
          
          SPI_test = SPI_aux[iy]
          tx_test = tx[iy]
          years_test = years[iy]
          
          mydata_train = data.frame(
            "y" = log(BAS_train),
            "x0" = (years_train),
            "x1" = (years_train^2),
            "x2" = (tx_train),
            "x3" = (SPI_train)
          )
          
          mydata_test = data.frame(
            "x0" = (years_test),
            "x1" = (years_test^2),
            "x2" = (tx_test),
            "x3" = (SPI_test)
          )
          
          fit1 <- lm(y ~ x0 + x1 , data = mydata_train)
          fit12 <- lm(y ~ x0 + x1 + x2 , data = mydata_train)
          fit13 <- lm(y ~ x0 + x1 + x3, data = mydata_train)
          fit123 <- lm(y ~ x0 + x1 + x2 + x3, data = mydata_train)
          
          pre1[iy] = fit1$coefficients[1] + fit1$coefficients[2] * mydata_test$x0 + fit1$coefficients[3] * mydata_test$x1 
          pre12[iy] = fit12$coefficients[1] + fit12$coefficients[2] * mydata_test$x0 + fit12$coefficients[3] * mydata_test$x1 + fit12$coefficients[4] * mydata_test$x2
          pre13[iy] = fit13$coefficients[1] + fit13$coefficients[2] * mydata_test$x0 + fit13$coefficients[3] * mydata_test$x1 + fit13$coefficients[4] * mydata_test$x3
          pre123[iy] = fit123$coefficients[1] + fit123$coefficients[2] * mydata_test$x0 + fit123$coefficients[3] * mydata_test$x1 + fit123$coefficients[4] * mydata_test$x2 + fit123$coefficients[5] * mydata_test$x3
          
          # best_x1_1[k, iy] = fit1$coefficients[2]
          # best_x1_12[k, iy] = fit12$coefficients[2]
          # best_x1_13[k, iy] = fit13$coefficients[2]
          # best_x1_123[k, iy] = fit123$coefficients[2]
          # 
          # best_x2_12[k, iy] = fit12$coefficients[3]
          # best_x2_123[k, iy] = fit123$coefficients[3]
          # 
          # best_x3_13[k, iy] = fit13$coefficients[3]
          # best_x3_123[k, iy] = fit123$coefficients[4]
          
          
        }
        
        rho1 = cor.test(log(BAS),
                        pre1,
                        use = "pairwise.complete.obs",
                        alternative = "greater")
        rho12 = cor.test(log(BAS),
                         pre12,
                         use = "pairwise.complete.obs",
                         alternative = "greater")
        rho13 = cor.test(log(BAS),
                         pre13,
                         use = "pairwise.complete.obs",
                         alternative = "greater")
        rho123 = cor.test(log(BAS),
                          pre123,
                          use = "pairwise.complete.obs",
                          alternative = "greater")
        if (rho1$estimate > rho12$estimate &
            rho1$estimate > rho13$estimate &
            rho1$estimate > rho123$estimate) {
          corr[k] = rho1$estimate
          sig[k] = rho1$p.value
          # best_x1[k, ] = best_x1_1[k, ]
        }
        if (rho12$estimate > rho1$estimate &
            rho12$estimate > rho13$estimate &
            rho12$estimate > rho123$estimate) {
          corr[k] = rho12$estimate
          sig[k] = rho12$p.value
          # best_x1[k, ] = best_x1_12[k, ]
          # best_x2[k, ] = best_x2_12[k, ]
          best_start_tx[k] = im
          best_stop_tx[k] = tmp_stop[istop]
        }
        if (rho13$estimate > rho1$estimate &
            rho13$estimate > rho12$estimate &
            rho13$estimate > rho123$estimate) {
          corr[k] = rho13$estimate
          sig[k] = rho13$p.value
          # best_x1[k, ] = best_x1_13[k, ]
          # best_x3[k, ] = best_x3_13[k, ]
          best_t_SPI[k] = isc
          best_m_SPI[k] = im_spi
        }
        if (rho123$estimate > rho1$estimate &
            rho123$estimate > rho12$estimate &
            rho123$estimate > rho13$estimate) {
          corr[k] = rho123$estimate
          sig[k] = rho123$p.value
          # best_x1[k, ] = best_x1_123[k, ]
          # best_x2[k, ] = best_x2_123[k, ]
          # best_x3[k, ] = best_x3_123[k, ]
          best_start_tx[k] = im
          best_stop_tx[k] = tmp_stop[istop]
          best_t_SPI[k] = isc
          best_m_SPI[k] = im_spi
          
        }
        
      }
    }
  }
}

end_time <- Sys.time()

end_time - start_time

sig = p.adjust(sig, method = "fdr")
corr[sig > 0.05] = NA



dum = max((corr), na.rm = TRUE)
idx = which(abs(corr) == dum, arr.ind = TRUE)

aux_corr=sort(corr, decreasing=TRUE,index.return=TRUE)
aux_corr$ix
aux_corr$x

# best_t_SPI[aux_corr$ix[1:10]]
# best_m_SPI[aux_corr$ix[1:10]]
# best_x3[aux_corr$ix[1:10]]
# aux_corr$x[1:10]
best_t_SPI[aux_corr$ix]
corr[aux_corr$ix]


best_corr = corr[idx]
best_m_fin = best_m_SPI[idx]
best_t_SPI_fin = best_t_SPI[idx]
best_start_tx_fin = best_start_tx[idx]
best_stop_tx_fin = best_stop_tx[idx]
best_sig = sig[idx]

# best_x1_fin = best_x1[idx,]
# best_x2_fin = best_x2[idx,]
# best_x3_fin = best_x3[idx,]


# save(best_x1_fin, file = file.path(dir_out, paste0('best_x1_fin_', model_name, '.RData')))
# save(best_x2_fin, file = file.path(dir_out, paste0('best_x2_fin_', model_name, '.RData')))
# save(best_x3_fin, file = file.path(dir_out, paste0('best_x3_fin_', model_name, '.RData')))

rm(tx)
tx = rep(NA, length(years))
for (iyear in 1:length(years)) {
  i1 = (iyear - 1) * 12 + best_start_tx_fin
  i2 = (iyear - 1) * 12 + best_stop_tx_fin
  tx[iyear] = mean(nclimgrid[12, i1:i2], na.rm = TRUE)
}
SPI_aux = nclimgrid[best_t_SPI_fin, seq(best_m_fin, dim(nclimgrid)[2], 12)]

## all the series
mydata_reg = data.frame(
  "y" = log(BAS/1000),
  "x1" = years_model,
  "x2" = years_model^2,
  "x3" = scale_1981_2000(tx,years_model),
  "x4" = scale_1981_2000(SPI_aux,years_model)
)
fit <- lm(y ~ x1 + x2 + x3 + x4, data = mydata_reg)
# confint.default(fit)
summary(fit)
pred = fit$coefficients[1] + fit$coefficients[2] *
  mydata_reg$x1 + fit$coefficients[3] *
  mydata_reg$x2 + fit$coefficients[4] *
  mydata_reg$x3 + fit$coefficients[5] *
  mydata_reg$x4

# dev.off()
# par(mfrow=c(1,3))
# plot(log(BAS))
# plot.ts(log(BAS))
plot.ts(log(BAS))
# lines(log(BAS))
lines((pred), col = "red")
print(cor.test(log(BAS), pred))




AIC(fit)


# save(best_corr, file = file.path(dir_out, paste0('corr_', model_name, '.RData')))
# save(best_sig, file = file.path(dir_out, paste0('sig_', model_name, '.RData')))
# save(best_m_fin, file = file.path(dir_out, paste0('best_m_', model_name, '.RData')))
# save(best_t_SPI_fin, file = file.path(dir_out, paste0('best_t_SPI_', model_name, '.RData')))
# # save(best_t_tmx_fin, file = file.path(dir_out, paste0('best_t_tmax_',model_name,'.RData')))
# save(best_x1_fin, file = file.path(dir_out, paste0('best_x1_', model_name, '.RData')))
# save(best_x2_fin, file = file.path(dir_out, paste0('best_x2_', model_name, '.RData')))
