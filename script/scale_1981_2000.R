scale_1981_2000 <- function(datain,years) {
  period=1981:2000
  iok=match(period, years)
  avg=mean(datain[iok],na.rm=TRUE)
  std=sd(datain[iok],na.rm=TRUE)
  dataout=(datain-avg)/std
}