rm(list = ls())
graphics.off()
gc()

require(rgdal)
library(rgeos)
library(secr)
library(R.matlab)
library(dplyr)
library(readr)
library(ggplot2)
library(ggthemes)
library(scales)
library(maps)
library(mapproj)

# color palette for major fire causes
cause_pal <- c("#ffff00","#d397fc","#ffffff")

# load and process data
## load data
mtbs <- readOGR("/Users/marco/Documents/dati/fire_us/fire20_1_baileys.shp", verbose = FALSE)
summary(mtbs)
mtbs = spTransform(mtbs,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
iok=!is.na(mtbs@data$ALARM_D)
mtbs=mtbs[iok,]
anni=as.numeric(format(as.Date(mtbs@data$ALARM_D),"%Y"))

iok=anni>=1951 & anni<=2020


mtbs=mtbs[iok,]

anni=as.numeric(format(as.Date(mtbs@data$ALARM_D),"%Y"))
mesi=as.numeric(format(as.Date(mtbs@data$ALARM_D),"%m"))
giorni=as.numeric(format(as.Date(mtbs@data$ALARM_D),"%d"))
fires=mtbs@data$GIS_ACR



BAS=fires
year_=anni

plot_date1 = as.Date(mtbs@data$ALARM_D)
plot_date=as.Date(format(plot_date1,"2020-%m-%d"))

calfire <- data.frame(year_, BAS,as.Date(plot_date, format =  "2020-%m-%d"))


name <- c("year_", "gis_acres", "plot_date")

colnames(calfire) <- name
print (calfire)

calfire[,1] <- as.numeric(calfire[,1])


lapply(calfire, class)

# pdf(file = "/Users/marco/Documents/dati/fire_us/calfire_frap.pdf")
# setEPS()
# postscript("/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/paper/figures_temp/annual_cycle.eps")
# pdf("/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/paper/figures_temp/annual_cycle.pdf")
cairo_ps(file = "/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/paper/figures_temp/annual_cycle.eps", onefile = FALSE, fallback_resolution = 600)
# plot template
plot_template <- ggplot(calfire, aes(y=year_)) +
  geom_hline(yintercept = seq(1951, 2020, by = 1), color = "gray", size = 0.05) +
  scale_size_area(max_size = 15, guide = FALSE) +
  scale_x_date(date_breaks = "months", date_labels = "%b") +
  # scale_y_reverse(limits = c(2017,1950), breaks = c(2010,1990,1970,1950)) +
  scale_y_continuous(limits = c(1951,2020), breaks = c(1951,1960,1970,1980,1990,2000,2010,2020)) +
  xlab("") +
  ylab("")  +
  # geom_vline(xintercept = as.numeric(as.Date("2020-05-01")), linetype=4)
# theme_hc(bgcolor = "darkunica", base_size = 20) +
  theme_hc(bgcolor = "darkunica", base_size = 20, base_family = "Arial") +
  # theme_hc() +
  theme(axis.text = element_text(color = "#ffffff"))

plot_template +
  geom_point(aes(size=gis_acres, x=plot_date), color="#ffa500", alpha=0.7)
# dev.copy2pdf(file="/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/paper/figures_temp/annual_cycle.pdf",out.type="cairo", width=10, height=7.5)
dev.off()

# dev.off()