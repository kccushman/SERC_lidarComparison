library(terra); library(lidR)

# set working directory
  setwd("/Volumes/KC_JPL/SERC_lidar/GEDI/")

  data2a <- read.csv("data_GEDI2_A.csv")
  data2b <- read.csv('data_GEDI2_B.csv')
  
  years <- unique(data2a$year)

#### data availability ####  
  
# For entire area, sort data into:
  
  # total shots
  # good shots
  # good leaf on
  # good leaf off
  
  # total shots per year
  # good shots per year
  # good leaf on per year
  # good leaf off per year
  
  gediSummary <- data.frame(year = c("All",years),
                            nAll = NA,
                            nGood = NA,
                            nGoodOn = NA)
  
  gediSummary[gediSummary$year=="All","nAll"] <- nrow(data2a)
  gediSummary[gediSummary$year=="All","nGood"] <- nrow(data2a[data2a$quality_flag==1,])
  gediSummary[gediSummary$year=="All","nGoodOn"] <- nrow(data2a[data2a$quality_flag==1 & data2a$leaf_off_flag=="00",])
  
  for(i in 1:length(years)){
    gediSummary[gediSummary$year==years[i],"nAll"] <- nrow(data2a[data2a$year==years[i],])
    gediSummary[gediSummary$year==years[i],"nGood"] <- nrow(data2a[data2a$year==years[i] & data2a$quality_flag==1,])
    gediSummary[gediSummary$year==years[i],"nGoodOn"] <- nrow(data2a[data2a$year==years[i] & data2a$quality_flag==1 & data2a$leaf_off_flag=="00",])
  }
  
  
#### plot data availability ####
  
  # turn data into spatial object
  data2aSp <- vect(data2a, geom=c("lon_lowestmode", "lat_lowestmode"), crs="epsg:4326", keepgeom=FALSE)
  
  # get ha 4 outline
  cat_ha4 <- catalog("/Volumes/KC_JPL/SERC_lidar/ha4/als_ha4_clean.laz")
  extent4_utm <- vect(ext(cat_ha4),crs="epsg:32618")
  extent4_latLon <- project(extent_ha4,"epsg:4326")
  
  # SERC NEON site
  par(mfrow=c(2,2),las=1)
  plot(data2aSp,
       pch=20,
       cex=0.3,
       col=adjustcolor("black",0.2),
       main = "All data")
  plot(data2aSp[data2aSp$quality_flag==1,],
       pch=20,
       cex=0.3,
       col=adjustcolor("black",0.2),
       main = "Good quality data")
  plot(data2aSp[data2aSp$quality_flag==1 & data2a$leaf_off_flag=="00",],
       pch=20,
       cex=0.3,
       col=adjustcolor("black",0.2),
       main = "Good quality leaf on data")
  plot(data2aSp[data2aSp$quality_flag==1 & data2a$leaf_off_flag=="00" & data2a$year=="2021",],
       pch=20,
       cex=0.3,
       col=adjustcolor("black",0.2),
       main = "2021 Good quality leaf on data")
  
  # ha 4
  
  data_ha4 <- mask(data2aSp,extent4_latLon)
  data_ha4_utm <- project(data_ha4, "epsg:32618")
  data_ha4_utmPoly <- buffer(data_ha4_utm,width=12.5)
  
  par(mfrow=c(1,2),las=1)
  plot(data_ha4_utmPoly,
       ext = extent4_utm,
       pch=20,
       cex=1,
       col=adjustcolor("black",0.5),
       main = "All data")
  plot(data_ha4_utmPoly[data_ha4_utmPoly$quality_flag==1 & data_ha4_utmPoly$leaf_off_flag=="00" & data_ha4_utmPoly$year=="2021",],
       ext = extent4_utm,
       pch=20,
       cex=1,
       col=adjustcolor("black",0.5),
       main = "2021 Good quality leaf on data")

  
  