library(terra); library(lidR); library(sf)

# set working directory
  setwd("/Volumes/KC_JPL/SERC_lidar/GEDI/")

  data2a <- read.csv("data_GEDI2_A.csv")
    data2a$shot_number <- read.csv("data_GEDI2_A_shot_number.csv",colClasses = "character")[,1]
  data2b <- read.csv('data_GEDI2_B.csv')
    data2b$shot_number <- read.csv("data_GEDI2_B_shot_number.csv",colClasses = "character")[,1]
  data1b <- read.csv('data_GEDI1_B.csv')
    data1b$shot_number <- read.csv("data_GEDI1_B_shot_number.csv",colClasses = "character")[,1]
  
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
  extent4_latLon <- project(extent4_utm,"epsg:4326")
  
  # SERC NEON site
mainCex=0.7
labCex=0.8

jpeg(filename = "GEDI_sercSite.jpeg",
     width = 1100, height = 1200, units = "px", pointsize = 36,
     quality = 300)
  par(mfrow=c(2,2),las=1,mar=c(0,0,0,0),oma=c(2,2,1,1))
  plot(data2aSp,
       mar=c(1,1,1,1),
       pax=list(lab=2),
       pch=20,
       cex=0.3,
       las=1,
       col=adjustcolor("black",0.2))
  mtext("All data", side=3, outer=F, cex = mainCex)
  plot(data2aSp[data2aSp$quality_flag==1,],
       mar=c(1,1,1,1),
       pax=list(lab=0),
       pch=20,
       cex=0.3,
       las=1,
       col=adjustcolor("black",0.2))
  mtext("Good quality data", side=3, outer=F, cex = mainCex)
  plot(data2aSp[data2aSp$quality_flag==1 & data2a$leaf_off_flag=="00",],
       mar=c(1,1,1,1),
       pax=list(lab=1:2),
       pch=20,
       cex=0.3,
       las=1,
       col=adjustcolor("black",0.2))
  mtext("Good quality leaf on data", side=3, outer=F, cex = mainCex)
  
  plot(data2aSp[data2aSp$quality_flag==1 & data2a$leaf_off_flag=="00" & data2a$year=="2021",],
       mar=c(1,1,1,1),
       pax=list(lab=1),
       pch=20,
       cex=0.3,
       las=1,
       col=adjustcolor("black",0.2))
  mtext("2021 Good quality leaf on data", side=3, outer=F, cex = mainCex)
  mtext("Latitude (deg)", side=2, outer=T, cex = labCex, las=0, line=1)
  mtext("Longitude (deg)", side=1, outer=T, cex = labCex, line=0)
  
dev.off()

  # ha 4
  
  data_ha4 <- mask(data2aSp,extent4_latLon)
  data_ha4_utm <- project(data_ha4, "epsg:32618")
  data_ha4_utmPoly <- buffer(data_ha4_utm,width=12.5)
  data_ha4_utmPoly_good <- data_ha4_utmPoly[data_ha4_utmPoly$quality_flag==1 & data_ha4_utmPoly$leaf_off_flag=="00" & data_ha4_utmPoly$year=="2021",]
  
jpeg(filename = "GEDI_ha4.jpeg",
     width = 1100, height = 500, units = "px", pointsize = 36,
     quality = 300)
  par(mfrow=c(1,2),las=1,mar=c(0,0,0,0),oma=c(2,2,0,0))
  plot(data_ha4_utmPoly,
       ext = extent4_utm,
       mar=c(1,1,1,1),
       pax=list(lab=1:2),
       las=1,
       pch=20,
       cex=1,
       col=adjustcolor("black",0.5))
  mtext("All data", side=3, outer=F, cex = mainCex)
  
  plot(data_ha4_utmPoly_good,
       ext = extent4_utm,
       mar=c(1,1,1,1),
       pax=list(lab=1),
       pch=20,
       cex=1,
       col=adjustcolor("black",0.5))
  mtext("2021 Good quality leaf on data", side=3, outer=F, cex = mainCex)
  mtext("Easting (m)", side=2, outer=T, cex = labCex, las=0, line=1)
  mtext("Northing (m)", side=1, outer=T, cex = labCex, line=0)
  
dev.off()  
#### clip drone lidar over the same area as GEDI data ####
  droneCtg <- catalog("/Volumes/KC_JPL/SERC_lidar/GatorEye/Original/")
  for(i in 1:length(data_ha4_utmPoly_good)){
    GEDI_clip <- clip_roi(droneCtg, st_as_sf(data_ha4_utmPoly_good[i,]))
    writeLAS(GEDI_clip,paste0("GEDI_example",i,".laz"))
  }
  
#### read and plot full waveform for comparison data ####
library(rGEDI)
  
goodShot <- data_ha4_utmPoly_good$shot_number[1]

# find correct file to read for 1b waveform
goodFile <- data1b[data1b$shot_number==goodShot,"fn"]
elevation <- data2a[data2a$shot_number==goodShot,"elev_lowestmode"]
  
data1B <- readLevel1B(goodFile) # both good ha 4 shots are in the same .h5 file

data1B_shot <- getLevel1BWF(data1B, shot_number=goodShot)    
  
data1B_amplitude <- data1B_shot@dt$rxwaveform
  # rescale amplitude values between 0 and 1
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  data1B_amplitudeScaled <- range01(data1B_amplitude)
  
data1B_elevation <- data1B_shot@dt$elevation
  # rescale elevation relative to ground height
  data1B_htAboveGround <- data1B_shot@dt$elevation - elevation 
  
# get min and max Z values from drone lidar
droneShot <- readLAS("GEDI_example1.laz")

lineWidth <- 4
jpeg(filename = "GEDI_example.jpeg",
     width = 1800, height = 1200, units = "px", pointsize = 36,
     quality = 300)
  par(mar=c(3,3,1,1), oma=c(2,2,0,0), mfrow=c(1,2), las=1)
  plot(x = data1B_amplitudeScaled,
       y = data1B_elevation,
       xlab=NA,ylab=NA,
       type = "l", lwd=lineWidth)
  abline(h=c(min(droneShot$Z),quantile(droneShot$Z,0.98)),
         lwd=lineWidth+2,col=adjustcolor("red",0.6))
  abline(h=c(elevation,data_ha4_utmPoly_good$rh_98[1] + elevation),
         lwd=lineWidth+2, col=adjustcolor("blue",0.6))
  plot(x = data1B_amplitudeScaled,
       y = data1B_elevation,
       ylim=c(-40,10),
       xlab=NA,ylab=NA,
       type = "l",lwd=lineWidth)
  abline(h=c(min(droneShot$Z),quantile(droneShot$Z,0.98)),
         lwd=lineWidth+2,col=adjustcolor("red",0.6))
  abline(h=c(elevation,data_ha4_utmPoly_good$rh_98[1] + elevation),
         lwd=lineWidth+2, col=adjustcolor("blue",0.6))
  mtext("Scaled waveform amplitude", side=1, outer=T, cex=1.5)
  mtext("Elevation (m)", side=2, outer=T, las=0, cex=1.5, line=0.5)
dev.off()  