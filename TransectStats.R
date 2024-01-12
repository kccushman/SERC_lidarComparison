library("lidR")
library("terra")

#### setup ####
rootDir <- "/Volumes/KC_JPL/SERC_lidar/"

# Define file path for MLS, ALS, and drone lidar files
mlsFile <- paste0(rootDir,"transect_mls.laz")
alsFile <- paste0(rootDir,"transect_als.laz")
droneFile <- paste0(rootDir,"transect_drone.laz")
tlsFile <- paste0(rootDir,"transect_tls.laz")

# Make lidR catalog objects
mlsCat <- catalog(mlsFile)
alsCat <- catalog(alsFile)
droneCat <- catalog(droneFile)
tlsCat <- catalog(tlsFile)

#### make terrain model ####

# use previously classified ALS data for whole area
  alsAll <- catalog(paste0(rootDir,"als_ha4_clean.laz"))
  dtm = rasterize_terrain(alsAll, res = 1, algorithm = knnidw(k = 6L, p = 2))

#### summary stats of transects ####
  
#  Make raster template for voxelized point count/density
  
  # define voxel size
  voxelSz <- 5
  
  # create template
  rastTemplate <- rast(nrows=45/voxelSz, ncols = 80/voxelSz,
                       xmin=0,xmax=80,
                       ymin=0,ymax=45)
  
  # define function to calculate PAI using MacArthur Horn method
  MH_PAI <- function(X,Y,Z,N,type="air"){
    
    # combine into data frame
    dataPoints <- data.frame(X,Y,Z,N)
    dataPoints$XY <- paste(X,Y,sep="-")
    dataPoints$ePAI <- NA
    
    # get unique X,Y columns of voxels
    uniqueXY <- unique(dataPoints$XY)
    
    for(i in 1:length(uniqueXY)){
      data_i <- dataPoints[dataPoints$XY==uniqueXY[i],]
      
      # order by height depending on observation type
      if(type=="air"){
        # "top to bottom" for airborne
        data_i <- data_i[order(-data_i$Z),]
      }
      if(type=="ground"){
        # "bottom to top" for ground
        data_i <- data_i[order(data_i$Z),]
      }
      
      for(j in 1:nrow(data_i)){
        pulseIn <- sum(data_i$N[j:nrow(data_i)])
        pulseOut <- sum(data_i$N[(j+1):nrow(data_i)])
        if(j==nrow(data_i)){
          pulseOut <- 0
        }
        dataPoints[dataPoints$XY==uniqueXY[i] & dataPoints$Z==data_i$Z[j],"ePAI"] <- log(pulseIn/pulseOut)
        if(pulseIn==0|pulseOut==0){
          dataPoints[dataPoints$XY==uniqueXY[i] & dataPoints$Z==data_i$Z[j],"ePAI"] <- 0
        }
      }
      
    }
    
    return(dataPoints)
  }
  
  
### ALS
  
    # read lidar data
    data <- readLAS(alsFile)
    # subtract ground height
    dataNorm <- normalize_height(data, dtm)

    # 2D POINT DENSITY
      
      # calculate voxel point density
      voxel_als <- voxel_metrics(dataNorm, ~list(N = length(Z)), voxelSz, all_voxels = T)
      # replace "NA" values with 0
      voxel_als[is.na(voxel_als$N),"N"] <- 0
      # creates 2 "layers", recombine into one
      voxel_als <- aggregate(N~X+Z, data = voxel_als, sum)
      # scale to pts/m3
      voxel_als$pts_m3 <- voxel_als$N/(voxelSz^3)
      # rescale X values
      voxel_als$Xplot <- voxel_als$X - 364560
      voxel_als$Yplot <- voxel_als$Z + 2.5
      voxel_als_vect <- vect(voxel_als,
                             geom=c("Xplot", "Yplot"))
      densRast_als <- terra::rasterize(voxel_als_vect,rastTemplate, field="pts_m2")
    
    # CANOPY HEIGHT
      chm_als <- grid_canopy(dataNorm, res=1,
                             algorithm = p2r(subcircle=0.01, 
                                             na.fill = tin()))
      
    # Voxelized effective PAI
      
      # only keep first returns, and remove points < 1 m height
      dataNormVeg <- dataNorm[dataNorm$Z>=2 & dataNorm$ReturnNumber==1,]
      # re-calculate voxel point density
      voxel_als_pai <- voxel_metrics(dataNormVeg, ~list(N = length(Z)), voxelSz, all_voxels = T)
      # replace "NA" values with 0
      voxel_als_pai[is.na(voxel_als_pai$N),"N"] <- 0
      # creates 2 "layers", recombine into one
      voxel_als_pai <- aggregate(N~X+Z, data = voxel_als_pai, sum)
      # scale to pts/m3
      voxel_als_pai[voxel_als_pai$Z>0,"pts_m3"] <- voxel_als_pai[voxel_als_pai$Z>0,"N"]/(voxelSz*voxelSz*voxelSz)
        # account for removing points < 1 m height
        voxel_als_pai[voxel_als_pai$Z==0,"pts_m3"] <- voxel_als_pai[voxel_als_pai$Z==0,"N"]/((voxelSz-1)*voxelSz*voxelSz)
      # rescale X values
        voxel_als_pai$Xplot <- voxel_als_pai$X - 364560
        voxel_als_pai$Yplot <- voxel_als_pai$Z + 2.5
      # Use McArthur Horn to estimate effective PAI
        ePAI_als <- MH_PAI(X=voxel_als_pai$Xplot,
                           Y=1,
                           Z=voxel_als_pai$Yplot,
                           N=voxel_als_pai$pts_m3,
                           type="air")
        pai_als_vect <- vect(ePAI_als,
                               geom=c("X", "Z"))
        paiRast_als <- terra::rasterize(pai_als_vect,rastTemplate, field="ePAI")
      
### Drone
      
      # read lidar data
      data <- readLAS(droneFile)
      # subtract ground height
      dataNorm <- normalize_height(data, dtm)
      
      # 2D POINT DENSITY
      # calculate voxel point density
      voxel_drone <- voxel_metrics(dataNorm, ~list(N = length(Z)), 5, all_voxels = T)
      # replace "NA" values with 0
      voxel_drone[is.na(voxel_drone$N),"N"] <- 0
      # creates 2 "layers", recombine into one
      voxel_drone <- aggregate(N~X+Z, data = voxel_drone, sum)
      # scale to pts/m3
      voxel_drone$pts_m3 <- voxel_drone$N/25
      # rescale X values
      voxel_drone$Xplot <- voxel_drone$X - 364560
      voxel_drone$Yplot <- voxel_drone$Z + 2.5
      voxel_drone_vect <- vect(voxel_drone,
                             geom=c("Xplot", "Yplot"))
      densRast_drone <- terra::rasterize(voxel_drone_vect,rastTemplate, field="pts_m3")
      
      # CANOPY HEIGHT
      chm_drone <- grid_canopy(dataNorm, res=1,
                             algorithm = p2r(subcircle=0.01, 
                                             na.fill = tin()))
      
      # Voxelized effective PAI
      
        # only keep first returns, and remove points < 1 m height
        dataNormVeg <- dataNorm[dataNorm$Z>=2 & dataNorm$ReturnNumber==1,]
        # re-calculate voxel point density
        voxel_drone_pai <- voxel_metrics(dataNormVeg, ~list(N = length(Z)), voxelSz, all_voxels = T)
        # replace "NA" values with 0
        voxel_drone_pai[is.na(voxel_drone_pai$N),"N"] <- 0
        # creates 2 "layers", recombine into one
        voxel_drone_pai <- aggregate(N~X+Z, data = voxel_drone_pai, sum)
        # scale to pts/m3
        voxel_drone_pai[voxel_drone_pai$Z>0,"pts_m3"] <- voxel_drone_pai[voxel_drone_pai$Z>0,"N"]/(voxelSz*voxelSz*voxelSz)
        # account for removing points < 1 m height
        voxel_drone_pai[voxel_drone_pai$Z==0,"pts_m3"] <- voxel_drone_pai[voxel_drone_pai$Z==0,"N"]/((voxelSz-1)*voxelSz*voxelSz)
        # rescale X values
        voxel_drone_pai$Xplot <- voxel_drone_pai$X - 364560
        voxel_drone_pai$Yplot <- voxel_drone_pai$Z + 2.5
        # Use McArthur Horn to estimate effective PAI
        ePAI_drone <- MH_PAI(X=voxel_drone_pai$Xplot,
                           Y=1,
                           Z=voxel_drone_pai$Yplot,
                           N=voxel_drone_pai$pts_m3,
                           type="air")
        pai_drone_vect <- vect(ePAI_drone,
                             geom=c("X", "Z"))
        paiRast_drone <- terra::rasterize(pai_drone_vect,rastTemplate, field="ePAI")
      

### MLS
      
      # read lidar data
      data <- readLAS(mlsFile)
      # subtract ground height
      dataNorm <- normalize_height(data, dtm)
      
      # 2D POINT DENSITY
      # calculate voxel point density
      voxel_mls <- voxel_metrics(dataNorm, ~list(N = length(Z)), 5, all_voxels = T)
      # replace "NA" values with 0
      voxel_mls[is.na(voxel_mls$N),"N"] <- 0
      # creates 2 "layers", recombine into one
      voxel_mls <- aggregate(N~X+Z, data = voxel_mls, sum)
      # scale to pts/m3
      voxel_mls$pts_m3 <- voxel_mls$N/25
      # rescale X values
      voxel_mls$Xplot <- voxel_mls$X - 364560
      voxel_mls$Yplot <- voxel_mls$Z + 2.5
      voxel_mls_vect <- vect(voxel_mls,
                             geom=c("Xplot", "Yplot"))
      densRast_mls <- terra::rasterize(voxel_mls_vect,rastTemplate, field="pts_m3")
      
      # CANOPY HEIGHT
      chm_mls <- grid_canopy(dataNorm, res=1,
                             algorithm = p2r(subcircle=0.01, 
                                             na.fill = tin()))
      
    # Voxelized effective PAI
      
      # only keep first returns, and remove points < 1 m height
      dataNormVeg <- dataNorm[dataNorm$Z>=2 & dataNorm$ReturnNumber==1,]
      # re-calculate voxel point density
      voxel_mls_pai <- voxel_metrics(dataNormVeg, ~list(N = length(Z)), voxelSz, all_voxels = T)
      # replace "NA" values with 0
      voxel_mls_pai[is.na(voxel_mls_pai$N),"N"] <- 0
      # creates 2 "layers", recombine into one
      voxel_mls_pai <- aggregate(N~X+Z, data = voxel_mls_pai, sum)
      # scale to pts/m3
      voxel_mls_pai[voxel_mls_pai$Z>0,"pts_m3"] <- voxel_mls_pai[voxel_mls_pai$Z>0,"N"]/(voxelSz*voxelSz*voxelSz)
      # account for removing points < 1 m height
      voxel_mls_pai[voxel_mls_pai$Z==0,"pts_m3"] <- voxel_mls_pai[voxel_mls_pai$Z==0,"N"]/((voxelSz-1)*voxelSz*voxelSz)
      # rescale X values
      voxel_mls_pai$Xplot <- voxel_mls_pai$X - 364560
      voxel_mls_pai$Yplot <- voxel_mls_pai$Z + 2.5
      # Use McArthur Horn to estimate effective PAI
      ePAI_mls <- MH_PAI(X=voxel_mls_pai$Xplot,
                         Y=1,
                         Z=voxel_mls_pai$Yplot,
                         N=voxel_mls_pai$pts_m3,
                         type="ground")
      pai_mls_vect <- vect(ePAI_mls,
                           geom=c("X", "Z"))
      paiRast_mls <- terra::rasterize(pai_mls_vect,rastTemplate, field="ePAI")    
      
      
### TLS
      
      # read lidar data
      data <- readLAS(tlsFile)
      # subtract ground height
      dataNorm <- normalize_height(data, dtm)
      
      # 2D POINT DENSITY
      # calculate voxel point density
      voxel_tls <- voxel_metrics(dataNorm, ~list(N = length(Z)), 5, all_voxels = T)
      # replace "NA" values with 0
      voxel_tls[is.na(voxel_tls$N),"N"] <- 0
      # creates 2 "layers", recombine into one
      voxel_tls <- aggregate(N~X+Z, data = voxel_tls, sum)
      # scale to pts/m3
      voxel_tls$pts_m3 <- voxel_tls$N/25
      # rescale X values
      voxel_tls$Xplot <- voxel_tls$X - 364560
      voxel_tls$Yplot <- voxel_tls$Z + 2.5
      voxel_tls_vect <- vect(voxel_tls,
                             geom=c("Xplot", "Yplot"))
      densRast_tls <- terra::rasterize(voxel_tls_vect,rastTemplate, field="pts_m3")
      
      # CANOPY HEIGHT
      chm_tls <- grid_canopy(dataNorm, res=1,
                             algorithm = p2r(subcircle=0.01, 
                                             na.fill = tin()))
      
    # Voxelized effective PAI
      
      # only keep first returns, and remove points < 2 m height
      dataNormVeg <- dataNorm[dataNorm$Z>=2 & dataNorm$ReturnNumber==1,]
      # re-calculate voxel point density
      voxel_tls_pai <- voxel_metrics(dataNormVeg, ~list(N = length(Z)), voxelSz, all_voxels = T)
      # replace "NA" values with 0
      voxel_tls_pai[is.na(voxel_tls_pai$N),"N"] <- 0
      # creates 2 "layers", recombine into one
      voxel_tls_pai <- aggregate(N~X+Z, data = voxel_tls_pai, sum)
      # scale to pts/m3
      voxel_tls_pai[voxel_tls_pai$Z>0,"pts_m3"] <- voxel_tls_pai[voxel_tls_pai$Z>0,"N"]/(voxelSz*voxelSz*voxelSz)
      # account for removing points < 1 m height
      voxel_tls_pai[voxel_tls_pai$Z==0,"pts_m3"] <- voxel_tls_pai[voxel_tls_pai$Z==0,"N"]/((voxelSz-1)*voxelSz*voxelSz)
      # rescale X values
      voxel_tls_pai$Xplot <- voxel_tls_pai$X - 364560
      voxel_tls_pai$Yplot <- voxel_tls_pai$Z + 2.5
      # Use McArthur Horn to estimate effective PAI
      ePAI_tls <- MH_PAI(X=voxel_tls_pai$Xplot,
                         Y=1,
                         Z=voxel_tls_pai$Yplot,
                         N=voxel_tls_pai$pts_m3,
                         type="ground")
      pai_tls_vect <- vect(ePAI_tls,
                           geom=c("X", "Z"))
      paiRast_tls <- terra::rasterize(pai_tls_vect,rastTemplate, field="ePAI")
      
#### Point cloud plots ####

      
      
jpeg(filename = "PointCloudPlot.jpeg",
     width = 1200, height = 3000, units = "px", pointsize = 36,
     quality = 300)

  par(mfrow=c(4,1), oma=c(2,2,1,0), las=1, mar=c(3,4,1,1))
  ptCex <- 0.05      
  
    #ALS      
      data <- readLAS(alsFile)   
      dataNorm <- normalize_height(data, dtm)
      plot(x=dataNorm$X- 364560,y=dataNorm$Z,
           cex=ptCex,
           pch=19,
           col = adjustcolor("black",0.5),
           ylim=c(0,45),
           ylab=NA,
           asp=1,
           axes=F)
      mtext("Airborne lidar",side=3,line=-3,outer=F)
      axis(side=2,at=seq(0,45,5),pos=-1)
      
    
    #Drone   
    data <- readLAS(droneFile)   
    dataNorm <- normalize_height(data, dtm)
    plot(x=dataNorm$X- 364560,y=dataNorm$Z,
         cex=ptCex,
         pch=19,
         col = adjustcolor("black",0.5),
         ylim=c(0,45),
         ylab=NA,
         asp=1,
         axes=F)
    mtext("Drone lidar",side=3,line=-3,outer=F)
    axis(side=2,at=seq(0,45,5),pos=-1)
    
    #MLS   
    data <- readLAS(mlsFile)   
    dataNorm <- normalize_height(data, dtm)
    plot(x=dataNorm$X- 364560,y=dataNorm$Z,
         cex=ptCex,
         pch=19,
         col = adjustcolor("black",0.5),
         ylim=c(0,45),
         ylab=NA,
         asp=1,
         axes=F)
    mtext("Mobile lidar",side=3,line=-3,outer=F)
    axis(side=2,at=seq(0,45,5),pos=-1)
    
    #TLS   
    data <- readLAS(tlsFile)   
    dataNorm <- normalize_height(data, dtm)
    plot(x=dataNorm$X- 364560,y=dataNorm$Z,
         cex=ptCex,
         pch=19,
         col = adjustcolor("black",0.5),
         ylim=c(0,45),
         ylab=NA,
         asp=1,
         axes=F)
    mtext("Terrestrial lidar",side=3,line=-3,outer=F)
    axis(side=2,at=seq(0,45,5),pos=-1)
    axis(side=1,at=seq(0,80,5),pos=-1)
    
    mtext("Ground distance (m)",side=1,line=0,outer=T)
    mtext("Height (m)",side=2,line=0,outer=T,las=0)
    
    
dev.off()  
  
  
#### Rasterized metric plots ####
  
  densRange <- range(values(densRast_als),
                     values(densRast_drone),
                     values(densRast_mls),
                     values(densRast_tls))
  densBreaks <- c(0,1e-8,1,10,20,40,80,160,320,640,1000,10000,densRange[2])
  paiRange <- range(values(paiRast_als),
                     values(paiRast_drone),
                     values(paiRast_mls),
                     values(paiRast_tls))
  paiBreaks <- c(0,1e-8,0.25,0.5,1:5,paiRange[2])
  cexLab <- 1.2
  
jpeg(filename = "VoxelMetricsPlot.jpeg",
     width = 2400, height = 3000, units = "px", pointsize = 36,
     quality = 300)
  
  par(mfrow=c(4,2),las=1, mar=c(2,2,1,1), oma=c(4,4,10,1), xpd=T)
  
  terra::plot(densRast_als,
       breaks= densBreaks,
       col = c("white",viridisLite::plasma(length(densBreaks)-2)),
       decreasing=F,
       box=F,
       las=1,
       axes=F,
       legend=F)
  axis(side=1, at=seq(0,80,5), pos=0,
       labels=F)
  axis(side=2, at=seq(0,45,5), pos=0,
       labels=T)
  legend(x = -10, y = 60, bty="n",
         c("0", "0 - 1", "1 - 10", "10 - 20","20 - 40","40 - 80",
           "80 - 160","160 - 320","320 - 640","640 - 1,000","1,000 - 10,000", "> 10,000"),
         x.intersp = 0.2,
         fill = c("white",viridisLite::plasma(length(densBreaks)-2)),
         xpd=NA,ncol=4,
         cex=1.5)
  mtext("Point density (points/m3)",side=3,line=8,outer=F, cex = cexLab)
  mtext("Airborne lidar",side=2,line=1,outer=F, las=0)
  terra::plot(paiRast_als,
              breaks= paiBreaks,
              col = c("white",viridisLite::viridis(length(paiBreaks)-2)),
              decreasing=F,
              box=F,
              las=1,
              axes=F,
              legend=F)
  axis(side=1, at=seq(0,80,5), pos=0,
       labels=F)
  axis(side=2, at=seq(0,45,5), pos=0,
       labels=F)
  legend(x = 0, y = 60, bty="n",
         c("0", "0 - 0.25", "0.25 - 0.5", "0.5 - 1","1 - 2","2 - 3",
           "3 - 4","4 - 5",">5"),
         x.intersp = 0.2,
         fill = c("white",viridisLite::viridis(length(paiBreaks)-2)),
         xpd=NA,ncol=4,
         cex=1.5)
  mtext("Effective PAI (m2/m2)",side=3,line=8,outer=F, cex = cexLab)
  
  terra::plot(densRast_drone,
              breaks= densBreaks,
              col = c("white",viridisLite::plasma(length(densBreaks)-2)),
              decreasing=F,
              box=F,
              las=1,
              axes=F,
              legend=F)
  axis(side=1, at=seq(0,80,5), pos=0,
       labels=F)
  axis(side=2, at=seq(0,45,5), pos=0,
       labels=T)
  mtext("Drone lidar",side=2,line=1,outer=F, las=0)
  terra::plot(paiRast_drone,
              breaks= paiBreaks,
              col = c("white",viridisLite::viridis(length(paiBreaks)-2)),
              decreasing=F,
              box=F,
              las=1,
              axes=F,
              legend=F)
  axis(side=1, at=seq(0,80,5), pos=0,
       labels=F)
  axis(side=2, at=seq(0,45,5), pos=0,
       labels=F)
  
  terra::plot(densRast_mls,
              breaks= densBreaks,
              col = c("white",viridisLite::plasma(length(densBreaks)-2)),
              decreasing=F,
              box=F,
              las=1,
              axes=F,
              legend=F)
  axis(side=1, at=seq(0,80,5), pos=0,
       labels=F)
  axis(side=2, at=seq(0,45,5), pos=0,
       labels=T)
  mtext("Mobile lidar",side=2,line=1,outer=F, las=0)
  terra::plot(paiRast_mls,
              breaks= paiBreaks,
              col = c("white",viridisLite::viridis(length(paiBreaks)-2)),
              decreasing=F,
              box=F,
              las=1,
              axes=F,
              legend=F)
  axis(side=1, at=seq(0,80,5), pos=0,
       labels=F)
  axis(side=2, at=seq(0,45,5), pos=0,
       labels=F)
  
  terra::plot(densRast_tls,
              breaks= densBreaks,
              col = c("white",viridisLite::plasma(length(densBreaks)-2)),
              decreasing=F,
              box=F,
              las=1,
              axes=F,
              legend=F)
  axis(side=1, at=seq(0,80,5), pos=0,
       labels=T)
  axis(side=2, at=seq(0,45,5), pos=0,
       labels=T)
  mtext("Terrestrial lidar",side=2,line=1,outer=F, las=0)

  terra::plot(paiRast_tls,
              breaks= paiBreaks,
              col = c("white",viridisLite::viridis(length(paiBreaks)-2)),
              decreasing=F,
              box=F,
              las=1,
              axes=F,
              legend=F)
  axis(side=1, at=seq(0,80,5), pos=0,
       labels=T)
  axis(side=2, at=seq(0,45,5), pos=0,
       labels=F)
  mtext("Ground distance (m)",side=1,line=1,outer=T, cex = cexLab)
  mtext("Height (m)",side=2,line=1,outer=T,las=0, cex = cexLab)
  
dev.off()  
  