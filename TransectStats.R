library("lidR")
library("terra")

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

# make terrain model
  # use previously classified ALS data for whole area
  alsAll <- catalog(paste0(rootDir,"als_ha4_clean.laz"))
  dtm = rasterize_terrain(alsAll, res = 1, algorithm = knnidw(k = 6L, p = 2))

# Calculate voxelized point density for 5 x 5 m voxels
  
  rastTemplate <- rast(nrows=9, ncols = 16,
                       xmin=0,xmax=80,
                       ymin=0,ymax=45)
  
  
  # ALS
  
    # read lidar data
    data <- readLAS(alsFile)
    # subtract ground height
    dataNorm <- normalize_height(data, dtm)
    
    # 2D POINT DENSITY
      # calculate voxel point density
      voxel_als <- voxel_metrics(dataNorm, ~list(N = length(Z)), 5, all_voxels = T)
      # replace "NA" values with 0
      voxel_als[is.na(voxel_als$N),"N"] <- 0
      # creates 2 "layers", recombine into one
      voxel_als <- aggregate(N~X+Z, data = voxel_als, sum)
      # scale to pts/m3
      voxel_als$pts_m3 <- voxel_als$N/25
      # rescale X values
      voxel_als$Xplot <- voxel_als$X - 364560
      voxel_als$Yplot <- voxel_als$Z + 2.5
      voxel_als_vect <- vect(voxel_als,
                             geom=c("Xplot", "Yplot"))
      densRast_als <- terra::rasterize(voxel_als_vect,rastTemplate, field="pts_m3")
    
    # CANOPY HEIGHT
      chm_als <- grid_canopy(dataNorm, res=1,
                             algorithm = p2r(subcircle=0.01, 
                                             na.fill = tin()))
      
    
    # drone
      
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
      
    # MLS
      
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
      
  # TLS
      
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
      
      
#### plots ####

      
      
jpeg(filename = "TestPlot.jpeg",
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
    mtext("Mobile laser scanning",side=3,line=-3,outer=F)
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
    mtext("Height (m)",side=2,line=0,outer=T)
    
    
dev.off()  
  
  plot(densRast_als)
  
  
  
  #Drone      
  data <- readLAS(droneFile)   
  plot(x=data$X,y=data$Z,
       cex=0.1,
       pch=19,
       main = "Drone lidar")