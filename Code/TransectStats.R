#### packages ####
library("lidR")
library("terra")

#### setup ####

  # Define file path for MLS, ALS, and drone lidar files
  mlsFile <- "Data/transect/transect_mls.laz"
  alsFile <- "Data/transect/transect_als.laz"
  droneFile <- "Data/transect/transect_drone.laz"
  droneFile_leafOff <- "Data/transect/transect_drone_leafoff.laz"
  tlsFile <- "Data/transect/tls/"
  
  # Make lidR catalog objects
  mlsCat <- catalog(mlsFile)
  alsCat <- catalog(alsFile)
  droneCat <- catalog(droneFile)
  droneCat_leafOff <- catalog(droneFile_leafOff)
  tlsCat <- catalog(tlsFile)

#### make terrain model ####

  # use previously classified ALS data for entire ha 4
  alsAll <- catalog("Data/ha4_data/als_ha4.laz")
  dtm <- rasterize_terrain(alsAll, res = 1, algorithm = knnidw(k = 6L, p = 2))

#### Calculate summary stats of transects: setup ####
  
#  Make raster template for voxelized point count/density
  
  # define voxel size
  voxelSz <- 5
  
  # create template
  rastTemplate <- rast(nrows=45/voxelSz, ncols = 80/voxelSz,
                       xmin=0,xmax=80,
                       ymin=0,ymax=45)
  
  # Define function to calculate PAI using MacArthur Horn method
  # Input data:
      # X: center coordinate of voxels in x-direction
      # Y: center coordinate of voxels in y-direction
      # Z: center coordinate of voxels in z-direction (vertical direction)
      # N: number coordinate of returns in voxels
  
  MH_PAI <- function(X,Y,Z,N,endN,type="air"){
    
    # combine into data frame
    dataPoints <- data.frame(X,Y,Z,N,endN) # make a data frame of input data
    dataPoints$XY <- paste(X,Y,sep="-") # make a unique identifier for each "column" of voxels
    dataPoints$ePAI <- NA # column to store ePAI
    
    # get unique X,Y columns of voxels
    uniqueXY <- unique(dataPoints$XY)
    
    for(i in 1:length(uniqueXY)){
      
      # get data for one "column" of voxels at a time
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
        
        # for the last voxel, use the adjust number of pulses in/out of the end of the column
        if(!is.na(data_i$endN[j])){
          pulseIn <- data_i$N[j] + data_i$endN[j]
          pulseOut <- data_i$endN[j]
        }
        
        dataPoints[dataPoints$XY==uniqueXY[i] & dataPoints$Z==data_i$Z[j],"ePAI"] <- log(pulseIn/pulseOut)
        
        if(pulseIn==0|pulseOut==0){
          dataPoints[dataPoints$XY==uniqueXY[i] & dataPoints$Z==data_i$Z[j],"ePAI"] <- 0
        }
      }
      
    }
    
    return(dataPoints)
  }
  
#### Calculate summary stats of transects: ALS ####
  
  # Read lidar data
    data <- readLAS(alsFile)
    # subtract ground height
    dataNorm <- normalize_height(data, dtm)

  # 2D POINT DENSITY
    # create a data frame to store results
    voxel_als <- data.frame(X = rep(seq(364560,(364640-5),5),each=length(seq(0,40,5))), # leftmost x-value for each voxel
                                Y = 1, # dummy y variable for this example because constant 5 m y interval for this transect
                                Z = rep(seq(0,40,5),length(seq(364560,(364640-5),5))),
                                N = NA) # number of returns to be calculated
    # Calculate voxel point density
    for(i in 1:nrow(voxel_als)){
      voxel_als$N[i] <- length(dataNorm$Z[dataNorm$X>=voxel_als$X[i] & dataNorm$X<(voxel_als$X[i]+voxelSz)
                                                 & dataNorm$Z>=voxel_als$Z[i] & dataNorm$Z<(voxel_als$Z[i]+voxelSz)])
    }
    # scale to pts/m3
    voxel_als$pts_m3 <- voxel_als$N/(voxelSz^3)
    # rescale X values
    voxel_als$Xplot <- voxel_als$X - 364560
    voxel_als$Yplot <- voxel_als$Z + 2.5
    voxel_als_vect <- vect(voxel_als,
                           geom=c("Xplot", "Yplot"))
    densRast_als <- terra::rasterize(voxel_als_vect,rastTemplate, field="pts_m3")
    
  # CANOPY HEIGHT
      chm_als <- rasterize_canopy(dataNorm, res=0.25,
                             algorithm = p2r(subcircle=0.01))
      
  # VOXELIZED EFFECTIVE PAI
      
      # create a data frame to store results
      voxel_als_pai <- data.frame(X = rep(seq(364560,(364640-5),5),each=length(seq(0,40,5))), # leftmost x-value for each voxel
                                  Y = 1, # dummy y variable for this example because constant 5 m y interval for this transect
                                  Z = rep(seq(0,40,5),length(seq(364560,(364640-5),5))),
                                  N = NA, # number of returns to be calculated
                                  nOut = NA)
      
      # Only keep first returns, and remove points < 1 m height
        htMin <- 1
        dataNormVeg <- dataNorm[dataNorm$Z>= htMin & dataNorm$ReturnNumber==1,]
        
        # Calculate voxel point density
        for(i in 1:nrow(voxel_als_pai)){
          voxel_als_pai$N[i] <- length(dataNormVeg$Z[dataNormVeg$X>=voxel_als_pai$X[i] & dataNormVeg$X<(voxel_als_pai$X[i]+voxelSz)
                                                     & dataNormVeg$Z>=voxel_als_pai$Z[i] & dataNormVeg$Z<(voxel_als_pai$Z[i]+voxelSz)])
        }
      
      # For each column, calculate the number of points below height threshold to use as the true number of beams "out" of the voxel
        voxel_als_pai$nOut <- NA
        voxel_als_pai$voxelOutHt <- NA
        xValues <- unique(voxel_als_pai$X)
        
        for(i in xValues){
          nOut_i <- nrow(dataNorm[dataNorm$Z < htMin & dataNorm$ReturnNumber==1 & dataNorm$X>=(i) & dataNorm$X < (i+voxelSz),])
          
          # if there is at least one point below the threshold, use the number of points as the number of lasers out of the voxel and adjust voxel height based on the minimim threshold
          if(nOut_i>0){
            voxel_als_pai[voxel_als_pai$X==i & voxel_als_pai$Z==0,"nOut"] <- nOut_i
            voxel_als_pai[voxel_als_pai$X==i & voxel_als_pai$Z==0,"voxelOutHt"] <- voxelSz - htMin
          }
          
          # if no points below the threshold, find the height of the lowest point
          if(nOut_i==0){
            minZ <- min(dataNormVeg$Z[dataNormVeg$X>=i & dataNormVeg$X < (i+voxelSz)])
            
          # find correct row to edit, assign the number out = 1, and change the voxel height to reflect new z distance to shortest point
            whichRow <- which(voxel_als_pai$X==i  & minZ>=voxel_als_pai$Z & minZ<(voxel_als_pai$Z+5))
            voxel_als_pai[whichRow,"nOut"] <- 1
            voxel_als_pai[whichRow,"voxelOutHt"] <- voxel_als_pai$Z[whichRow] + 5 - minZ
          }
          
        }
        
      # scale to pts/m3
      voxel_als_pai[,"pts_m3"] <- voxel_als_pai[,"N"]/(voxelSz*voxelSz*voxelSz)
        # account for edits to last return in each column
        voxel_als_pai[!is.na(voxel_als_pai$nOut),"pts_m3"] <- voxel_als_pai[!is.na(voxel_als_pai$nOut),"N"]/(voxel_als_pai[!is.na(voxel_als_pai$nOut),"voxelOutHt"]*voxelSz*voxelSz)
        voxel_als_pai[!is.na(voxel_als_pai$nOut),"nOut_m3"]<-  voxel_als_pai[!is.na(voxel_als_pai$nOut),"nOut"]/(voxel_als_pai[!is.na(voxel_als_pai$nOut),"voxelOutHt"]*voxelSz*voxelSz)
        
      # rescale X values
        voxel_als_pai$Xplot <- voxel_als_pai$X - 364560 + 2.5
        voxel_als_pai$Yplot <- voxel_als_pai$Z + 2.5
      
        # Use McArthur Horn to estimate effective PAI
        ePAI_als <- MH_PAI(X=voxel_als_pai$Xplot,
                           Y=1, # dummy y value because we have a transect that is exactly 5 m in the y direction
                           Z=voxel_als_pai$Yplot,
                           N=voxel_als_pai$pts_m3,
                           endN=voxel_als_pai$nOut_m3,
                           type="air")
        pai_als_vect <- vect(ePAI_als,
                               geom=c("X", "Z"))
        paiRast_als <- terra::rasterize(pai_als_vect,rastTemplate, field="ePAI")
      
#### Calculate summary stats of transects: ULS leaf on ####
        
  # Read lidar data
    data <- readLAS(droneFile)
    # subtract ground height
    dataNorm <- normalize_height(data, dtm)
      
  # 2D POINT DENSITY
    # create a data frame to store results
    voxel_drone <- data.frame(X = rep(seq(364560,(364640-5),5),each=length(seq(0,40,5))), # leftmost x-value for each voxel
                            Y = 1, # dummy y variable for this example because constant 5 m y interval for this transect
                            Z = rep(seq(0,40,5),length(seq(364560,(364640-5),5))),
                            N = NA) # number of returns to be calculated
    # Calculate voxel point density
    for(i in 1:nrow(voxel_drone)){
      voxel_drone$N[i] <- length(dataNorm$Z[dataNorm$X>=voxel_drone$X[i] & dataNorm$X<(voxel_drone$X[i]+voxelSz)
                                          & dataNorm$Z>=voxel_drone$Z[i] & dataNorm$Z<(voxel_drone$Z[i]+voxelSz)])
    }
      
      # scale to pts/m3
      voxel_drone$pts_m3 <- voxel_drone$N/25
      # rescale X values
      voxel_drone$Xplot <- voxel_drone$X - 364560
      voxel_drone$Yplot <- voxel_drone$Z + 2.5
      voxel_drone_vect <- vect(voxel_drone,
                             geom=c("Xplot", "Yplot"))
      densRast_drone <- terra::rasterize(voxel_drone_vect,rastTemplate, field="pts_m3")
      
      # CANOPY HEIGHT
      chm_drone <- rasterize_canopy(dataNorm, res=0.25,
                             algorithm = p2r(subcircle=0.01))
      
      # Voxelized effective PAI
      
      # create a data frame to store results
      voxel_drone_pai <- data.frame(X = rep(seq(364560,(364640-5),5),each=length(seq(0,40,5))), # leftmost x-value for each voxel
                                  Y = 1, # dummy y variable for this example because constant 5 m y interval for this transect
                                  Z = rep(seq(0,40,5),length(seq(364560,(364640-5),5))),
                                  N = NA, # number of returns to be calculated
                                  nOut = NA)
      
      # Only keep first returns, and remove points < 1 m height
      htMin <- 1
      dataNormVeg <- dataNorm[dataNorm$Z>= htMin & dataNorm$ReturnNumber==1,]
      
      # Calculate voxel point density
      for(i in 1:nrow(voxel_drone_pai)){
        voxel_drone_pai$N[i] <- length(dataNormVeg$Z[dataNormVeg$X>=voxel_drone_pai$X[i] & dataNormVeg$X<(voxel_drone_pai$X[i]+voxelSz)
                                                   & dataNormVeg$Z>=voxel_drone_pai$Z[i] & dataNormVeg$Z<(voxel_drone_pai$Z[i]+voxelSz)])
      }
      
      # For each column, calculate the number of points below height threshold to use as the true number of beams "out" of the voxel
      voxel_drone_pai$nOut <- NA
      voxel_drone_pai$voxelOutHt <- NA
      xValues <- unique(voxel_drone_pai$X)
      
      for(i in xValues){
        nOut_i <- nrow(dataNorm[dataNorm$Z < htMin & dataNorm$ReturnNumber==1 & dataNorm$X>=(i) & dataNorm$X < (i+voxelSz),])
        
        # if there is at least one point below the threshold, use the number of points as the number of lasers out of the voxel and adjust voxel height based on the minimim threshold
        if(nOut_i>0){
          voxel_drone_pai[voxel_drone_pai$X==i & voxel_drone_pai$Z==0,"nOut"] <- nOut_i
          voxel_drone_pai[voxel_drone_pai$X==i & voxel_drone_pai$Z==0,"voxelOutHt"] <- voxelSz - htMin
        }
        
          # find correct row to edit, assign the number out = 1, and change the voxel height to reflect new z distance to shortest point
        if(nOut_i==0){
          minZ <- min(dataNormVeg$Z[dataNormVeg$X>=i & dataNormVeg$X < (i+voxelSz)])
          
          # find correct row to edit
          whichRow <- which(voxel_drone_pai$X==i  & minZ>=voxel_drone_pai$Z & minZ<(voxel_drone_pai$Z+5))
          voxel_drone_pai[whichRow,"nOut"] <- 1
          voxel_drone_pai[whichRow,"voxelOutHt"] <- voxel_drone_pai$Z[whichRow] + 5 - minZ
        }
        
      }
      
      # scale to pts/m3
      voxel_drone_pai[,"pts_m3"] <- voxel_drone_pai[,"N"]/(voxelSz*voxelSz*voxelSz)
      # account for edits to last return in each column
      voxel_drone_pai[!is.na(voxel_drone_pai$nOut),"pts_m3"] <- voxel_drone_pai[!is.na(voxel_drone_pai$nOut),"N"]/(voxel_drone_pai[!is.na(voxel_drone_pai$nOut),"voxelOutHt"]*voxelSz*voxelSz)
      voxel_drone_pai[!is.na(voxel_drone_pai$nOut),"nOut_m3"]<-  voxel_drone_pai[!is.na(voxel_drone_pai$nOut),"nOut"]/(voxel_drone_pai[!is.na(voxel_drone_pai$nOut),"voxelOutHt"]*voxelSz*voxelSz)
      
      # rescale X values
      voxel_drone_pai$Xplot <- voxel_drone_pai$X - 364560 + 2.5
      voxel_drone_pai$Yplot <- voxel_drone_pai$Z + 2.5
      
      # Use McArthur Horn to estimate effective PAI
      ePAI_drone <- MH_PAI(X=voxel_drone_pai$Xplot,
                         Y=1, # dummy y value because we have a transect that is exactly 5 m in the y direction
                         Z=voxel_drone_pai$Yplot,
                         N=voxel_drone_pai$pts_m3,
                         endN=voxel_drone_pai$nOut_m3,
                         type="air")
      pai_drone_vect <- vect(ePAI_drone,
                           geom=c("X", "Z"))
      paiRast_drone <- terra::rasterize(pai_drone_vect,rastTemplate, field="ePAI")
      
      
#### Calculate summary stats of transects: ULS leaf off ####
      
        # read lidar data
        data <- readLAS(droneFile_leafOff)
        # subtract ground height
        dataNorm <- normalize_height(data, dtm)
        
    # 2D POINT DENSITY
      # create a data frame to store results
      voxel_droneLO <- data.frame(X = rep(seq(364560,(364640-5),5),each=length(seq(0,40,5))), # leftmost x-value for each voxel
                              Y = 1, # dummy y variable for this example because constant 5 m y interval for this transect
                              Z = rep(seq(0,40,5),length(seq(364560,(364640-5),5))),
                              N = NA) # number of returns to be calculated
      # Calculate voxel point density
      for(i in 1:nrow(voxel_droneLO)){
        voxel_droneLO$N[i] <- length(dataNorm$Z[dataNorm$X>=voxel_droneLO$X[i] & dataNorm$X<(voxel_droneLO$X[i]+voxelSz)
                                            & dataNorm$Z>=voxel_droneLO$Z[i] & dataNorm$Z<(voxel_droneLO$Z[i]+voxelSz)])
      }
      
      # scale to pts/m3
      voxel_droneLO$pts_m3 <- voxel_droneLO$N/25
      # rescale X values
      voxel_droneLO$Xplot <- voxel_droneLO$X - 364560
      voxel_droneLO$Yplot <- voxel_droneLO$Z + 2.5
      voxel_droneLO_vect <- vect(voxel_droneLO,
                               geom=c("Xplot", "Yplot"))
      densRast_droneLO <- terra::rasterize(voxel_droneLO_vect,rastTemplate, field="pts_m3")
        
    # CANOPY HEIGHT
        chm_droneLO <- rasterize_canopy(dataNorm, res=0.25,
                                      algorithm = p2r(subcircle=0.01))
        
  # Don't calculate ePAI for leaf-off data
        
#### Calculate summary stats of transects: MLS ####
        
    # Read lidar data
      data <- readLAS(mlsFile)
      # subtract ground height
      dataNorm <- normalize_height(data, dtm)
      
    # 2D POINT DENSITY
      # create a data frame to store results
      voxel_mls <- data.frame(X = rep(seq(364560,(364640-5),5),each=length(seq(0,40,5))), # leftmost x-value for each voxel
                              Y = 1, # dummy y variable for this example because constant 5 m y interval for this transect
                              Z = rep(seq(0,40,5),length(seq(364560,(364640-5),5))),
                              N = NA) # number of returns to be calculated
      # Calculate voxel point density
      for(i in 1:nrow(voxel_mls)){
        voxel_mls$N[i] <- length(dataNorm$Z[dataNorm$X>=voxel_mls$X[i] & dataNorm$X<(voxel_mls$X[i]+voxelSz)
                                            & dataNorm$Z>=voxel_mls$Z[i] & dataNorm$Z<(voxel_mls$Z[i]+voxelSz)])
      }
      # scale to pts/m3
      voxel_mls$pts_m3 <- voxel_mls$N/25
      # rescale X values
      voxel_mls$Xplot <- voxel_mls$X - 364560
      voxel_mls$Yplot <- voxel_mls$Z + 2.5
      voxel_mls_vect <- vect(voxel_mls,
                             geom=c("Xplot", "Yplot"))
      densRast_mls <- terra::rasterize(voxel_mls_vect,rastTemplate, field="pts_m3")
      
      # CANOPY HEIGHT
      chm_mls <- rasterize_canopy(dataNorm, res=0.25,
                             algorithm = p2r(subcircle=0.01))
      
  # VOXELIZED EFFECTIVE PAI
      
      # create a data frame to store results
      voxel_mls_pai <- data.frame(X = rep(seq(364560,(364640-5),5),each=length(seq(0,40,5))), # leftmost x-value for each voxel
                                  Y = 1, # dummy y variable for this example because constant 5 m y interval for this transect
                                  Z = rep(seq(0,40,5),length(seq(364560,(364640-5),5))),
                                  N = NA, # number of returns to be calculated
                                  nOut = NA)
      
      # Only keep first returns, and remove points < 1 m height
        htMin <- 1
        dataNormVeg <- dataNorm[dataNorm$Z>= htMin & dataNorm$ReturnNumber==1,]
        
        # Calculate voxel point density
        for(i in 1:nrow(voxel_mls_pai)){
          voxel_mls_pai$N[i] <- length(dataNormVeg$Z[dataNormVeg$X>=voxel_mls_pai$X[i] & dataNormVeg$X<(voxel_mls_pai$X[i]+voxelSz)
                                                     & dataNormVeg$Z>=voxel_mls_pai$Z[i] & dataNormVeg$Z<(voxel_mls_pai$Z[i]+voxelSz)])
        }
      
      # For each column, calculate the highest point and use it to correct for the number of outgoing pulses in the tallest voxel
        voxel_mls_pai$nOut <- NA
        voxel_mls_pai$voxelOutHt <- NA
        xValues <- unique(voxel_mls_pai$X)
        
        for(i in xValues){
        
          # if no points below the threshold, find the height of the tallest point
            maxZ <- max(dataNormVeg$Z[dataNormVeg$X>=i & dataNormVeg$X < (i+voxelSz)])
            
          # find correct row to edit, assign the number out = 1, and change the voxel height to reflect new z distance to tallest point
            whichRow <- which(voxel_mls_pai$X==i  & maxZ>=voxel_mls_pai$Z & maxZ<(voxel_mls_pai$Z+5))
            voxel_mls_pai[whichRow,"nOut"] <- 1
            voxel_mls_pai[whichRow,"voxelOutHt"] <- maxZ - voxel_mls_pai$Z[whichRow]

        }
        
      # scale to pts/m3
      voxel_mls_pai[,"pts_m3"] <- voxel_mls_pai[,"N"]/(voxelSz*voxelSz*voxelSz)
        # account for edits to last return in each column
        voxel_mls_pai[voxel_mls_pai$Z==0,"pts_m3"] <- voxel_mls_pai[voxel_mls_pai$Z==0,"N"]/((voxelSz-htMin)*voxelSz*voxelSz)
        voxel_mls_pai[!is.na(voxel_mls_pai$nOut),"pts_m3"] <- voxel_mls_pai[!is.na(voxel_mls_pai$nOut),"N"]/(voxel_mls_pai[!is.na(voxel_mls_pai$nOut),"voxelOutHt"]*voxelSz*voxelSz)
        voxel_mls_pai[!is.na(voxel_mls_pai$nOut),"nOut_m3"]<-  voxel_mls_pai[!is.na(voxel_mls_pai$nOut),"nOut"]/(voxel_mls_pai[!is.na(voxel_mls_pai$nOut),"voxelOutHt"]*voxelSz*voxelSz)
        
      # rescale X values
        voxel_mls_pai$Xplot <- voxel_mls_pai$X - 364560 + 2.5
        voxel_mls_pai$Yplot <- voxel_mls_pai$Z + 2.5
      
        # Use McArthur Horn to estimate effective PAI
        ePAI_mls <- MH_PAI(X=voxel_mls_pai$Xplot,
                           Y=1, # dummy y value because we have a transect that is exactly 5 m in the y direction
                           Z=voxel_mls_pai$Yplot,
                           N=voxel_mls_pai$pts_m3,
                           endN=voxel_mls_pai$nOut_m3,
                           type="ground")
        pai_mls_vect <- vect(ePAI_mls,
                               geom=c("X", "Z"))
        paiRast_mls <- terra::rasterize(pai_mls_vect,rastTemplate, field="ePAI")
      
#### Calculate summary stats of transects: TLS ####
      
    # Read lidar data
      data <- readLAS(tlsCat)
      # subtract ground height
      dataNorm <- normalize_height(data, dtm)
      
    # 2D POINT DENSITY
      # create a data frame to store results
      voxel_tls <- data.frame(X = rep(seq(364560,(364640-5),5),each=length(seq(0,40,5))), # leftmost x-value for each voxel
                              Y = 1, # dummy y variable for this example because constant 5 m y interval for this transect
                              Z = rep(seq(0,40,5),length(seq(364560,(364640-5),5))),
                              N = NA) # number of returns to be calculated
      # Calculate voxel point density
      for(i in 1:nrow(voxel_tls)){
        voxel_tls$N[i] <- length(dataNorm$Z[dataNorm$X>=voxel_tls$X[i] & dataNorm$X<(voxel_tls$X[i]+voxelSz)
                                            & dataNorm$Z>=voxel_tls$Z[i] & dataNorm$Z<(voxel_tls$Z[i]+voxelSz)])
      }
      # scale to pts/m3
      voxel_tls$pts_m3 <- voxel_tls$N/25
      # rescale X values
      voxel_tls$Xplot <- voxel_tls$X - 364560
      voxel_tls$Yplot <- voxel_tls$Z + 2.5
      voxel_tls_vect <- vect(voxel_tls,
                             geom=c("Xplot", "Yplot"))
      densRast_tls <- terra::rasterize(voxel_tls_vect,rastTemplate, field="pts_m3")
      
    # CANOPY HEIGHT
      chm_tls <- rasterize_canopy(dataNorm, res=0.25,
                             algorithm = p2r(subcircle=0.01))
      
  # VOXELIZED EFFECTIVE PAI
      
      # create a data frame to store results
      voxel_tls_pai <- data.frame(X = rep(seq(364560,(364640-5),5),each=length(seq(0,40,5))), # leftmost x-value for each voxel
                                  Y = 1, # dummy y variable for this example because constant 5 m y interval for this transect
                                  Z = rep(seq(0,40,5),length(seq(364560,(364640-5),5))),
                                  N = NA, # number of returns to be calculated
                                  nOut = NA)
      
      # Only keep first returns, and remove points < 1 m height
        htMin <- 1
        dataNormVeg <- dataNorm[dataNorm$Z>= htMin & dataNorm$ReturnNumber==1,]
        
        # Calculate voxel point density
        for(i in 1:nrow(voxel_tls_pai)){
          voxel_tls_pai$N[i] <- length(dataNormVeg$Z[dataNormVeg$X>=voxel_tls_pai$X[i] & dataNormVeg$X<(voxel_tls_pai$X[i]+voxelSz)
                                                     & dataNormVeg$Z>=voxel_tls_pai$Z[i] & dataNormVeg$Z<(voxel_tls_pai$Z[i]+voxelSz)])
        }
      
      # For each column, calculate the highest point and use it to correct for the number of outgoing pulses in the tallest voxel
        voxel_tls_pai$nOut <- NA
        voxel_tls_pai$voxelOutHt <- NA
        xValues <- unique(voxel_tls_pai$X)
        
        for(i in xValues){
        
          # if no points below the threshold, find the height of the tallest point
            maxZ <- max(dataNormVeg$Z[dataNormVeg$X>=i & dataNormVeg$X < (i+voxelSz)])
            
          # find correct row to edit, assign the number out = 1, and change the voxel height to reflect new z distance to tallest point
            whichRow <- which(voxel_tls_pai$X==i  & maxZ>=voxel_tls_pai$Z & maxZ<(voxel_tls_pai$Z+5))
            voxel_tls_pai[whichRow,"nOut"] <- 1
            voxel_tls_pai[whichRow,"voxelOutHt"] <- maxZ - voxel_tls_pai$Z[whichRow]

        }
        
      # scale to pts/m3
      voxel_tls_pai[,"pts_m3"] <- voxel_tls_pai[,"N"]/(voxelSz*voxelSz*voxelSz)
        # account for edits to last return in each column
        voxel_tls_pai[voxel_tls_pai$Z==0,"pts_m3"] <- voxel_tls_pai[voxel_tls_pai$Z==0,"N"]/((voxelSz-htMin)*voxelSz*voxelSz)
        voxel_tls_pai[!is.na(voxel_tls_pai$nOut),"pts_m3"] <- voxel_tls_pai[!is.na(voxel_tls_pai$nOut),"N"]/(voxel_tls_pai[!is.na(voxel_tls_pai$nOut),"voxelOutHt"]*voxelSz*voxelSz)
        voxel_tls_pai[!is.na(voxel_tls_pai$nOut),"nOut_m3"]<-  voxel_tls_pai[!is.na(voxel_tls_pai$nOut),"nOut"]/(voxel_tls_pai[!is.na(voxel_tls_pai$nOut),"voxelOutHt"]*voxelSz*voxelSz)
        
      # rescale X values
        voxel_tls_pai$Xplot <- voxel_tls_pai$X - 364560 + 2.5
        voxel_tls_pai$Yplot <- voxel_tls_pai$Z + 2.5
      
        # Use McArthur Horn to estimate effective PAI
        ePAI_tls <- MH_PAI(X=voxel_tls_pai$Xplot,
                           Y=1, # dummy y value because we have a transect that is exactly 5 m in the y direction
                           Z=voxel_tls_pai$Yplot,
                           N=voxel_tls_pai$pts_m3,
                           endN=voxel_tls_pai$nOut_m3,
                           type="ground")
        pai_tls_vect <- vect(ePAI_tls,
                               geom=c("X", "Z"))
        paiRast_tls <- terra::rasterize(pai_tls_vect,rastTemplate, field="ePAI")
      
#### Make table of summary stats ####
      
  transectSummary <- data.frame(instrument = c("Airplane lidar",
                                               "Drone lidar (leaf on)",
                                               "Drone lidar (leaf off)",
                                               "Mobile lidar",
                                               "Terrestrial lidar"),
                                pointDensity = NA,
                                meanCanopyHeight = NA,
                                canopyRugosity = NA, # calculated a la https://www.sciencedirect.com/science/article/abs/pii/S0378112713001254?via%3Dihub
                                meanLAI = NA)

  # Airplane
      transectSummary$pointDensity[1] <- round(alsCat$Number.of.point.records/(80*5))
      transectSummary$meanCanopyHeight[1] <- round(mean(values(chm_als),na.rm=T),1)
      transectSummary$canopyRugosity[1] <- round(sd(aggregate(ePAI_als$ePAI, by=list(ePAI_als$X), FUN="sd")[,2]),2)
      transectSummary$meanLAI[1] <- round(mean(aggregate(ePAI_als$ePAI, by=list(ePAI_als$X), FUN="sum")[,2]),2)
     
  # Drone - leaf on
      transectSummary$pointDensity[2] <- round(droneCat$Number.of.point.records/(80*5))
      transectSummary$meanCanopyHeight[2] <- round(mean(values(chm_drone),na.rm=T),1)
      transectSummary$canopyRugosity[2] <- round(sd(aggregate(ePAI_drone$ePAI, by=list(ePAI_drone$X), FUN="sd")[,2]),2)
      transectSummary$meanLAI[2] <- round(mean(aggregate(ePAI_drone$ePAI, by=list(ePAI_drone$X), FUN="sum")[,2]),2)
    
  # Drone - leaf off
      transectSummary$pointDensity[3] <- round(droneCat_leafOff$Number.of.point.records/(80*5))
      transectSummary$meanCanopyHeight[3] <- round(mean(values(chm_droneLO),na.rm=T),1)
      transectSummary$canopyRugosity[3] <- "--"
      transectSummary$meanLAI[3] <- "--" 
  # Mobile
      transectSummary$pointDensity[4] <- round(mlsCat$Number.of.point.records/(80*5))
      transectSummary$meanCanopyHeight[4] <- round(mean(values(chm_mls),na.rm=T),1)
      transectSummary$canopyRugosity[4] <- round(sd(aggregate(ePAI_mls$ePAI, by=list(ePAI_mls$X), FUN="sd")[,2]),2)
      transectSummary$meanLAI[4] <- round(mean(aggregate(ePAI_mls$ePAI, by=list(ePAI_mls$X), FUN="sum")[,2]),2)
      
  # TLS
      transectSummary$pointDensity[5] <- round(sum(tlsCat$Number.of.point.records)/(80*5))
      transectSummary$meanCanopyHeight[5] <- round(mean(values(chm_tls),na.rm=T),1)
      transectSummary$canopyRugosity[5] <- round(sd(aggregate(ePAI_tls$ePAI, by=list(ePAI_tls$X), FUN="sd")[,2]),2)
      transectSummary$meanLAI[5] <- round(mean(aggregate(ePAI_tls$ePAI, by=list(ePAI_tls$X), FUN="sum")[,2]),2)
      
      write.csv(transectSummary,"Results/transectSummaryStats.csv",row.names = F)
      
      
      transectSummary$canopyRugosity <- as.numeric(transectSummary$canopyRugosity)
      transectSummary$meanLAI <- as.numeric(transectSummary$meanLAI)
      
     # calculate % variation in mean canopy height
     round(100*(max(transectSummary$meanCanopyHeight)-min(transectSummary$meanCanopyHeight))/min(transectSummary$meanCanopyHeight),1)
     # not inlcuding leaf-off ULS
     round(100*(max(transectSummary$meanCanopyHeight[-3])-min(transectSummary$meanCanopyHeight[-3]))/min(transectSummary$meanCanopyHeight[-3]),1)
     
     # calculate % variation in mean LAI
     round(100*(max(transectSummary$meanLAI,na.rm=T)-min(transectSummary$meanLAI,na.rm=T))/min(transectSummary$meanLAI,na.rm=T),1)
     
     # calculate % variation in canopy rugosity
     round(100*(max(transectSummary$canopyRugosity,na.rm=T)-min(transectSummary$canopyRugosity,na.rm=T))/min(transectSummary$canopyRugosity,na.rm=T),1)
     

#### Figure 2. Point cloud plot: four discrete return platforms, leaf-on ####
     
tiff(filename = "Results/Figure2.tiff",
     width = 4080, height = 2720, units = "px", pointsize = 72)

  par(mfrow=c(2,2), oma=c(4,4,1,1), las=1, mar=c(1,1,1,0))
  ptCex <- 0.1      
  
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
      mtext("a. ALS",side=3,line=-1,outer=F)
      #text("a", x = 2.5, y = 45)
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
    mtext("b. ULS",side=3,line=-1,outer=F)
    #text("b", x = 2.5, y = 45)
    
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
    mtext("c. MLS",side=3,line=-1,outer=F)
    #text("c", x = 2.5, y = 45)
    axis(side=2,at=seq(0,45,5),pos=-1)
    axis(side=1,at=seq(0,80,5),pos=-1)
    
    #TLS   
    data <- readLAS(tlsCat)   
    dataNorm <- normalize_height(data, dtm)
    plot(x=dataNorm$X- 364560,y=dataNorm$Z,
         cex=ptCex,
         pch=19,
         col = adjustcolor("black",0.5),
         ylim=c(0,45),
         ylab=NA,
         asp=1,
         axes=F)
    mtext("d. TLS",side=3,line=-1,outer=F)
    #text("d", x = 2.5, y = 45)
    axis(side=1,at=seq(0,80,5),pos=-1)
    
    mtext("Ground distance (m)",side=1,line=1,outer=T)
    mtext("Height (m)",side=2,line=1,outer=T,las=0)
    
dev.off()  
  
  
#### Figure 4. Rasterized metric plots ####

densRange <- range(values(densRast_als),
                   values(densRast_drone),
                   values(densRast_mls),
                   values(densRast_tls))
densBreaks <- c(0,1e-8,1,10,20,40,80,160,320,640,1000,10000,densRange[2])
paiRange <- range(values(paiRast_als),
                  values(paiRast_drone),
                  values(paiRast_mls),
                  values(paiRast_tls))
paiBreaks <- c(0,1e-8,0.5,1,2,4,6,8,10,paiRange[2])
cexLab <- 0.9
cexLetter <- 1.2
cexAxis <- 1.2
rastMar <- c(1, 1, 1, 1)


pdf(file = "Results/Figure4.pdf",
     width = 4.33, height = 6, pointsize = 12)

par(mfrow=c(4,2),las=1, mar=c(0,0,0,0), oma=c(4,4,6,0), xpd=T)

terra::plot(densRast_als,
            breaks= densBreaks,
            col = c("white",viridisLite::plasma(length(densBreaks)-2)),
            mar = rastMar,
            decreasing=F,
            box=F,
            las=1,
            axes=F,
            legend=F)
axis(side=1, at=seq(0,80,5), pos=0,
     labels=F)
axis(side=2, at=seq(0,45,5), pos=0,
     labels=F)
axis(side=2, at=seq(0,45,20), pos=0,
     labels=T, cex.axis = cexAxis)

legend(x = 0, y = 82, bty="n",
       c("0", "< 1", "1 - 10", "10 - 20","20 - 40","40 - 80",
         "80 - 160","160 - 320","320 - 640","640 - 1,000","1,000 - 10,000", "> 10,000"),
       x.intersp = 0.4,y.intersp = 1,text.width = 30,
       fill = c("white",viridisLite::plasma(length(densBreaks)-2)),
       pt.cex = 1.2, pt.lwd = 0.5,
       xpd=NA,ncol=2,
       cex=1)

mtext("Point density (points/m3)",side=3,line=4.5,outer=F, cex = 0.7)
mtext("ALS",side=2,line=1,outer=F, las=0, cex = cexLab)
text("a", x = 5, y = 42.5, cex = cexLetter)

terra::plot(paiRast_als,
            breaks= paiBreaks,
            col = c("white",viridisLite::viridis(length(paiBreaks)-2)),
            mar = rastMar,
            decreasing=F,
            box=F,
            las=1,
            axes=F,
            legend=F)
axis(side=1, at=seq(0,80,5), pos=0,
     labels=F)
axis(side=2, at=seq(0,45,5), pos=0,
     labels=F)
legend(x = 0, y = 82, bty="n",
       c("0", "< 0.5", "0.5 - 1", "1 - 2","2 - 4","4 - 6",
         "6 - 8","8 - 10","> 10"),
       x.intersp = 0.4,y.intersp = 1,text.width = 30,
       fill = c("white",viridisLite::viridis(length(paiBreaks)-2)),
       pt.cex = 1.2, pt.lwd = 0.5,
       xpd=NA,ncol=2,
       cex=1)

mtext("Effective PAI (m2/m2)",side=3,line=4.5,outer=F, cex = 0.7)
text("b", x = 5, y = 42.5, cex = cexLetter)

terra::plot(densRast_drone,
            breaks= densBreaks,
            col = c("white",viridisLite::plasma(length(densBreaks)-2)),
            mar = rastMar,
            decreasing=F,
            box=F,
            las=1,
            axes=F,
            legend=F)
axis(side=1, at=seq(0,80,5), pos=0,
     labels=F)
axis(side=2, at=seq(0,45,5), pos=0,
     labels=F)
axis(side=2, at=seq(0,45,20), pos=0,
     labels=T, cex.axis = cexAxis)

mtext("ULS",side=2,line=1,outer=F, las=0, cex = cexLab)
text("c", x = 5, y = 42.5, cex = cexLetter)

terra::plot(paiRast_drone,
            breaks= paiBreaks,
            col = c("white",viridisLite::viridis(length(paiBreaks)-2)),
            mar = rastMar,
            decreasing=F,
            box=F,
            las=1,
            axes=F,
            legend=F)
axis(side=1, at=seq(0,80,5), pos=0,
     labels=F)
axis(side=2, at=seq(0,45,5), pos=0,
     labels=F)
text("d", x = 5, y = 42.5, cex = cexLetter)

terra::plot(densRast_mls,
            breaks= densBreaks,
            col = c("white",viridisLite::plasma(length(densBreaks)-2)),
            mar = rastMar,
            decreasing=F,
            box=F,
            las=1,
            axes=F,
            legend=F)
axis(side=1, at=seq(0,80,5), pos=0,
     labels=F)
axis(side=2, at=seq(0,45,5), pos=0,
     labels=F)
axis(side=2, at=seq(0,45,20), pos=0,
     labels=T, cex.axis = cexAxis)
mtext("MLS",side=2,line=1,outer=F, las=0, cex = cexLab)
text("e", x = 5, y = 42.5, cex = cexLetter)

terra::plot(paiRast_mls,
            breaks= paiBreaks,
            col = c("white",viridisLite::viridis(length(paiBreaks)-2)),
            mar = rastMar,
            decreasing=F,
            box=F,
            las=1,
            axes=F,
            legend=F)
axis(side=1, at=seq(0,80,5), pos=0,
     labels=F)
axis(side=2, at=seq(0,45,5), pos=0,
     labels=F)
text("f", x = 5, y = 42.5, cex = cexLetter)

terra::plot(densRast_tls,
            breaks= densBreaks,
            col = c("white",viridisLite::plasma(length(densBreaks)-2)),
            mar = rastMar,
            decreasing=F,
            box=F,
            las=1,
            axes=F,
            legend=F)
axis(side=1, at=seq(0,80,5), pos=0,
     labels=F)
axis(side=2, at=seq(0,45,5), pos=0,
     labels=F)
axis(side=1, at=seq(0,80,20), pos=0,
     labels=T, cex.axis = cexAxis)
axis(side=2, at=seq(0,45,20), pos=0,
     labels=T, cex.axis = cexAxis)
mtext("TLS",side=2,line=1,outer=F, las=0, cex = cexLab)
text("g", x = 5, y = 42.5, cex = cexLetter)

terra::plot(paiRast_tls,
            breaks= paiBreaks,
            col = c("white",viridisLite::viridis(length(paiBreaks)-2)),
            mar = rastMar,
            decreasing=F,
            box=F,
            las=1,
            axes=F,
            legend=F)
axis(side=1, at=seq(0,80,5), pos=0,
     labels=F)
axis(side=2, at=seq(0,45,5), pos=0,
     labels=F)
axis(side=1, at=seq(0,80,20), pos=0,
     labels=T, cex.axis = cexAxis)
mtext("Ground distance (m)",side=1,line=2,outer=T, cex = cexLab)
mtext("Height (m)",side=2,line=2,outer=T,las=0, cex = cexLab)
text("h", x = 5, y = 42.5, cex = cexLetter)

dev.off()  

#### Figure 5. Point cloud plot: trunk cross section ####

# These point clouds were manually subsetted in CloudCompare and saved as
# separate .laz files for ease of plotting
trunkDrone <- readLAS("Data/trunk/trunk_drone.laz")
trunkMLS <- readLAS("Data/trunk/trunk_mls.laz")
trunkTLS <- readLAS("Data/trunk/trunk_tls.laz")


zMin <- 7.97
zMax <- 8.00
xRange <- c(364623.6,364624.7)
yRange <- c(4305790.9,4305791.6)
ptCex <- 0.4

pdf(file = "Results/Figure5.pdf",
     width = 6.81, height = 2.3, pointsize = 12)

  par(mfrow=c(1,3), mar=c(0,0,0,0),oma=c(2,2,2,2))
  
  plot(x = trunkTLS$X[trunkTLS$Z>zMin & trunkTLS$Z<zMax],
       y = trunkTLS$Y[trunkTLS$Z>zMin & trunkTLS$Z<zMax],
       pch=19,
       xlim=xRange,
       ylim=yRange,
       cex=ptCex,
       axes=F,ylab=NA,xlab=NA,
       col = adjustcolor("black",0.5),
       asp=1)
  mtext("a. TLS",side=3,line=-1,outer=F)
  #text("a", x= xRange[1], y = yRange[2])
  
  lines(x=c(364623.6,364623.8),
        y=c(4305790.88,4305790.88),
        lwd=2)
  text(" 20 cm", x=364623.6,y=4305790.95, adj=0, cex=1.4)

  plot(x = trunkMLS$X[trunkMLS$Z>zMin & trunkMLS$Z<zMax],
       y = trunkMLS$Y[trunkMLS$Z>zMin & trunkMLS$Z<zMax],
       pch=19,
       xlim=xRange,
       ylim=yRange,
       cex=ptCex,
       axes=F,ylab=NA,xlab=NA,
       col = adjustcolor("black",0.5),
       asp=1)
  mtext("b. MLS",side=3,line=-1,outer=F)
  #text("b", x= xRange[1], y = yRange[2])
  
  plot(x = trunkDrone$X[trunkDrone$Z>zMin & trunkDrone$Z<zMax],
       y = trunkDrone$Y[trunkDrone$Z>zMin & trunkDrone$Z<zMax],
       pch=19,
       xlim=xRange,
       ylim=yRange,
       cex=ptCex,
       axes=F,ylab=NA,xlab=NA,
       col = adjustcolor("black",0.5),
       asp=1)
  mtext("c. ULS (leaf off)",side=3,line=-1,outer=F)
  #text("c", x= xRange[1], y = yRange[2])
  
dev.off()


#### Figures 5 and 6 are in GEDI_analyzeData.R #### 
#### Figure S2. Point cloud plot: drone lidar, leaf-on vs leaf-off ####

jpeg(filename = "Results/FigureS2.jpeg",
     width = 1800, height = 750, units = "px", pointsize = 36,
     quality = 300)

par(mfrow=c(1,2), oma=c(2,2,1,1), las=1, mar=c(1,1,1,0))
ptCex <- 0.1      

#Drone: leaf on   
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
mtext("a. Leaf on (growing season)",side=3,line=-1,outer=F)
axis(side=2,at=seq(0,45,5),pos=-1)
axis(side=1,at=seq(0,80,5),pos=-1)

#Drone: leaf off   
data <- readLAS(droneFile_leafOff)   
dataNorm <- normalize_height(data, dtm)
plot(x=dataNorm$X- 364560,y=dataNorm$Z,
     cex=ptCex,
     pch=19,
     col = adjustcolor("black",0.5),
     ylim=c(0,45),
     ylab=NA,
     asp=1,
     axes=F)
mtext("b. Leaf off",side=3,line=-1,outer=F)
axis(side=1,at=seq(0,80,5),pos=-1)

mtext("Ground distance (m)",side=1,line=1,outer=T)
mtext("Height (m)",side=2,line=1,outer=T,las=0)
dev.off()  

