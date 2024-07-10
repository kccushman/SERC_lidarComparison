# Subset a common "slice" from all point clouds

library("lidR")

rootDir <- "/Volumes/KC_JPL/SERC_lidar/"

# Define file path for MLS, ALS, and drone lidar files
  mlsFile <- paste0(rootDir,"mls_ha4_alignedToALS.laz")
  alsFile <- paste0(rootDir,"als_ha4_clean.laz")
  droneFile <- paste0(rootDir,"drone_ha4_alignedToALS.laz")
  tlsFile <- paste0(rootDir,"TLS_alignedToALS/")
  
# Make lidR catalog objects
  mlsCat <- catalog(mlsFile)
  alsCat <- catalog(alsFile)
  droneCat <- catalog(droneFile)
  tlsCat <- catalog(tlsFile)

# Define transect ends and width (in m)
  transectP1 <- c(364560,4305790)
  transectP2 <- c(364640,4305790)
  transectWidth <- 5

# Subset point clouds
  mlsSub <-   clip_transect(mlsCat, p1 = transectP1, p2 = transectP2, width = transectWidth)
  writeLAS(mlsSub,"/Volumes/KC_JPL/SERC_lidar/transect_mls.laz")
  
  alsSub <-   clip_transect(alsCat, p1 = transectP1, p2 = transectP2, width = transectWidth)
  writeLAS(alsSub,"/Volumes/KC_JPL/SERC_lidar/transect_als.laz")
  
  droneSub <-   clip_transect(droneCat, p1 = transectP1, p2 = transectP2, width = transectWidth)
  writeLAS(droneSub,"/Volumes/KC_JPL/SERC_lidar/transect_drone.laz")
  
  tlsSub <-   clip_transect(tlsCat, p1 = transectP1, p2 = transectP2, width = transectWidth)
  writeLAS(tlsSub,"/Volumes/KC_JPL/SERC_lidar/transect_tls.laz")
  

  