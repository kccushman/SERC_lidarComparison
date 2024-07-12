# Subset a common "slice" from all point clouds

library("lidR")

# NOTE: before running this section, download ha 4 data from the following Google Drive link:
# https://drive.google.com/drive/folders/1zTAWz94pnOWlbFfPlgONtNq5IFj58dgW?usp=sharing
# and unzip the "TLS" folder

# Define file path for MLS, ALS, and uls lidar files
  mlsFile <- "Data/ha4/mls_ha4.laz"
  alsFile <- "Data/ha4/als_ha4.laz"
  ulsFile <- "Data/ha4/drone_ha4.laz"
  ulsOffFile <- "Data/ha4/drone_leafoff_ha4.las"
  tlsFile <- "Data/ha4/TLS/"
  
# Make lidR catalog objects
  mlsCat <- catalog(mlsFile)
  alsCat <- catalog(alsFile)
  ulsCat <- catalog(ulsFile)
  ulsOffCat <- catalog(ulsOffFile)
  tlsCat <- catalog(tlsFile)

# Define transect ends and width (in m)
  transectP1 <- c(364560,4305790)
  transectP2 <- c(364640,4305790)
  transectWidth <- 5

# Subset point clouds
  mlsSub <-   clip_transect(mlsCat, p1 = transectP1, p2 = transectP2, width = transectWidth)
  writeLAS(mlsSub,"Data/transect/transect_mls.laz")
  
  alsSub <-   clip_transect(alsCat, p1 = transectP1, p2 = transectP2, width = transectWidth)
  writeLAS(alsSub,"Data/transect/transect_als.laz")
  
  ulsSub <-   clip_transect(ulsCat, p1 = transectP1, p2 = transectP2, width = transectWidth)
  writeLAS(ulsSub,"Data/transect/transect_drone.laz")
  
  ulsOffSub <-   clip_transect(ulsOffCat, p1 = transectP1, p2 = transectP2, width = transectWidth)
  writeLAS(ulsOffSub,"Data/transect/transect_drone_leafoff.laz")
  
  tlsSub <-   clip_transect(tlsCat, p1 = transectP1, p2 = transectP2, width = transectWidth)
  writeLAS(tlsSub,"Data/transect/transect_tls.laz")
    # rewrite TLS transect into 40 m tiled data so that files are small enough for GitHub
    tlsTransectCat <- catalog("Data/transect/transect_tls.laz")
    opt_chunk_buffer(tlsTransectCat) <- 0 # set 0 m buffer
    opt_chunk_size(tlsTransectCat) <- 40 # set 40 m max tile length
    opt_laz_compression(tlsTransectCat) <- T
    opt_output_files(tlsTransectCat) <- "Data/transect/tls/transect_tls_{XLEFT}_{YBOTTOM}"
    tlsTransectCat <- catalog_retile(tlsTransectCat)
    # remove single TLS file
    file.remove("Data/transect/transect_tls.laz")
