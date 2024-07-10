library("lidR")
library("sf")

  # # running the whole point cloud crashed my ORNL laptop
  # 
  #   tls <- catalog("/Volumes/KC_JPL/SERC_lidar/TLS_ha4.laz")
  #   crs(tls) <- "epsg:26985"
  #   tls_data <- readLAS(tls)
  #   crs(tls_data) <- "epsg:26985"
  #   tls_utm <- st_transform(tls_data, crs("epsg:32618"))
  #   writeLAS(tls_utm,"/Volumes/KC_JPL/SERC_lidar/TLS_ha4_utm.laz")
  
  
  # # tile TLS data into smaller files
  #   opt_chunk_buffer(tls) <- 0
  #   opt_chunk_size(tls) <- 50
  #   opt_output_files(tls) <- "/Volumes/KC_JPL/SERC_lidar/TLS_tiled/TLS_tiled{XLEFT}_{YBOTTOM}"
  #   opt_laz_compression(tls) <- T
  #   catalog_retile(tls)
    
# transform TLS point cloud to match others
  # TLSfiles <- list.files("/Volumes/KC_JPL/SERC_lidar/TLS_tiled/",full.names = T)
  # 
  # for(i in 1:length(TLSfiles)){
  #   tls_data <- readLAS(TLSfiles[i])
  #   crs(tls_data) <- "epsg:26985"
  #   tls_utm <- st_transform(tls_data, crs("epsg:32618"))
  #   writeLAS(tls_utm,gsub("TLS_tiled//","TLS_tiled_utm/utm_",TLSfiles[i]))
  # }
  # 
  # TLSfiles_utm <- list.files("/Volumes/KC_JPL/SERC_lidar/TLS_tiled_utm/",full.names = T)
  # for(i in 1:length(TLSfiles_utm)){
  #   cat_utm <- catalog(TLSfiles_utm[i])
  #   plot(extent(cat_utm),add=T,col="red")
  # }
  

  # drone_ha4 <-   clip_roi(drone,extent(tls))
  # writeLAS(drone_ha4,"/Volumes/KC_JPL/SERC_lidar/drone_ha4.laz")
  # 
  # als_ha4 <-   clip_roi(als,extent(tls))
  # writeLAS(als_ha4,"/Volumes/KC_JPL/SERC_lidar/als_ha4.laz")
  # 
    mls_ha4 <- clip_roi(mls,extent(tls))
    writeLAS(mls_ha4,"/Volumes/KC_JPL/SERC_lidar/mls_ha4.laz")
    
    
drone <- catalog("/Volumes/KC_JPL/SERC_lidar/drone_ha4.laz")
  crs(drone) <- "epsg:32618"
  
als <- catalog("/Volumes/KC_JPL/SERC_lidar/als_ha4.laz")
  crs(als) <- "epsg:32618"
  
tls <- catalog("/Volumes/KC_JPL/SERC_lidar/TLS_tiled_utm/")
  crs(tls) <- "epsg:32618"
  
mls <- catalog("/Volumes/KC_JPL/SERC_lidar/mls_ha4.laz")
  crs(mls) <- "epsg:32618"

  
# make decimated versions for alignment

opt_output_files(drone) <- "/Volumes/KC_JPL/SERC_lidar/drone_ha4_align"
opt_laz_compression(drone) <- T
drone_align <- decimate_points(drone, algorithm=highest(res=0.5))

opt_output_files(als) <- "/Volumes/KC_JPL/SERC_lidar/als_ha4_align"
opt_laz_compression(als) <- T
als_align <- decimate_points(als, algorithm=highest(res=0.5))  

opt_output_files(mls) <- "/Volumes/KC_JPL/SERC_lidar/mls_ha4_align"
opt_laz_compression(mls) <- T
mls_align <- decimate_points(mls, algorithm=highest(res=0.5)) 



tlsFiles <- list.files("/Volumes/KC_JPL/SERC_lidar/TLS_tiled_utm/",full.names = T)
for(i in 1:length(tlsFiles)){
  tlsData <- readLAS(tlsFiles[i])
  tls_align <- decimate_points(tlsData, algorithm=highest(res=0.5))
  writeLAS(tls_align, paste0("/Volumes/KC_JPL/SERC_lidar/TLS_tiled_toAlign/tls_ha4_align_",i,".laz"))
}

