library(rhdf5)
library(bit64)
library(sp)
library(terra)

# NOTE: before running this section, download all folders from the following Google Drive link:
# https://drive.google.com/drive/folders/11y4RBRfPJSHqWiekDEiomRX5UlqNJeF0?usp=sharing
# Keep organization in "GEDI01_B", "GEDI02_A", "GEDI02_B"

#### Get 2_A data ####
files <- unlist(list.files(path="Data/GEDI/GEDI02_A", pattern="h5", recursive=T, full.names=T))

temp_i <- vector("list", length(files))

for (i in 1:length(files)){

   year <- as.integer(substr(files[i],42,45))
  
   temp <- h5ls(files[i])
   temp <- temp[which(temp$dim != "1"),]
   temp <- temp[-which(temp$dim==""),]
   beams <- unique(unlist(lapply(strsplit(temp$group, "/"), function (x) x[2])))

   fn <- files[i] 
   temp_j <- vector("list", length(beams))

   for (j in 1:length(beams)){

      temp <- data.frame(

         beam = as.numeric(h5read(fn, paste("/",beams[j],"/beam",sep=""))),

         channel = as.numeric(h5read(fn, paste("/",beams[j],"/channel",sep=""))),
         
         lat_lowestmode = as.numeric(h5read(fn, paste("/",beams[j],"/lat_lowestmode",sep=""))),

         lon_lowestmode = as.numeric(h5read(fn, paste("/",beams[j],"/lon_lowestmode",sep=""))),
         
         elev_lowestmode = as.numeric(h5read(fn, paste("/",beams[j],"/elev_lowestmode",sep=""))),
         
         sensitivity = as.numeric(h5read(fn, paste("/",beams[j],"/sensitivity",sep=""))),

         quality_flag = as.numeric(h5read(fn, paste("/",beams[j],"/quality_flag",sep=""))),
         
         shot_number = as.character(h5read(fn, paste("/",beams[j],"/shot_number",sep=""), bit64conversion="bit64")),

         solar_azimuth = as.numeric(h5read(fn, paste("/",beams[j],"/solar_azimuth",sep=""))),

         solar_elevation = as.numeric(h5read(fn, paste("/",beams[j],"/solar_elevation",sep=""))),
      
         landsat_treecover = h5read(fn, paste("/",beams[j],"/land_cover_data/landsat_treecover",sep="")),

         landsat_water_persistence = h5read(fn, paste("/",beams[j],"/land_cover_data/landsat_water_persistence",sep="")),

         leaf_off_doy = h5read(fn, paste("/",beams[j],"/land_cover_data/leaf_off_doy",sep="")),

         leaf_off_flag = h5read(fn, paste("/",beams[j],"/land_cover_data/leaf_off_flag",sep="")),

         leaf_on_cycle = h5read(fn, paste("/",beams[j],"/land_cover_data/leaf_on_cycle",sep="")),

         leaf_on_doy = h5read(fn, paste("/",beams[j],"/land_cover_data/leaf_on_doy",sep="")),

         modis_nonvegetated = h5read(fn, paste("/",beams[j],"/land_cover_data/modis_nonvegetated",sep="")),

         modis_nonvegetated_sd = h5read(fn, paste("/",beams[j],"/land_cover_data/modis_nonvegetated_sd",sep="")),

         modis_treecover= h5read(fn, paste("/",beams[j],"/land_cover_data/modis_treecover",sep="")),

         modis_treecover_sd = h5read(fn, paste("/",beams[j],"/land_cover_data/modis_treecover_sd",sep="")),

         pft_class = h5read(fn, paste("/",beams[j],"/land_cover_data/pft_class",sep="")),
  
         region_class = h5read(fn, paste("/",beams[j],"/land_cover_data/region_class",sep=""))

      )

      rh = t(h5read(fn, paste("/",beams[j],"/rh",sep="")))
      dimnames(rh)[[2]] <- paste("rh_",seq(1,dim(rh)[2],1),sep="")
      temp_j[[j]] <- cbind(temp, rh, fn)

   }

   temp_j <- do.call(rbind, temp_j)
   if(nrow(temp_j>0)){
     temp_j$year <- year
   }
   temp_i[[i]] <- temp_j

}

temp_i <- do.call(rbind, temp_i)
write.csv(temp_i, "Data/GEDI/data_GEDI2_A.csv", row.names = F, quote=T)
write.csv(temp_i$shot_number, "Data/GEDI/data_GEDI2_A_shot_number.csv", row.names = F)

#### Get 2_B data (only pulls useful metadata, vertical PAI profiles, and vertical PAVD profiles) ####
files <- unlist(list.files(path="Data/GEDI/GEDI02_B", pattern="h5", recursive=T, full.names=T))

temp_i <- vector("list", length(files))

for (i in 1:length(files)){
  
   year <- as.integer(substr(files[i],42,45))
  
   temp <- h5ls(files[i])
   temp <- temp[which(temp$dim != "1"),]
   temp <- temp[-which(temp$dim==""),]
   beams <- unique(unlist(lapply(strsplit(temp$group, "/"), function (x) x[2])))
   
   fn <- files[i] 
   temp_j <- vector("list", length(beams))
   
   for (j in 1:length(beams)){
      
      temp <- data.frame(
         
         beam = as.numeric(h5read(fn, paste("/",beams[j],"/beam",sep=""))),
         
         channel = as.numeric(h5read(fn, paste("/",beams[j],"/channel",sep=""))),
         
         lat_lowestmode = as.numeric(h5read(fn, paste("/",beams[j],"/geolocation/lat_lowestmode",sep=""))),
         
         lon_lowestmode = as.numeric(h5read(fn, paste("/",beams[j],"/geolocation/lon_lowestmode",sep=""))),
         
         sensitivity = as.numeric(h5read(fn, paste("/",beams[j],"/sensitivity",sep=""))),
         
         shot_number = as.character(h5read(fn, paste("/",beams[j],"/shot_number",sep=""), bit64conversion="bit64")),
         
         solar_azimuth = as.numeric(h5read(fn, paste("/",beams[j],"/geolocation/solar_azimuth",sep=""))),
         
         solar_elevation = as.numeric(h5read(fn, paste("/",beams[j],"/geolocation/solar_elevation",sep=""))),
         
         landsat_treecover = h5read(fn, paste("/",beams[j],"/land_cover_data/landsat_treecover",sep="")),
         
         landsat_water_persistence = h5read(fn, paste("/",beams[j],"/land_cover_data/landsat_water_persistence",sep="")),
         
         leaf_off_doy = h5read(fn, paste("/",beams[j],"/land_cover_data/leaf_off_doy",sep="")),
         
         leaf_off_flag = h5read(fn, paste("/",beams[j],"/land_cover_data/leaf_off_flag",sep="")),
         
         leaf_on_cycle = h5read(fn, paste("/",beams[j],"/land_cover_data/leaf_on_cycle",sep="")),
         
         leaf_on_doy = h5read(fn, paste("/",beams[j],"/land_cover_data/leaf_on_doy",sep="")),
         
         modis_nonvegetated = h5read(fn, paste("/",beams[j],"/land_cover_data/modis_nonvegetated",sep="")),
         
         modis_nonvegetated_sd = h5read(fn, paste("/",beams[j],"/land_cover_data/modis_nonvegetated_sd",sep="")),
         
         modis_treecover= h5read(fn, paste("/",beams[j],"/land_cover_data/modis_treecover",sep="")),
         
         modis_treecover_sd = h5read(fn, paste("/",beams[j],"/land_cover_data/modis_treecover_sd",sep="")),
         
         pft_class = h5read(fn, paste("/",beams[j],"/land_cover_data/pft_class",sep="")),
         
         region_class = h5read(fn, paste("/",beams[j],"/land_cover_data/region_class",sep="")),
         
         pai = h5read(fn, paste("/",beams[j],"/pai",sep=""))
         
      )
      
      dz = h5read(fn, paste("/",beams[j],"/ancillary/dz",sep=""))
      
      pai = t(h5read(fn, paste("/",beams[j],"/pai_z",sep="")))
      dimnames(pai)[[2]] <- paste("pai_",seq(dz,(dim(pai)[2])*dz,dz),sep="")
      
      pavd = t(h5read(fn, paste("/",beams[j],"/pavd_z",sep="")))
      dimnames(pavd)[[2]] <- paste("pavd_",seq(dz,(dim(pavd)[2])*dz,dz),sep="")
      
      temp_j[[j]] <- cbind(temp, pai, pavd, fn)
      
   }
   
   temp_j <- do.call(rbind, temp_j)
   if(nrow(temp_j>0)){
     temp_j$year <- year
   }
   temp_i[[i]] <- temp_j
   
}

temp_i <- do.call(rbind, temp_i)
write.csv(temp_i, "Data/GEDI/data_GEDI2_B.csv", row.names = F)
write.csv(temp_i$shot_number, "Data/GEDI/data_GEDI2_B_shot_number.csv", row.names = F)

#### Get 1_B data (only pulls useful metadata) ####
files <- unlist(list.files(path="Data/GEDI//GEDI01_B", pattern="h5", recursive=T, full.names=T))

temp_i <- vector("list", length(files))

for (i in 1:length(files)){
  
  year <- as.integer(substr(files[i],42,45))
  
  temp <- h5ls(files[i])
  temp <- temp[which(temp$dim != "1"),]
  temp <- temp[-which(temp$dim==""),]
  beams <- unique(unlist(lapply(strsplit(temp$group, "/"), function (x) x[2])))
  
  fn <- files[i] 
  temp_j <- vector("list", length(beams))
  
  for (j in 1:length(beams)){
    
    temp <- data.frame(
      
      beam = as.numeric(h5read(fn, paste("/",beams[j],"/beam",sep=""))),
      
      channel = as.numeric(h5read(fn, paste("/",beams[j],"/channel",sep=""))),
      
      latitude_lastbin = as.numeric(h5read(fn, paste("/",beams[j],"/geolocation/latitude_lastbin",sep=""))),
      
      longitude_lastbin = as.numeric(h5read(fn, paste("/",beams[j],"/geolocation/longitude_lastbin",sep=""))),
      
      digital_elevation_model = as.numeric(h5read(fn, paste("/",beams[j],"/geolocation/digital_elevation_model",sep=""))),
      
      shot_number = as.character(h5read(fn, paste("/",beams[j],"/shot_number",sep=""), bit64conversion="bit64"))

    )
    
    temp_j[[j]] <- cbind(temp,fn)
    
  }
  
  temp_j <- do.call(rbind, temp_j)
  if(nrow(temp_j>0)){
    temp_j$year <- year
  }
  temp_i[[i]] <- temp_j
  
}

temp_i <- do.call(rbind, temp_i)
write.csv(temp_i, "Data/GEDI/data_GEDI1_B.csv", row.names = F)
write.csv(temp_i$shot_number, "Data/GEDI/data_GEDI1_B_shot_number.csv", row.names = F)
