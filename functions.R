sti <- function(species, distri, temperature){
  sti_list <- c()
  for(i in species){
    cat(paste("Species:", i, "\n"))
    nb.sp <- grep(i, distri)
    
    if(length(nb.sp)>0){
      cat(paste("Found in shapefiles", "\n"))
      distri_shp <- readOGR(distri[nb.sp][1], verbose = F) # load shapefile
      distri_shp <- distri_shp[distri_shp@data$SEASONAL %in% season,] # select resident and breeding ranges
      
      sti <- mean(unlist(extract(temperature, SpatialPolygons(distri_shp@polygons))), na.rm=T) # calculate STI
      cat(paste("STI =", round(sti, 2), "deg. C", "\n\n"))
      sti_list <- rbind.data.frame(sti_list, cbind.data.frame(Species = i, STI = sti))
    }else{
      cat(paste("Not found in shapefiles", "\n\n"))
      sti_list <- rbind.data.frame(sti_list, cbind.data.frame(Species = i, STI = NA))
    }
  }
  return(sti_list)
}
