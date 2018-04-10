##############################################
##############################################
######                                  ######
######  Various functions required for  ######
######        FORMAS project            ######
######     about climate change         ######
######        and fragmentation         ######
######                                  ######
##############################################
##############################################


fragstat <- function(points, LC_SNH, width, cores){
  require(rgeos)
  require(raster)
  require(gdalUtils)
  require(doSNOW)
  
  buf <- gBuffer(spgeom = points, width = width, byid=T)
  
  # listBuf <- replicate(length(buf), tempfile(pattern = "", fileext = ".shp"))
  # 
  # for(i in c(1:length(buf))){
  #   shapefile(split(buf, buf$Site)[[i]], listBuf[i], overwrite = T)
  # }
  
  cropped_ras <- tempfile(pattern = "", fileext = ".tif")
  
  gdalwarp(srcfile = LC_SNH, 
           dstfile = cropped_ras, 
           te = extend(extent(buf), 5000)[c(1,3,2,4)],
           overwrite = T)
  
  landFunc <- function(buffer, LC){
    
    output_ras <- tempfile(pattern = "", fileext = ".tif")
    output_shape <- tempfile(pattern = "", fileext = ".shp")
    
    shapefile(buffer, output_shape, overwrite = T )
    
    gdalwarp(srcfile = LC, 
             dstfile = output_ras, 
             crop_to_cutline = T, 
             cutline = output_shape, 
             dstalpha = T,
             overwrite = T)
    
    r <- raster(output_ras)
    
    ai = length(Which(r == 1, cells = TRUE))
    
    if(length(r[!is.na(r)]) - ai > 1 & ai > 0){
      
      a <- adjacent(r, 1:ncell(r), 4, pairs=TRUE)
      tb <- table(r[a[,1]], r[a[,2]])
      
      gii = tb[2,2]
      gik = sum(tb[2,])
      
      P = ai/length(r[!is.na(r)])
      
      # n <- sqrt(fun(ai))
      # 
      # m = ai - n^2
      # 
      # minp = ifelse(m == 0, 4*n, ifelse(n^2 < ai & ai <= n*(1+n), 4*n + 2, 4*n+4))
      
      G = gii / gik
      
      CLUMPY = ifelse(G < P & P < .5,(G - P)/P, (G - P)/(1 - P))
      return(cbind.data.frame(Site = buffer$Site, CLUMPY = CLUMPY, PLAND = P))
    }
  }
  
  
  cl <- makeCluster(cores, type = "SOCK")
  clusterExport(cl, list("landFunc", "shapefile", "raster", "gdalwarp", "cropped_ras", "Which", "adjacent", "ncell"), 
                envir=environment())
  registerDoSNOW(cl)
  
  iterations <- length(split(buf, buf$Site))
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  frag_res <- foreach(x = split(buf, buf$Site), .combine = rbind, .options.snow = opts) %dopar% {
    landFunc(x, cropped_ras)
  }
  close(pb)
  stopCluster(cl)
  
  return(frag_res)
}


sti <- function(distri, temperature, output){
  require(raster)
  require(rgdal)
  sti_list <- c()
  for(i in 1:length(distri)){
    
    cat(paste("Species:", distri@data$AltName[i], "\n"))
    
    distri_shp <- distri[distri@data$AltName[i],] # select species
    
    sti.temp <- mean(unlist(extract(temperature, distri_shp)), na.rm =T) # calculate STI
    sti_sd.temp <- sd(unlist(extract(temperature, distri_shp)), na.rm =T) # calculate SD of STI
    
    cat(paste("STI =", round(sti.temp, 2), "deg. C", "\n\n"))
    sti_list <- rbind.data.frame(sti_list, cbind.data.frame(Species =  distri@data$AltName[i], STI = sti.temp, STI_sd = sti_sd.temp))
    
    write.csv(sti_list, output, row.names = F, append = T)
    
    rm(distri_shp)
    gc()
  }
  return(sti_list)
}


extractRTS <- function(rts, sites, fun, id.col){
  require(qdapRegex)
  require(rts)
  res <- c()
  for (i in 1:nlayers(rts@raster)){
    if(is.data.frame(sites)){
      site.name <- sites[,id.col]
    }else{
      site.name <- names(sites)
    }
    temp <- extract(rts[[i]], sites)
    res <- rbind.data.frame(res,
                            cbind.data.frame(Site = site.name, 
                                             Year = as.numeric(ex_between(names(rts[[i]]), "X", ".")[[1]]), 
                                             x = unlist(lapply(temp, fun))))
  }
  return(res)
}


sp_site_occupancy_trend <- function(data){
  res <- c()
  for(i in 1:length(unique(data$Site))){
    # cat(paste("\r", unique(data$Site)[i]))
    cat(paste("\r",round(i / length(unique(data$Site)), 2)*100,"%"))
    data.temp <- data %>% dplyr:::filter(Site == unique(data$Site)[i])
    for(j in 1:length(unique(data.temp$Species))){
      cat(".")
      if(any((data.temp %>% dplyr:::filter(Species == unique(data.temp$Species)[j]))$n != 1)){
        m <- glm(n ~ Year ,
                 family = binomial, 
                 data = data.temp %>% dplyr:::filter(Species == unique(data.temp$Species)[j]))
        
        res <- rbind.data.frame(res,
                                cbind.data.frame(Site = unique(data$Site)[i], Species = unique(data.temp$Species)[j], 
                                                 t(summary(m)$coefficients[2,c(1,2,4)]),
                                                 nYear = length(unique((data.temp %>% dplyr:::filter(Species == unique(data.temp$Species)[j]))$Year))))
      } else{
        res <- rbind.data.frame(res,
                                cbind.data.frame(Site = unique(data$Site)[i], Species = unique(data.temp$Species)[j], 
                                                 Estimate = 0, 'Std. Error' = 0, 'Pr(>|z|)' = 0,
                                                 nYear = length(unique((data.temp %>% dplyr:::filter(Species == unique(data.temp$Species)[j]))$Year))))
      }
    }
  }
  return(res)
}


turnover2 <- function(dat){
  output <- c()
  for(i in 1:length(unique(dat$Site))){
    data.temp <- dat %>% dplyr::filter(Site == unique(dat$Site)[i])
    print(paste(i, "/", length(unique(dat$Site))))
    if(length(unique(data.temp$Year)) > 1){
      for(j in 2:length(unique(data.temp$Year))){
        
        nSp <- union(data.temp %>% ungroup () %>% dplyr::filter(Year == unique(data.temp$Year)[j]) %>% dplyr::select(Species),
                     data.temp %>% ungroup () %>% dplyr::filter(Year == unique(data.temp$Year)[j-1]) %>% dplyr::select(Species)) %>% 
          summarise(n = length(unique(Species)))
        
        nSp.dis <- setdiff(data.temp %>% ungroup () %>% dplyr::filter(Year == unique(data.temp$Year)[j]) %>% dplyr::select(Species),
                           data.temp %>% ungroup () %>% dplyr::filter(Year == unique(data.temp$Year)[j-1])%>% dplyr::select(Species)) %>% 
          summarise(n = length(unique(Species)))
        
        nSp.app <- setdiff(data.temp %>% ungroup () %>% dplyr::filter(Year == unique(data.temp$Year)[j-1]) %>% dplyr::select(Species),
                           data.temp %>% ungroup () %>% dplyr::filter(Year == unique(data.temp$Year)[j]) %>% dplyr::select(Species)) %>% 
          summarise(n = length(unique(Species)))
        
        output <- rbind.data.frame(output,
                                   cbind.data.frame(Site = unique(dat$Site)[i], Year = unique(data.temp$Year)[j], 
                                                    total.turnover = as.vector((nSp.dis + nSp.app)/nSp)$n,
                                                    appearance = as.vector(nSp.app / nSp)$n,
                                                    disappearance = as.vector(nSp.dis / nSp)$n, 
                                                    total.species = as.vector(nSp)$n))
      }
    }
  }
  return(output)
}


predict_raster <- function(model, scaleList, n = 100){
  require(alphahull)
  require(raster)
  
  newdata <- cbind.data.frame(X = median(model$data$X), Y = median(model$data$Y), Habitat = c("Open", "Forest"),
                              expand.grid(PLAND = seq(min(model$data$PLAND), max(model$data$PLAND), 
                                                      length.out = n),
                                          CLUMPY = seq(min(model$data$CLUMPY), max(model$data$CLUMPY), 
                                                       length.out = n),
                                          Year = seq(min(model$data$Year), max(model$data$Year), 
                                                     length.out = 10)))
  newdata$Habitat <- as.character(newdata$Habitat)
  
  pred <- predict(model, newdata, level = 0)
  pred <- cbind.data.frame(newdata, pred = pred) %>% dplyr::filter(Habitat == "Open")
  pred <- pred %>% group_by(PLAND, CLUMPY) %>% 
    mutate(Year = Year * scaleList$scale["Year"] + scaleList$center["Year"]) %>%
    do(trend = lm(pred ~ Year, data = .)) %>% tidy(trend) %>% dplyr::filter(term == "Year")
  pred <- rasterFromXYZ(pred[,c("PLAND", "CLUMPY", "estimate")])
  
  ah <- ahull(alpha = 2.5, x = unique(model$data[,c("PLAND", "CLUMPY")]))
  ah <- a2shp(ah)
  
  pred <- mask(pred, ah)
  
  pred <- as.data.frame(rasterToPoints(pred))
  pred$x <- pred$x * scaleList$scale["PLAND"] + scaleList$center["PLAND"]
  pred$y <- pred$y * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"]
  
  return(pred)
}

predict_raster2 <- function(model, scaleList, n = 100){
  require(alphahull)
  require(raster)
  
  newdata <- expand.grid(X = median(model$frame$X), Y = median(model$frame$Y),
                         STI_rel = c(- 1, 1), 
                         PC11 =  c(model$frame[which.min(abs(model$frame$PC11 - (-1))),"PC11"],
                                   model$frame[which.min(abs(model$frame$PC11 - 0)),"PC11"],
                                   model$frame[which.min(abs(model$frame$PC11 - 1)),"PC11"]),
                         Habitat = "Open",country = "NL", 
                         PLAND = seq(min(model$frame$PLAND), max(model$frame$PLAND), 
                                     length.out = n),
                         CLUMPY = seq(min(model$frame$CLUMPY), max(model$frame$CLUMPY), 
                                      length.out = n),
                         Year = seq(min(model$frame$Year), max(model$frame$Year), 
                                    length.out = 10))
  
  newdata <- newdata %>% mutate(PC12 = ifelse(PC11 == model$frame[which.min(abs(model$frame$PC11 - (-1))),"PC11"],
                                              model$frame[which.min(abs(model$frame$PC11 - (-1))),"PC12"],
                                              ifelse(PC11 == model$frame[which.min(abs(model$frame$PC11)),"PC11"],
                                                     model$frame[which.min(abs(model$frame$PC11)),"PC12"],
                                                     model$frame[which.min(abs(model$frame$PC11 - 1)),"PC12"])))
  
  pred <- predict(model, newdata, re.form = NA, type = "response")
  pred <- cbind.data.frame(newdata, pred = pred)
  pred <- pred %>% group_by(PLAND, CLUMPY, STI_rel, PC11, PC12) %>% 
    mutate(Year = Year * scaleList$scale["Year"] + scaleList$center["Year"] ) %>%
    do(trend = lm(pred ~ Year, data = .)) %>% tidy(trend) %>% dplyr::filter(term == "Year")
  
  ah <- ahull(alpha = 2.5, x = unique(model$frame[,c("PLAND", "CLUMPY")]))
  ah <- a2shp(ah)
  
  pred <- foreach(i = unique(pred$STI_rel), .combine = rbind) %:%
    foreach(j = unique(pred$PC11), .combine = rbind) %do% {
      pred_temp <- pred %>% dplyr::filter(STI_rel == i, PC11 == j)
      res_pred <- rasterFromXYZ(pred_temp[,c("PLAND", "CLUMPY", "estimate")])
      res_pred <- mask(res_pred, ah)
      res_pred <- cbind.data.frame(STI_rel = i, PC1 = j, as.data.frame(rasterToPoints(res_pred)))
      return(res_pred)
    }
  
  pred$x <- pred$x * scaleList$scale["PLAND"] + scaleList$center["PLAND"]
  pred$y <- pred$y * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"]
  pred$STI_rel <- ifelse(pred$STI_rel == -1, "Low STI", "High STI")
  pred$PC1 <- ifelse(pred$PC1 == model$frame[which.min(abs(model$frame$PC11 - (-1))),"PC11"], 
                     "Low disperal", 
                     ifelse(pred$PC1 == model$frame[which.min(abs(model$frame$PC11)),"PC11"],
                            "Medium dispersal", "High disperal"))
  
  return(pred)
}


label_wrap <- function(x) {
  lapply(x,function(x){paste(x/1000, "km")})
}  


a2shp <- function(x, increment=360, rnd=10, proj4string=CRS(as.character(NA))){
  require(alphahull)
  require(maptools)
  if (class(x) != "ahull"){
    stop("x needs to be an ahull class object")
  }
  # Extract the edges from the ahull object as a dataframe
  xdf <- as.data.frame(x$arcs)
  # Remove all cases where the coordinates are all the same      
  xdf <- subset(xdf,xdf$r > 0)
  res <- NULL
  if (nrow(xdf) > 0){
    # Convert each arc to a line segment
    linesj <- list()
    prevx<-NULL
    prevy<-NULL
    j<-1
    for(i in 1:nrow(xdf)){
      rowi <- xdf[i,]
      v <- c(rowi$v.x, rowi$v.y)
      theta <- rowi$theta
      r <- rowi$r
      cc <- c(rowi$c1, rowi$c2)
      # Arcs need to be redefined as strings of points. Work out the number of points to allocate in this arc segment.
      ipoints <- 2 + round(increment * (rowi$theta / 2),0)
      # Calculate coordinates from arc() description for ipoints along the arc.
      angles <- anglesArc(v, theta)
      seqang <- seq(angles[1], angles[2], length = ipoints)
      x <- round(cc[1] + r * cos(seqang),rnd)
      y <- round(cc[2] + r * sin(seqang),rnd)
      # Check for line segments that should be joined up and combine their coordinates
      if (is.null(prevx)){
        prevx<-x
        prevy<-y
      } else if (x[1] == round(prevx[length(prevx)],rnd) && y[1] == round(prevy[length(prevy)],rnd)){
        if (i == nrow(xdf)){
          #We have got to the end of the dataset
          prevx<-append(prevx,x[2:ipoints])
          prevy<-append(prevy,y[2:ipoints])
          prevx[length(prevx)]<-prevx[1]
          prevy[length(prevy)]<-prevy[1]
          coordsj<-cbind(prevx,prevy)
          colnames(coordsj)<-NULL
          # Build as Line and then Lines class
          linej <- Line(coordsj)
          linesj[[j]] <- Lines(linej, ID = as.character(j))
        } else {
          prevx<-append(prevx,x[2:ipoints])
          prevy<-append(prevy,y[2:ipoints])
        }
      } else {
        # We have got to the end of a set of lines, and there are several such sets, so convert the whole of this one to a line segment and reset.
        prevx[length(prevx)]<-prevx[1]
        prevy[length(prevy)]<-prevy[1]
        coordsj<-cbind(prevx,prevy)
        colnames(coordsj)<-NULL
        # Build as Line and then Lines class
        linej <- Line(coordsj)
        linesj[[j]] <- Lines(linej, ID = as.character(j))
        j<-j+1
        prevx<-NULL
        prevy<-NULL
      }
    }
    # Promote to SpatialLines
    lspl <- SpatialLines(linesj)
    # Convert lines to polygons
    # Pull out Lines slot and check which lines have start and end points that are the same
    lns <- slot(lspl, "lines")
    polys <- sapply(lns, function(x) { 
      crds <- slot(slot(x, "Lines")[[1]], "coords")
      identical(crds[1, ], crds[nrow(crds), ])
    }) 
    # Select those that do and convert to SpatialPolygons
    polyssl <- lspl[polys]
    list_of_Lines <- slot(polyssl, "lines")
    sppolys <- SpatialPolygons(list(Polygons(lapply(list_of_Lines, function(x) { Polygon(slot(slot(x, "Lines")[[1]], "coords")) }), ID = "1")), proj4string=proj4string)
    # Create a set of ids in a dataframe, then promote to SpatialPolygonsDataFrame
    hid <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "ID"))
    areas <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "area"))
    df <- data.frame(hid,areas)
    names(df) <- c("HID","Area")
    rownames(df) <- df$HID
    res <- SpatialPolygonsDataFrame(sppolys, data=df)
    res <- res[which(res@data$Area > 0),]
  }  
  return(res)
}
