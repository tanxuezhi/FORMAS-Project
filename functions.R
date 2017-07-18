sti <- function(distri, temperature, output){
  sti_list <- c()
  for(i in 1:length(distri)){
    
    cat(paste("Species:", distri@data$AltName[i], "\n"))
    
    distri_shp <- distri[distri@data$AltName[i],] # select species
    
    sti.temp <- mean(unlist(extract(temperature, distri_shp)), na.rm =T) # calculate STI
    cat(paste("STI =", round(sti.temp, 2), "deg. C", "\n\n"))
    sti_list <- rbind.data.frame(sti_list, cbind.data.frame(Species =  distri@data$AltName[i], STI = sti.temp))
    
    write.csv(sti_list, output, row.names = F, append = T)
    
    rm(distri_shp)
    gc()
  }
  return(sti_list)
}

extractRTS <- function(rts, sites, fun, id.col){
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

extractFrag <- function(folder, dup.sites, sites){
  files_frag <- list.files(folder, pattern = ".class", full.names = T)
  frag_all_scales <- c()
  for(i in 1:length(files_frag)){
    frag <- readLines(files_frag[i])
    frag <- gsub(pattern = "N/A", replace = " ", x = frag)
    writeLines(frag, con=files_frag[i])
    
    frag <- read.table(files_frag[i], h=T, sep = ",", dec=".")
    frag[,1] <- as.numeric(gsub("point_", "", frag[,1]))
    frag <- frag[,c("LID","LSI","PLAND", "CLUMPY")]
    
    if(identical(sites, sites_SWE)){
      corres <- cbind.data.frame(Site = sites$Site, SiteID = 1:nrow(sites))
      frag[,1] <- merge(frag, corres, by.x = "LID", by.y = "SiteID")[,"Site"]
    }
    
    for(j in unique(sites$Site)){
      if(all(!j == frag$LID)){
        if(any(j == dup.sites[,1])){
          replace.site <- dup.sites[dup.sites$from == j,2]
          if(any(frag$LID == replace.site)){
            frag <- rbind.data.frame(frag, cbind.data.frame(LID = j, frag[frag$LID == replace.site,-1]))
          }else{
            frag <- rbind.data.frame(frag, cbind.data.frame(LID = j, LSI = 1, PLAND = 0, CLUMPY = 1))
          }
        }else{
          if(any(j == dup.sites[,2])){
            replace.site <- dup.sites[dup.sites$to == j,1]
            if(any(frag$LID == replace.site)){
              frag <- rbind.data.frame(frag, cbind.data.frame(LID = j, frag[frag$LID == replace.site,-1]))
            }else{
              frag <- rbind.data.frame(frag, cbind.data.frame(LID = j, LSI = 1, PLAND = 0, CLUMPY = 1))
            }
          }else{
            frag <- rbind.data.frame(frag, cbind.data.frame(LID = j, LSI = 1, PLAND = 0, CLUMPY = 1))
          }
        }
      }
    }
    frag_all_scales <- rbind.data.frame(frag_all_scales, cbind.data.frame(Scale = str_sub(files_frag[i], -11, -7), frag))
  }
  frag_all_scales$Scale <- gsub("_","",frag_all_scales$Scale)
  return(frag_all_scales)
}


scaleTest <- function(dataToUse, plot = T){
  scale.dat <- c()
  for(i in unique(dataToUse$Scale)){
    if(any(names(dataToUse) == "country")){
      m_frag_std <- lmer(cti ~ CLUMPY * PLAND * Year + X*Y + LABEL3 + (1|country/Site), data = stdize(subset(dataToUse, dataToUse$Scale == i), prefix = F), na.action = na.fail)
    }else{
      m_frag_std <- lmer(cti ~ CLUMPY * PLAND * Year + X*Y + LABEL3 + (1|Site), data = stdize(subset(dataToUse, dataToUse$Scale == i), prefix = F), na.action = na.fail)
    }
    scale.dat <- rbind.data.frame(scale.dat, cbind.data.frame(Scale = i, AICc = AICc(m_frag_std), coefInter = fixef(m_frag_std)["CLUMPY:PLAND:Year"]))
  }
  rownames(scale.dat) <- NULL
  if(plot == T){
    plot(scale.dat[order(scale.dat$Scale),-3], xlab = "Spatial scale (m)", ylab = "AICc", type = "o", pch = 16)
  }
  return(scale.dat)
}


vis.2d <- function(model, var1, var2, var3, n = 10, plot = T){
  data.mod <- model@frame
  
  var2_val <- seq(from = min(data.mod[,var2]), to = max(data.mod[,var2]), length.out = n)
  var3_val <- seq(from = min(data.mod[,var3]), to = max(data.mod[,var3]), length.out = n)
  
  res <- c()
  n.max <- length(var2_val) * length(var3_val)
  n <- 0
  for(i in var2_val){
    for (j in var3_val){
      n = n+1
      cat(paste(round(n/n.max*100, 2), "%, parameters :", var2, "=", round(i,2), ",", var3, "=", round(j,2), "\n"))
      p <- visreg(model, xvar = var1, cond = setNames(list(i,j), c(var2, var3)), plot = F)
      slope <- coefficients(lm(visregFit ~ get(var1), data = p$fit))[2]
      res <- bind_rows(res, tibble(var2 = i, var3 = j, slope = slope))
    }
  }
  

  if(plot == T){
    
    quilt.plot(res,  nx = length(var3_val), ny = length(var2_val), 
               xlab = var3, ylab = var2,
               main = "")
    points(get(var2) ~ get(var3), data = data.mod, pch = 16, cex = .4)
  }
  return(res)
  
}