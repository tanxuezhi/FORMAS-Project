sti <- function(distri, temperature, output){
  require(raster)
  require(rgdal)
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

extractFrag <- function(folder, dup.sites, sites, verbose = F){
  require(stringr)
  files_frag <- list.files(folder, pattern = ".class", full.names = T)
  frag_all_scales <- c()
  for(i in 1:length(files_frag)){
    frag <- readLines(files_frag[i])
    frag <- gsub(pattern = "N/A", replace = " ", x = frag)
    writeLines(frag, con=files_frag[i])
    
    frag <- read.table(files_frag[i], h=T, sep = ",", dec=".")
    frag[,1] <- as.numeric(gsub("point_", "", frag[,1]))
    frag <- frag[,c("LID","NLSI","PLAND", "CLUMPY", "PAFRAC", "PD")]
    
    if(identical(sites, sites_SWE)){
      corres <- cbind.data.frame(Site = sites$Site, SiteID = 1:nrow(sites))
      frag[,1] <- merge(frag, corres, by.x = "LID", by.y = "SiteID")[,"Site"]
    }
    for(j in unique(sites$Site)){
      if(verbose == T){
        print(paste("i =",i,"; j =", j))
      }
      if(all(!j == frag$LID)){
        if(any(j == dup.sites[,1])){
          replace.site <- dup.sites[dup.sites$from == j,2]
          if(any(frag$LID %in% replace.site)){
            frag <- rbind.data.frame(frag, cbind.data.frame(LID = j, frag[frag$LID == replace.site,-1]))
          }else{
            frag <- rbind.data.frame(frag, cbind.data.frame(LID = j, NLSI = NA, PLAND = 0, CLUMPY = NA, PAFRAC = 0, PD = 0))
          }
        }else{
          if(any(j == dup.sites[,2])){
            replace.site <- dup.sites[dup.sites$to == j,1]
            if(any(frag$LID %in% replace.site)){
              frag <- rbind.data.frame(frag, cbind.data.frame(LID = j, frag[frag$LID == replace.site,-1]))
            }else{
              frag <- rbind.data.frame(frag, cbind.data.frame(LID = j, NLSI = NA, PLAND = 0, CLUMPY = NA, PAFRAC = 0, PD = 0))
            }
          }else{
            frag <- rbind.data.frame(frag, cbind.data.frame(LID = j, NLSI = NA, PLAND = 0, CLUMPY = NA, PAFRAC = 0, PD = 0))
          }
        }
      }
    }
    frag_all_scales <- rbind.data.frame(frag_all_scales, cbind.data.frame(Scale = str_sub(files_frag[i], -11, -7), frag))
  }
  frag_all_scales$Scale <- gsub("_","",frag_all_scales$Scale)
  return(frag_all_scales)
}


scaleTest <- function(dataToUse, AIC = T, plot = T, main = NULL){
  require(lmerTest)
  require(lme4)
  require(MuMIn)
  require(ggplot2)
  
  if(AIC == T){
    AIC.tab <- c()
    for(i in unique(dataToUse$Scale)){
      if(any(names(dataToUse) == "country")){
        m_frag_std <- gamm(cti ~ CLUMPY * PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                           random = list(country = ~1, Site = ~1), correlation = corCAR1(form = ~Year|Site),
                           data = data.frame(stdize(subset(dataToUse, dataToUse$Scale == i), prefix = F)))
      }else{
        m_frag_std <- gamm(cti ~ CLUMPY * PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                           random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                           data = data.frame(stdize(subset(dataToUse, dataToUse$Scale == i), prefix = F)))
      }
      AIC.tab <- rbind.data.frame(AIC.tab, cbind.data.frame(Scale = i, AICc = AICc(m_frag_std)))
    }
    
    if(plot == T){
      p <- ggplot(data = AIC.tab,aes(x= Scale, y = AICc)) +
        geom_point(size = 2) + 
        geom_line() + ggtitle(main) + 
        scale_x_continuous("Spatial scale (m)")
      
      print(p)
      
    }
    return(AIC.tab)
  }else{
    scale.dat <- c()
    for(i in unique(dataToUse$Scale)){
      if(any(names(dataToUse) == "country")){
        m_frag_std <- gamm(cti ~ CLUMPY * PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                           random = list(country = ~1, Site = ~1), correlation = corCAR1(form = ~Year|Site),
                           data = data.frame(stdize(subset(dataToUse, dataToUse$Scale == i), prefix = F)))
      }else{
        m_frag_std <- gamm(cti ~ CLUMPY * PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                           random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                           data = data.frame(stdize(subset(dataToUse, dataToUse$Scale == i), prefix = F)))
      }
      sum.m <- summary(m_frag_std$gam)
      sum.m <- cbind(sum.m$p.coeff, sum.m$se, sum.m$p.pv)
      sum.m <- sum.m[-grep("X,Y", rownames(sum.m)),]
      sum.m <- sum.m[-grep("LABEL", rownames(sum.m)),]
      sum.m <- as.data.frame(sum.m)
      sum.m$Variable <- rownames(sum.m)
      rownames(sum.m) <- NULL
      colnames(sum.m)[1:3] <- c("Estimate", "SE", "pvalue")
      
      scale.dat <- rbind.data.frame(scale.dat, cbind.data.frame(Scale = i, 
                                                                sum.m))
    }
    scale.dat$Significance <- ifelse(scale.dat$pvalue < .001, "***", ifelse(scale.dat$pvalue < .01, "**", ifelse(scale.dat$pvalue < .05, "*", ifelse(scale.dat$pvalue < .1, ".", ""))))
    
    if(plot == T){
      p <- ggplot(data = scale.dat[scale.dat$Variable %in% c("PLAND:Year", "CLUMPY:Year", "CLUMPY:PLAND:Year"),],
                  aes(x= Scale, y = Estimate, color = Variable)) +
        # geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE), width = 0) +
        geom_point(size = 2) + 
        geom_text(aes(x= Scale, y = (Estimate + .04 * diff(range(Estimate))),
                      label = Significance))+
        geom_hline(yintercept=0, lty = 2) +
        geom_line() + ggtitle(main) + 
        scale_color_manual("Interactions", 
                           values = c("#EE7600", "#698B69", "#00B2EE"), 
                           labels = c("Clumpiness x % SNH x Year", 
                                      "Clumpiness x Year",
                                      "% SNH x Year")) +
        scale_x_continuous("Spatial scale (m)")
      
      print(p)
    }
    
    return(scale.dat)
  }
}


vis.2d <- function(model, var1, var2, var3, n = 10, plot = T, origin, ...){
  require(dplyr)
  require(fields)
  
  data.mod <- model$model
  
  var1_val <- seq(from = min(data.mod[,var1]), to = max(data.mod[,var1]), length.out = n)
  var2_val <- seq(from = min(data.mod[,var2]), to = max(data.mod[,var2]), length.out = n)
  var3_val <- seq(from = min(data.mod[,var3]), to = max(data.mod[,var3]), length.out = n)
  
  newdata = expand.grid(var1_val, var2_val, var3_val)
  colnames(newdata) <- c(var1, var2, var3)
  
  vars <- all.vars(formula(model))
  vars <- vars[!vars %in% as.character(model$terms[[2]])]
  vars <- vars[!vars %in% colnames(newdata)]
  
  l = length(names(newdata))
  for(i in 1:length(vars)){
    l = l + 1
    if(is.character(model$model[,vars[i]])){
      newdata <- cbind.data.frame(newdata, names(which.max(table(model$model[,vars[i]]))))
      names(newdata)[l] <- vars[i]
    }else{
      newdata <- cbind.data.frame(newdata, median(model$model[,vars[i]]))
      names(newdata)[l] <- vars[i]
    }
  }
  
  pred <- predict(model, newdata = newdata)
  pred <- cbind.data.frame(pred, newdata)
  
  scaleList <- list(scale = attr(origin, "scaled:scale")[c(as.character(model$terms[[2]]), names(newdata))],
                    center = attr(origin, "scaled:center")[c(as.character(model$terms[[2]]), names(newdata))])
  
  pred[,"pred"] <- pred[,"pred"] * scaleList$scale[as.character(model$terms[[2]])] + scaleList$center[as.character(model$terms[[2]])]
  pred[,var1] <- pred[,var1] * scaleList$scale[var1] + scaleList$center[var1]
  pred[,var2] <- pred[,var2] * scaleList$scale[var2] + scaleList$center[var2]
  pred[,var3] <- pred[,var3] * scaleList$scale[var3] + scaleList$center[var3]
  
  pred <- pred %>% group_by(get(var2), get(var3)) %>% summarise(fit = lm(pred ~ Year)$coef[2])
  names(pred)[1:2] <- c(var2,var3)
  
  dat <- cbind.data.frame(data.mod[,var2]* scaleList$scale[var2] + scaleList$center[var2],
                          data.mod[,var3]* scaleList$scale[var3] + scaleList$center[var3])
  names(dat) <- c(var2, var3)
  
  if(plot == T){
    
    quilt.plot(pred,  nx = nrow(pred)/n, ny = nrow(pred)/n, 
               xlab = var2, ylab = var3, ...)
    points(dat, pch = 16, cex = .4)
  }
  return(pred)
}

vis.1d <- function(model, var1, var2, var3, n = 10, plot = T, origin){
  require(ggplot2)
  require(dplyr)
  
  data.mod <- model@frame
  
  var1_val <- seq(from = min(data.mod[,var1]), to = max(data.mod[,var1]), length.out = n)
  var2_val <- c(mean(data.mod[,var2]) - sd(data.mod[,var2]), median(data.mod[,var2]), mean(data.mod[,var2]) + sd(data.mod[,var2]))
  var3_val <- c(mean(data.mod[,var3]) - sd(data.mod[,var3]), mean(data.mod[,var3]) + sd(data.mod[,var3]))
  
  newdata = expand.grid(var1_val, var2_val, var3_val)
  colnames(newdata) <- c(var1, var2, var3)
  
  vars <- all.vars(formula(model))
  vars <- vars[!vars %in% as.character(attr(data.mod, "terms")[[2]])]
  vars <- vars[!vars %in% colnames(newdata)]
  
  l = length(names(newdata))
  for(i in 1:length(vars)){
    l = l + 1
    if(is.character(data.mod[,vars[i]]) | is.factor(data.mod[,vars[i]])){
      newdata <- cbind.data.frame(newdata, names(which.max(table(data.mod[,vars[i]]))))
      names(newdata)[l] <- vars[i]
    }else{
      newdata <- cbind.data.frame(newdata, median(data.mod[,vars[i]]))
      names(newdata)[l] <- vars[i]
    }
  }
  
  pred <- predict(model, newdata = newdata, se.fit = F, re.form = NA)
  pred <- cbind.data.frame(fit = pred, newdata)

  scaleList <- list(scale = attr(origin, "scaled:scale")[c(as.character(attr(data.mod, "terms")[[2]]), names(newdata))],
                    center = attr(origin, "scaled:center")[c(as.character(attr(data.mod, "terms")[[2]]), names(newdata))])
  
  pred[,"fit"] <- pred[,"fit"] * scaleList$scale[as.character(attr(data.mod, "terms")[[2]])] + scaleList$center[as.character(attr(data.mod, "terms")[[2]])]
  pred[,var1] <- pred[,var1] * scaleList$scale[var1] + scaleList$center[var1]
  pred[,var2] <- pred[,var2] * scaleList$scale[var2] + scaleList$center[var2]
  pred[,var3] <- pred[,var3] * scaleList$scale[var3] + scaleList$center[var3]
  
  if(plot == T){
    var3_names <- c("Little area of SNH","Large area of SNH"); names(var3_names) <- unique(pred[,var3])
    var2_names <- c("Highly fragmentated","Moderately fragmentated","Little fragmentated"); names(var2_names) <- unique(pred[,var2])
    
    col.var2 <- c("firebrick", "gold2", "dodgerblue"); names(col.var2) <- unique(pred[,var2])
    
    p <- ggplot(data = pred, aes(x = get(var1), y = fit, color = as.factor(get(var2)))) + geom_line() + 
      facet_grid( ~ get(var3), labeller = as_labeller(var3_names)) + 
      scale_y_continuous("Community temperature Index") +
      scale_x_continuous(var1) +
      scale_color_manual("Habitat fragmentation", labels = var2_names, values=col.var2) +
      scale_fill_manual("Habitat fragmentation", labels = var2_names, values=col.var2) +
      theme(legend.position = c(1, 1), 
            legend.justification = c(1, 1))
    print(p)
  }
  return(as.tbl(pred))
}
