library(dplyr)

#########################
### remove duplicates ###
#########################

folder <- "../Data/Butterflies - Netherlands/"

# load
sites <- read.csv(paste0(folder, "Sites_NL_ETRS89_landcover.csv"))
sites$uniqueXY <- paste0(sites$X, "_", sites$Y)

# write duplicates
sites.dup <- sites %>% group_by(uniqueXY) %>% filter(n()>1)
list.sites.dup <- split(sites.dup,sites.dup$uniqueXY)

sink(paste0(folder, "Duplicated_sites.txt"))
for(i in unique(sites.dup$uniqueXY)){
  tempDup <- list.sites.dup[[i]]
  cat(as.character(tempDup$Site))
  cat("\n")
}
sink()

dup.sites <- read.table(paste0(folder, "Duplicated_sites.txt"), sep = " ", fill = TRUE)
replace.sites <- c()
for(k in 1:nrow(dup.sites)){
  for(l in 2:ncol(dup.sites)){
    if(!is.na(dup.sites[k,l])){
      replace.sites <- rbind.data.frame(replace.sites, cbind.data.frame(from = dup.sites[k,l], to = dup.sites[k,1]))
    }
  }
}


#optional if Finland (very close sites)
replace.sites <- rbind.data.frame(replace.sites, cbind.data.frame(from = 265, to = 80))

write.table(replace.sites, paste0(folder, "Duplicated_sites.txt"), row.names = F, quote = F)
