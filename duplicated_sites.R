library(dplyr)

#########################
### remove duplicates ###
#########################

# load
sites_NL <- read.csv(file = "../Data/Butterflies - Netherlands/Sites_NL_ETRS89_landcover.csv")
sites_NL$uniqueXY <- paste0(sites_NL$X, "_", sites_NL$Y)

# write duplicates
sites_NL.dup <- sites_NL %>% group_by(uniqueXY) %>% filter(n()>1)
list.sites_NL.dup <- split(sites_NL.dup,sites_NL.dup$uniqueXY)

sink("../Data/Butterflies - Netherlands/Duplicated_NL2.txt")
for(i in unique(sites_NL.dup$uniqueXY)){
  tempDup <- list.sites_NL.dup[[i]]
  cat(as.character(tempDup$Site))
  cat("\n")
}
sink()

# write file without duplicates
sites_NL <- sites_NL[!duplicated(sites_NL$uniqueXY),-5]
write.csv(sites_NL, file = "../Data/Butterflies - Netherlands/Sites_NL_ETRS89_landcover.csv", row.names = F)

########################
### modify CTI files ###
###     DONT'T!!     ###
########################

# load
cti_AB_NL <- read.csv2("../Data/Butterflies - Netherlands/Original/cti_site_year_abundance_2017-05-19.csv", dec = ".")[,-3]
cti_AB_NL_mod <- cti_AB_NL

duplicates <- read.table("../Data/Butterflies - Netherlands/Duplicated_NL.txt", sep = "\t")

for(i in 1:nrow(duplicates)){
  dup.temp <- unlist(strsplit(as.character(duplicates[i,]), " ")) 
  for(j in 2:length(dup.temp)){
    cti_AB_NL_mod[cti_AB_NL_mod$Site %in% dup.temp[j],"Site"] <- dup.temp[j-1]
  }
}

cti_AB_NL_mod[cti_AB_NL_mod$Site %in% unique(sites_NL.dup$Site),]
cti_AB_NL_mod[cti_AB_NL_mod$Site %in% 13,]


con <- file("../Data/Butterflies - Netherlands/Duplicated_NL.txt")
open(con)
l <- length(readLines(con))
for(i in 1:l) {
  temp <- readLines(con,  n = i)
  print(temp)
}
close(con)

