library(adespatial)

# sites data
sites_NL <- read.csv(file = "../Data/Butterflies - Netherlands/Sites_NL_ETRS89_landcover.csv")
sites_FIN <- read.csv(file = "../Data/Butterflies - Finland/Sites_FIN_ETRS89_landcover.csv")

# merge sites
sites <- rbind(cbind(country = "NL", sites_NL), cbind(country = "FIN", sites_FIN))
sites <- subset(sites, sites$Scale == 1000)

# compute dbMEM
dbMEM.sites <- dbmem(sites[,c("X", "Y")])
plot(dbMEM.sites, sites[,c("X", "Y")])

dbMEM <- cbind.data.frame(Site = paste0(sites$Site, "_", sites$country), dbMEM = dbMEM.sites$MEM1)
write.csv(dbMEM, "../Data/dbMEM.csv", row.names = F)
