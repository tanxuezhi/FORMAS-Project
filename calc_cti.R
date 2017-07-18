library(dplyr)

########################
##### Butterflies ######
########################

##### load STI data #####
butterflies_sti <- rio::import(file = "//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Pyrgus Armoricanus/Various documents/CLIMBER database/CLIMBER - Climatic niche characteristics of the butterflies in Europe.xlsx", which = 1L)
butterflies_sti <- butterflies_sti[butterflies_sti$measurementType == "Temperature (STI)",c(1,4)]
colnames(butterflies_sti)[2] <- "STI"
butterflies_sti$STI <- as.numeric(butterflies_sti$STI)


##### load count data #####
countsFIN1 <- read.table("../Data/Butterflies - Finland/FINLAND_Records_1999-2015.txt", sep = ";", h=T)
countsFIN2 <- read.table("../Data/Butterflies - Finland/FINLAND_Records_2016.txt", sep = ";", h=T)[,-6]

countsFIN <- rbind(countsFIN1, countsFIN2)

countsFIN$Species_Faunaeur <- gsub("Polygonia c-album", "Nymphalis c-album", countsFIN$Species_Faunaeur)
countsFIN$Species_Faunaeur <- gsub("Cyaniris semiargus", "Polyommatus semiargus", countsFIN$Species_Faunaeur)
countsFIN$Species_Faunaeur <- gsub("Leptidea juvernica", "Leptidea sinapis", countsFIN$Species_Faunaeur)

##### Extract STI #####
speciesFIN <- unique(countsFIN$Species_Faunaeur)

stiFIN <-  butterflies_sti[butterflies_sti$id %in% speciesFIN,]


##### Calculate CTI #####

maxCounts_Site_year <- countsFIN %>% group_by(Species_Faunaeur, Year, Site) %>% summarise(Individuals = max(Individuals))
maxCounts_Site_year <- left_join(maxCounts_Site_year, stiFIN, by = c("Species_Faunaeur" = "id") )

cti_AB <- maxCounts_Site_year %>% group_by(Year, Site) %>% summarise(cti = weighted.mean(STI, Individuals))
cti_PA <- maxCounts_Site_year %>% group_by(Year, Site) %>% summarise(cti = mean(STI))

##### Write CTI data #####
write.csv(cti_AB, "../Data/Butterflies - Finland/CTI_Abundance_FINLAND_1999-2016.csv", row.names = F)
write.csv(cti_PA, "../Data/Butterflies - Finland/CTI_Presence_FINLAND_1999-2016.csv", row.names = F)


##### plot trend of CTI #####

par(mfrow=c(1,2))

mean.cti.ab <- cti_AB %>% group_by(Year) %>% summarise(CTI = mean(cti))
plot(CTI~Year, mean.cti.ab, type = "o", pch = 16, main = "Abundance weighted")
abline(lm(CTI~Year, mean.cti.ab))
summary(lm(CTI~Year, mean.cti.ab))

mean.cti.pa <- cti_PA %>% group_by(Year) %>% summarise(CTI = mean(cti))
plot(CTI~Year, mean.cti.pa, type = "o", pch = 16, main = "Presence only")
abline(lm(CTI~Year, mean.cti.pa))
summary(lm(CTI~Year, mean.cti.pa))



########################
######## Birds #########
########################

##### load STI data #####
birds_sti <- read.csv("../BOTW/STI.csv")
list.spSWE <- read.csv("../Data/Birds - Sweden/Species_list.csv")
birds_sti <- merge(birds_sti, list.spSWE, by.x = "Species", by.y = "AltName")[,1:3]

##### load count data #####
countsSWE1 <- rio::import("../Data/Birds - Sweden/public_totalstandard_Ia.xlsx")
countsSWE2 <- rio::import("../Data/Birds - Sweden/public_totalstandard_Ib.xlsx")
countsSWE3 <- rio::import("../Data/Birds - Sweden/public_totalstandard_IIa.xlsx")
countsSWE4 <- rio::import("../Data/Birds - Sweden/public_totalstandard_IIb.xlsx")

countsSWE <- as.tbl(bind_rows(countsSWE1, countsSWE2, countsSWE3, countsSWE4))
rm(list = c("countsSWE1", "countsSWE2", "countsSWE3", "countsSWE4"))
countsSWE$art <- as.numeric(countsSWE$art)

##### Extract STI #####
countsSWE <- countsSWE %>% filter(art %in% 1:645)
countsSWE <- countsSWE[,c(2,4,5,22,23)]

par(mfrow=c(1,2))
plot(countsSWE %>% group_by(Year = yr) %>% summarise('No. monitored sites' = length(unique(karta))))
plot(countsSWE %>% group_by(Year = yr) %>% summarise("No. species detected" = length(unique(art))))

# countsSWE.sel <- filter(countsSWE, karta == "02C2H")[,-c(6:21)]

countsSWE_sum_year <- left_join(countsSWE, birds_sti, by = c("art" = "art") )

countsSWE_sum_year %>% filter(is.na(STI))

cti_AB_SWE <- countsSWE_sum_year %>% group_by(Year = yr, Site = karta) %>% summarise(cti = weighted.mean(STI, pkind + lind, na.rm = T))
cti_PA_SWE <- countsSWE_sum_year %>% group_by(Year = yr, Site = karta) %>% summarise(cti = mean(STI, na.rm = T))

##### Write CTI data #####
write.csv(cti_AB_SWE, "../Data/Birds - Sweden/CTI_abundance_Sweden_1996-2016.csv", row.names = F)
write.csv(cti_PA_SWE, "../Data/Birds - Sweden/CTI_presence_Sweden_1996-2016.csv", row.names = F)


###### tests plots ######
cti_change <- cti_AB_SWE %>%  group_by(Site) %>%  do(le_lin_fit(.))
cti_change <- merge(cti_change, sites_SWE)

library(gstat)
library(raster)
library(rgdal)
CLC_SNH_SWE <- raster("../Landcover/Sweden/SNH_SWE.tif")
CLC_SNH_SWE <- aggregate(CLC_SNH_SWE, 20)

download.file(url = 'http://biogeo.ucdavis.edu/data/diva/adm/SWE_adm.zip', 
              destfile = 'sweden.zip')
unzip(zipfile = 'sweden.zip')
sweden <- readOGR('SWE_adm0.shp')
sweden <- spTransform(sweden, CRS("+init=epsg:3035"))

mg <- gstat(id = "Site", formula = slope~1, locations = ~x+y, data=cti_change, nmax = 4)
z <- interpolate(CLC_SNH_SWE, mg)
z <- mask(z, sweden)
spplot(z)
