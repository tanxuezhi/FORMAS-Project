library(dplyr)

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


##### plot evolution of CTI #####
par(mfrow=c(1,2))

mean.cti.ab <- cti_AB %>% group_by(Year) %>% summarise(CTI = mean(cti))
plot(CTI~Year, mean.cti.ab, type = "o", pch = 16, main = "Abundance weighted")
abline(lm(CTI~Year, mean.cti.ab))
summary(lm(CTI~Year, mean.cti.ab))

mean.cti.pa <- cti_PA %>% group_by(Year) %>% summarise(CTI = mean(cti))
plot(CTI~Year, mean.cti.pa, type = "o", pch = 16, main = "Presence only")
abline(lm(CTI~Year, mean.cti.pa))
summary(lm(CTI~Year, mean.cti.pa))

