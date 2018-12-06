library(tidyverse)

##### load data #####

## Site
sites_FIN <- read_csv(file = "../Data/Butterflies - Finland/Sites_FIN_ETRS89_landcover.csv")
sites_NL <- as.tbl(rio:::import(file = "../Data/Butterflies - Netherlands/Sites_NL.xlsx"))

## Traits
traits.butterflies <- as.tbl(read_csv("../Data/Butterflies - Netherlands/SpeciesTraits_WDV2014.csv"))
traits.butterflies[grep("Colias croceus", traits.butterflies$Scientific_name), "Scientific_name"] <- "Colias crocea"
traits.butterflies[grep("walbum", traits.butterflies$Scientific_name), "Scientific_name"] <- "Satyrium w-album"
traits.butterflies[grep("Neozephyrus quercus", traits.butterflies$Scientific_name), "Scientific_name"] <- "Favonius quercus"
traits.butterflies[grep("lycaon", traits.butterflies$Scientific_name), "Scientific_name"] <- "Hyponephele lycaon"
traits.butterflies[grep("Polygonia", traits.butterflies$Scientific_name), "Scientific_name"] <- "Nymphalis c-album"
traits.butterflies[grep("Inachis io", traits.butterflies$Scientific_name), "Scientific_name"] <- "Aglais io"
traits.butterflies[grep("tithonus", traits.butterflies$Scientific_name), "Scientific_name"] <- "Pyronia tithonus"
traits.butterflies[grep("alcon", traits.butterflies$Scientific_name), "Scientific_name"] <- "Phengaris alcon"

## Netherlands
data_NL <- as.tbl(readRDS("../Data/Butterflies - Netherlands/AllSpecies_reg_gam_ind_20171206_algroutes.rds")) %>%
  filter(prop_pheno_sampled > .5, YEAR > 1991)  %>% 
  left_join(traits.butterflies %>% mutate(Species = tolower(Species)), by = c("SPECIES" = "Species")) %>% 
  left_join(sites_NL, by = c("SITE" = "Site")) %>% mutate(Site = paste0(SITE, "_NL")) %>%
  mutate(X = X_km, Y = Y_km) %>%
  select(Scientific_name,Site,YEAR,regional_gam,X,Y) %>% 
  rename(Species = Scientific_name, Year = YEAR, n = regional_gam) %>% 
  filter(!is.na(X), !is.na(Species))
  
write_csv(data_NL, "../Data/butterfly_data_invasion_NL.csv")


## Finland
countsFIN1 <- read.table("../Data/Butterflies - Finland/FINLAND_Records_1999-2015.txt", sep = ";", h=T)
countsFIN2 <- read.table("../Data/Butterflies - Finland/FINLAND_Records_2016.txt", sep = ";", h=T)[,-6]
data_FIN <- rbind(countsFIN1, countsFIN2)
data_FIN$Species <- gsub("Polygonia c-album", "Nymphalis c-album", data_FIN$Species)
data_FIN$Species <- gsub("Cyaniris semiargus", "Polyommatus semiargus", data_FIN$Species)
data_FIN$Species <- gsub("Leptidea juvernica", "Leptidea sinapis", data_FIN$Species)

data_FIN <- data_FIN %>%
  group_by(Site, Species_Faunaeur, Year) %>% 
  summarise(n = max(Individuals)) %>%  filter(n > 0) %>% group_by(Site) %>%
  complete(Year, nesting(Site, Species_Faunaeur), fill = list(n = 0)) %>%
  ungroup() %>%
  left_join(sites_FIN, by = c("Site" = "Site")) %>% 
  mutate(Site = paste0(Site, "_FIN")) %>% 
  select(Species_Faunaeur,Site,Year,n,X,Y) %>% rename(Species = Species_Faunaeur)


##### merge and write #####
data <- rbind.data.frame(data_NL, data_FIN)
write_csv(data, "../Data/butterfly_data_invasion.csv")
