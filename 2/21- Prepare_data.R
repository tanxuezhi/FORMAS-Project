library(tidyverse)

#########################
######    Birds    ######
#########################

## load data ##
# monitoring data
list.spSWE <- read_csv("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Data/Birds - Sweden/Species_list.csv")

data.1 <- rio::import("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Data/Birds - Sweden/public_totalstandard_Ia.xlsx")
data.2 <- rio::import("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Data/Birds - Sweden/public_totalstandard_Ib.xlsx")
data.3 <- rio::import("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Data/Birds - Sweden/public_totalstandard_IIa.xlsx")
data.4 <- rio::import("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Data/Birds - Sweden/public_totalstandard_IIb.xlsx")

data <- as.tbl(bind_rows(data.1, data.2, data.3, data.4))
rm(data.1, data.2, data.3, data.4)


data$art <- as.numeric(data$art)
data <- data %>% filter(art < 645)
data <- left_join(data, list.spSWE, by= "art")

data <- data %>% dplyr::select(2,4,22,23,25) %>% filter(!is.na(latin), lind + pkind > 0) %>% mutate(presence = 1, abundance = pkind + lind) %>% dplyr::select(- pkind, -lind)


data <- data %>% group_by(karta) %>%
  complete(yr, latin) %>%
  mutate(presence = ifelse(is.na(presence), 0, presence),
         abundance = ifelse(is.na(abundance), 0, abundance))

sites.SWE <- rio::import("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Data/Birds - Sweden/Koordinater standardrutternas samtliga punkter.xlsx")
sites.SWE <- sites.SWE %>% group_by(karta) %>% summarise(x = mean(rt90_o), y = mean(rt90_n)) 

data <- data %>% left_join(sites.SWE, by = "karta")

write_csv(data, "../SWE_PA_data.csv")


#traits data
traits <- as.tbl(read.table("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Birds trait database/Life-history characteristics of European birds.txt", sep = "\t", h = T))

data %>% filter(!latin %in% traits$Species) %>% distinct(latin) %>% print(n = 30)

traits$Species <- gsub("Spatula clypeata", "Anas clypeata", traits$Species)
traits$Species <- gsub("Mareca penelope", "Anas penelope", traits$Species)
traits$Species <- gsub("Spatula querquedula", "Anas querquedula", traits$Species)
traits$Species <- gsub("Mareca strepera", "Anas strepera", traits$Species)
traits$Species <- gsub("Linaria cannabina", "Carduelis cannabina", traits$Species)
traits$Species <- gsub("Acanthis flammea hornemanni", "Carduelis hornemanni", traits$Species)
traits$Species <- gsub("Acanthis flammea", "Carduelis flammea", traits$Species)
traits$Species <- gsub("Linaria flavirostris", "Carduelis flavirostris", traits$Species)
traits$Species <- gsub("Spinus spinus", "Carduelis spinus", traits$Species)
traits$Species <- gsub("Eudromias morinellus", "Charadrius morinellus", traits$Species)
traits$Species <- gsub("Corvus corone", "Corvus corone cornix", traits$Species)
traits$Species <- gsub("Dryobates minor", "Dendrocopos minor", traits$Species)
traits$Species <- gsub("Lagopus muta", "Lagopus mutus", traits$Species)
traits$Species <- gsub("Hydrocoloeus minutus", "Larus minutus", traits$Species)
traits$Species <- gsub("Calidris falcinellus", "Limicola falcinellus", traits$Species)
traits$Species <- gsub("Cyanecula svecica", "Luscinia svecica", traits$Species)
traits$Species <- gsub("Mergellus albellus", "Mergus albellus", traits$Species)
traits$Species <- gsub("Poecile cinctus", "Parus cinctus", traits$Species)
traits$Species <- gsub("Poecile montanus", "Parus montanus", traits$Species)
traits$Species <- gsub("Poecile palustris", "Parus palustris", traits$Species)
traits$Species <- gsub("Calidris pugnax", "Philomachus pugnax", traits$Species)
traits$Species <- gsub("Regulus ignicapilla", "Regulus ignicapillus", traits$Species)
traits$Species <- gsub("Sternula albifrons", "Sterna albifrons", traits$Species)
traits$Species <- gsub("Hydroprogne caspia", "Sterna caspia", traits$Species)
traits$Species <- gsub("Thalasseus sandvicensis", "Sterna sandvicensis", traits$Species)
traits$Species <- gsub("Morus bassanus", "Sula bassana", traits$Species)
traits$Species <- gsub("Lyrurus tetrix", "Tetrao tetrix", traits$Species)

traits <- bind_rows(traits, traits %>% filter(grepl("Loxia", traits$Species)) %>% 
                      summarise_all(mean) %>% mutate(Species = "Loxia species"))

data %>% filter(!latin %in% traits$Species) %>% distinct(latin) %>% print(n = 30)

write_csv(traits, "//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Birds trait database/traits_birds_corrected.csv")



# phylogenetic distance
# from https://birdtree.org and Ericson all species dataset
tree <- read.nexus("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Birds phylogeny/output.nex")

phyl_dist <- lapply(tree, cophenetic)
phyl_dist <- apply(simplify2array(phyl_dist), 1:2, mean)
phyl_dist <- subset(reshape2::melt(phyl_dist), value!=0) %>% rename(Gained = Var1, Lost = Var2, phyl_dist = value) %>% as.tbl
phyl_dist <- bind_rows(phyl_dist, phyl_dist %>% rename(Gained = Lost, Lost = Gained))
phyl_dist <- phyl_dist %>% mutate(Gained = gsub("_", " ", Gained), Lost = gsub("_", " ", Lost))

phyl_dist$Gained <- gsub("Corvus corone", "Corvus corone cornix", phyl_dist$Gained); phyl_dist$Lost <- gsub("Corvus corone", "Corvus corone cornix", phyl_dist$Lost)
phyl_dist$Gained <- gsub("Lagopus muta", "Lagopus mutus", phyl_dist$Gained); phyl_dist$Lost <- gsub("Lagopus muta", "Lagopus mutus", phyl_dist$Lost)
phyl_dist$Gained <- gsub("Mergellus albellus", "Mergus albellus", phyl_dist$Gained); phyl_dist$Lost <- gsub("Mergellus albellus", "Mergus albellus", phyl_dist$Lost)
phyl_dist$Gained <- gsub("Morus bassanus", "Sula bassana", phyl_dist$Gained); phyl_dist$Lost <- gsub("Morus bassanus", "Sula bassana", phyl_dist$Lost)
phyl_dist$Gained <- gsub("Eudromias morinellus", "Charadrius morinellus", phyl_dist$Gained); phyl_dist$Lost <- gsub("Eudromias morinellus", "Charadrius morinellus", phyl_dist$Lost)
phyl_dist$Gained <- gsub("Regulus ignicapilla", "Regulus ignicapillus", phyl_dist$Gained); phyl_dist$Lost <- gsub("Regulus ignicapilla", "Regulus ignicapillus", phyl_dist$Lost)
phyl_dist$Gained <- gsub("Regulus ignicapilla", "Regulus ignicapillus", phyl_dist$Gained); phyl_dist$Lost <- gsub("Regulus ignicapilla", "Regulus ignicapillus", phyl_dist$Lost)
phyl_dist$Gained <- gsub("Parus caeruleus", "Cyanistes caeruleus", phyl_dist$Gained); phyl_dist$Lost <- gsub("Parus caeruleus", "Cyanistes caeruleus", phyl_dist$Lost)
phyl_dist$Gained <- gsub("Carduelis chloris", "Chloris chloris", phyl_dist$Gained); phyl_dist$Lost <- gsub("Carduelis chloris", "Chloris chloris", phyl_dist$Lost)
phyl_dist$Gained <- gsub("Parus ater", "Periparus ater", phyl_dist$Gained); phyl_dist$Lost <- gsub("Parus ater", "Periparus ater", phyl_dist$Lost)
phyl_dist$Gained <- gsub("Parus cristatus", "Lophophanes cristatus", phyl_dist$Gained); phyl_dist$Lost <- gsub("Parus cristatus", "Lophophanes cristatus", phyl_dist$Lost)

phyl_dist <- phyl_dist %>% bind_rows(rbind.data.frame(cbind.data.frame(Gained = "Loxia species", phyl_dist %>% filter(grepl("Loxia", Gained)) %>% group_by(Lost) %>% summarise(phyl_dist = mean(phyl_dist))),
                                                      cbind.data.frame(Lost = "Loxia species", phyl_dist %>% filter(grepl("Loxia", Lost)) %>% group_by(Gained) %>% summarise(phyl_dist = mean(phyl_dist)))))

phyl_dist %>% filter(!Gained %in% data$latin) %>% distinct(Gained) %>% print(n = 100)
data %>% filter(!latin %in% phyl_dist$Gained) %>% distinct(latin) %>% print(n = 100)

write_csv(phyl_dist, "//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Birds phylogeny/phylo_dist.csv")


###############################
######    Butterflies    ######
###############################

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
