library(data.table)
library(tidyverse)
library(raster)
library(sp)
library(landscapemetrics)
library(rgeos)
library(randomForest)
library(cooccur)
library(lmerTest)

#####################
## background info ##
#####################

expanding_species <- countsFIN %>% group_by(Year) %>% mutate(nSites_tot = length(unique(Site))) %>% group_by(Species, Year) %>%
  summarise(nSites = length(unique(Site)), nSites_tot = unique(nSites_tot)) %>% group_by(Species) %>%
  do(trend = glm(cbind(nSites,nSites_tot) ~ Year, family = binomial, data = .)) %>% 
  broom::tidy(trend) %>% filter(term == "Year") %>% filter(p.value < .05, estimate > 0) %>% arrange(desc(estimate)) %>%
  pull(Species)


#################
### load data ###
#################

##### butterfly data ######
# traits to merge names
traits.butterflies <- as.tbl(read_csv("../Data/Traits/SpeciesTraits_WDV2014.csv"))
traits.butterflies[grep("Colias croceus", traits.butterflies$Scientific_name), "Scientific_name"] <- "Colias crocea"
traits.butterflies[grep("walbum", traits.butterflies$Scientific_name), "Scientific_name"] <- "Satyrium w-album"
traits.butterflies[grep("Neozephyrus quercus", traits.butterflies$Scientific_name), "Scientific_name"] <- "Favonius quercus"
traits.butterflies[grep("lycaon", traits.butterflies$Scientific_name), "Scientific_name"] <- "Hyponephele lycaon"
traits.butterflies[grep("Polygonia", traits.butterflies$Scientific_name), "Scientific_name"] <- "Nymphalis c-album"
traits.butterflies[grep("Inachis io", traits.butterflies$Scientific_name), "Scientific_name"] <- "Aglais io"
traits.butterflies[grep("tithonus", traits.butterflies$Scientific_name), "Scientific_name"] <- "Pyronia tithonus"
traits.butterflies[grep("alcon", traits.butterflies$Scientific_name), "Scientific_name"] <- "Phengaris alcon"


# FIN
countsFIN1 <- fread("../Data/Butterflies - Finland/FINLAND_Records_1999-2015.txt", sep = ";", h=T)
countsFIN2 <- fread("../Data/Butterflies - Finland/FINLAND_Records_2016.txt", sep = ";", h=T)[,-6]
countsFIN <- rbind(countsFIN1, countsFIN2)
colnames(countsFIN) <- gsub ("Species_Faunaeur", "Species", colnames(countsFIN))
colnames(countsFIN) <- gsub ("Individuals", "n", colnames(countsFIN))

countsFIN <- countsFIN %>% 
  group_by(Site, Species, Year) %>% 
  summarise(n = max(n)) %>% ungroup %>% as.tbl # keep max no. individuals

countsFIN$Site <- paste0(countsFIN$Site, "_FIN")

countsFIN <- countsFIN %>% filter(n > 0)

# NL
countsNL1 <- readRDS("../Data/Butterflies - Netherlands/AllSpecies_reg_gam_ind_20171206_algroutes.rds") %>% as.tbl %>%
  filter(prop_pheno_sampled > 0.5) %>%
  dplyr::select(1,2,3,4) %>% rename(n = regional_gam)
countsNL2 <- rio::import("../Data/Butterflies - Netherlands/MissingSpecies.xlsx") %>% as.tbl %>%
  dplyr::select(1,3,4,5) %>% rename(n = Ntot, SITE = Site) %>% mutate(n = ifelse(n == -1, 0 , n))
countsNL <- bind_rows(countsNL1, countsNL2)

countsNL <- countsNL %>% rename(Site = SITE, Species = SPECIES, Year = YEAR)
countsNL$Site <- paste0(countsNL$Site, "_NL")
countsNL <- countsNL %>% filter(n > 0)
countsNL$Species <- gsub("       ", "", countsNL$Species)


countsNL <- left_join(countsNL,
                       traits.butterflies %>% mutate(Species = casefold(Species)) %>% dplyr:::select(1,2)) %>%
  dplyr::select(-2) %>%
  rename(Species = Scientific_name)

countsNL[is.na(countsNL$Species), "Species"] <- "Melitaea aurelia"


##### birds data ######
# SWE
countsSWE1 <- rio::import("../Data/Birds - Sweden/public_totalstandard_Ia.xlsx")
countsSWE2 <- rio::import("../Data/Birds - Sweden/public_totalstandard_Ib.xlsx")
countsSWE3 <- rio::import("../Data/Birds - Sweden/public_totalstandard_IIa.xlsx")
countsSWE4 <- rio::import("../Data/Birds - Sweden/public_totalstandard_IIb.xlsx")
countsSWE <- as.tbl(bind_rows(countsSWE1, countsSWE2, countsSWE3, countsSWE4))

rm(list = c("countsSWE1", "countsSWE2", "countsSWE3", "countsSWE4"))

countsSWE$art <- as.numeric(countsSWE$art)
countsSWE <- countsSWE %>% dplyr:::filter(art %in% 1:645)
countsSWE <- countsSWE[,c(2,4,5,22,23)]

countsSWE$n <- countsSWE$pkind + countsSWE$lind
countsSWE <- countsSWE[,-c(4,5)]

colnames(countsSWE) <- gsub ("karta", "Site", colnames(countsSWE))
colnames(countsSWE) <- gsub ("art", "Species", colnames(countsSWE))
colnames(countsSWE) <- gsub ("yr", "Year", colnames(countsSWE))

countsSWE <- countsSWE %>% filter(n > 0)


### merge ###
data_counts <- bind_rows("NL" = countsNL %>% filter(Year > 1991 & Year < 2017) %>% mutate(Site = as.character(Site)),
                         "FIN" = countsFIN %>% mutate(Site = as.character(Site)),
                         "SWE" = countsSWE %>% mutate(Species = as.character(Species)),
                         .id = "Data")

#############################
##### extract site data #####
#############################

## load climate data ##
pre <- stack(lapply(list.files("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/CRU_data\\pre", full.names = T), function(x)(stack(x))))
tmp <- stack(lapply(list.files("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/CRU_data\\tmp", full.names = T), function(x)(stack(x))))
tmn <- stack(lapply(list.files("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/CRU_data\\tmn", full.names = T), function(x)(stack(x))))
tmx <- stack(lapply(list.files("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/CRU_data\\tmx", full.names = T), function(x)(stack(x))))

## site data ##
sitesFIN <- read_csv(file = "../Data/Butterflies - Finland/Sites_FIN_ETRS89_landcover.csv")
sitesFIN$Site <- paste0(sitesFIN$Site, "_FIN")
sitesNL <- read_csv(file = "../Data/Butterflies - Netherlands/Sites_NL_ETRS89_landcover.csv")
sitesNL$Site <- paste0(sitesNL$Site, "_NL")
sitesSWE <- read_csv(file = "../Data/Birds - Sweden/Sites_SWE_ETRS89_landcover.csv")

sites <- bind_rows("NL" = sitesNL,
                   "FIN" = sitesFIN,
                   "SWE" = sitesSWE,
                   .id = "Data")


sites_3035 <- SpatialPoints(sites[,c("X", "Y")], proj4string=CRS("+init=epsg:3035"))
sites_wgs84 <- spTransform(sites_3035, crs(pre))

## extract temperature data
sites_temp <- c()
for(z in unique(data_counts$Year)){
  cat(z);cat("\n")
  pre_site_temp <- extract(mean(subset(pre, grep(z, names(pre)))), sites_wgs84)
  tmp_site_temp <- extract(mean(subset(tmp, grep(z, names(tmp)))), sites_wgs84)
  tmn_site_temp <- extract(mean(subset(tmn, grep(z, names(tmn)))), sites_wgs84)
  tmx_site_temp <- extract(mean(subset(tmx, grep(z, names(tmx)))), sites_wgs84)
  
  sites_temp <- rbind.data.frame(sites_temp,
                                 cbind.data.frame(sites[,1], Year = z,
                                                  pre = pre_site_temp, 
                                                  tmp = tmp_site_temp, 
                                                  tmn = tmn_site_temp, 
                                                  tmx = tmx_site_temp)
  )
}
sites_temp$Site <- rep(sites$Site, length(unique(data_counts$Year)))

###########################
### compute annual SDMs ###
###########################

dist_sites <- as.matrix(dist(sites_3035@coords))
colnames(dist_sites) <- sites$Site
rownames(dist_sites) <- sites$Site


suitability <- c()
for(j in unique(data_counts$Data)){
  dat_temp1 <- data_counts %>% filter(Data == j)
  
  
  for(i in unique(data_counts$Year)){
    
    dat_temp2 <- dat_temp1 %>% filter(Year == i)
    
    for(k in unique(dat_temp2$Species)){
      
      sp_obs <- dat_temp2 %>% filter(Species == k)
      sp_No_obs <- dat_temp2 %>% filter(!Site %in% sp_obs$Site) %>% mutate(n = 0, Species = unique(sp_obs$Species)) %>% 
        unique()
      
      sp_dat <- bind_rows(sp_obs, sp_No_obs) %>% left_join(sites_temp) %>% na.omit()
      
      m <- NULL
      tryCatch({
        m <- mahal(sp_dat %>% filter(n > 0) %>% dplyr::select(6:9))
        pred <- predict(m, sp_dat %>% dplyr::select(6:9))
        pred <- 1/(-pred + 1 + 1)
        
        suitability <-  rbind.data.frame(suitability,
                                         cbind.data.frame(Data = j, 
                                                          Year = i, 
                                                          Species = k, 
                                                          Site = sp_dat$Site, 
                                                          suitability = pred)
        )
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
      if(is.null(m)){
        
        centroid = colMeans(sp_dat %>% filter(n > 0) %>% dplyr::select(6:9))

        pred <- c(rep(1, nrow(sp_dat %>% filter(n > 0))),
                      1/(apply(sp_dat %>% filter(n == 0) %>% dplyr::select(6:9), 1, function(x)dist(rbind(x, centroid)))+1))
        suitability <-  rbind.data.frame(suitability,
                                         cbind.data.frame(Data = j, 
                                                          Year = i, 
                                                          Species = k, 
                                                          Site = sp_dat$Site, 
                                                          suitability = pred)
        )
      }
      
    }
  }
}


#######################################
### compute non-random associations ###
#######################################

assoc.Sp <- c()
for(n in unique(data_counts$Data)){
  data_counts_temp <- data_counts[data_counts$Data == n,]
  for(m in c(min(data_counts_temp$Year) : max(data_counts_temp$Year))){
    countsFIN_comm <- data_counts_temp %>% ungroup() %>% filter(Year == m) %>% mutate(n = 1) %>% 
      spread(Species, n, fill = 0) %>% as.data.frame %>% 
      mutate(id = paste0(Site,"_",Year)) %>% dplyr::select(-1,-2) %>% column_to_rownames("id")
    countsFIN_comm <- countsFIN_comm[,colSums(countsFIN_comm) > (0.01 * nrow(countsFIN_comm)) & 
                                       colSums(countsFIN_comm) < (0.99 * nrow(countsFIN_comm))]
    
    assoc.Sp_temp <- cooccur(t(countsFIN_comm), type = "spp_site", thresh = F, spp_names = T,
                             prob = "comb")
    assoc.Sp_temp <- assoc.Sp_temp$results
    assoc.Sp_temp <- assoc.Sp_temp %>% mutate(association = ifelse(p_gt < 0.05, "Aggregated", 
                                                                   ifelse(p_lt < 0.05, "Segregated", "Random")),
                                              Data = n,
                                              Year = m)
    assoc.Sp <- rbind.data.frame(assoc.Sp, assoc.Sp_temp)
  }
}
assoc.Sp <- as.tbl(assoc.Sp)

assoc.Sp %>% group_by(Data, association, Year) %>% summarise(n = n())

## with species-site mask ##
# 
# sp_site_suitability2 <- countsFIN %>% dplyr::select(-4) %>% ungroup %>% 
#   complete(Species, Year, Site) %>%
#   left_join(sp_site_suitability) %>% mutate(suitable = ifelse(suitable == 1, 1, 0)) %>% 
#   filter(Species %in% colnames(countsFIN_comm)) %>%
#   group_by(Species, Site) %>% spread(Species, suitable) %>% as.data.frame %>% 
#   mutate(id = paste0(Site,"_",Year)) %>% dplyr::select(-1,-2)%>% 
#   filter(id %in% rownames(countsFIN_comm))  %>% column_to_rownames("id")
# 
# assoc.Sp <- cooccur(t(countsFIN_comm), type = "spp_site", thresh = F, spp_names = T,
#                     site_mask = t(sp_site_suitability2),
#                     prob = "comb")
# plot(assoc.Sp)
# pair.profile(assoc.Sp)
# obs.v.exp(assoc.Sp)


## filter interactions ##

res_assoc <- c()
for(l in 1:nrow(assoc.Sp2)){
  
  diff_env <- NULL
  diff_dist <- NULL
  
  sp1 <- c(assoc.Sp2 %>% slice(l) %>% pull(sp1_name) %>% as.character)
  sp2 <- c(assoc.Sp2 %>% slice(l) %>% pull(sp2_name) %>% as.character)
  
  
  if(assoc.Sp2 %>% slice(l) %>% pull(association) == "Segregated"){
    
    assoc_cond_10 <- countsFIN %>% ungroup() %>% filter(Species %in% sp1, 
                                                        Year == assoc.Sp2 %>% slice(l) %>% pull(Year)) %>% 
      filter(!Species %in% sp2) %>% 
      dplyr::select(Site, Year) %>% distinct %>%
      left_join(sites_clim_pca)
    
    assoc_cond_01 <- countsFIN %>% ungroup() %>% filter(Species %in% sp2, 
                                                        Year == assoc.Sp2 %>% slice(l) %>% pull(Year)) %>% 
      filter(!Species %in% sp1) %>% 
      dplyr::select(Site, Year) %>% distinct %>%
      left_join(sites_clim_pca)
    
    assoc_cond_merged <- bind_rows("10" = assoc_cond_10, "00" = assoc_cond_01, .id = "assoc")
    
  } else if(assoc.Sp2 %>% slice(l) %>% pull(association) == "Aggregated"){
    
    assoc_cond_11 <- countsFIN %>% ungroup() %>% filter(Species %in% c(sp1, sp2), 
                                                        Year == assoc.Sp2 %>% slice(l) %>% pull(Year)) %>% 
      group_by(Site, Year) %>% summarise(n = length(unique(Species))) %>% filter(n > 1) %>%
      dplyr::select(Site, Year) %>% distinct %>%
      left_join(sites_clim_pca)
    
    assoc_cond_00 <- countsFIN %>% ungroup() %>% filter(!Species %in% c(sp1, sp2), 
                                                        !Site %in% assoc_cond_11$Site, 
                                                        Year == assoc.Sp2 %>% slice(l) %>% pull(Year)) %>% 
      dplyr::select(Site, Year) %>% distinct %>%
      left_join(sites_clim_pca)
    
    assoc_cond_merged <- bind_rows("00" = assoc_cond_00, "11" = assoc_cond_11, .id = "assoc")
    
  }
  
  association <- assoc.Sp2 %>% slice(l) %>% pull(association)
  
  if(association != "Random"){
    diff_env <- summary(manova(cbind(Axis1, Axis2) ~ assoc, data = assoc_cond_merged))$stat[1,6]
    diff_dist <- summary(manova(cbind(X, Y) ~ assoc, data = assoc_cond_merged))$stat[1,6]
  }
  
  
  if(association == "Random"){
    explanation <- "Random"
  } else if(association == "Aggregated"){
    if(diff_dist < 0.05){
      if(diff_env < 0.05){
        explanation <- "Dispersal_limitation_or_environmental_filtering"
      } else {
        explanation <- "Dispersal_limitation"
      }
    } else {
      if(diff_env < 0.05){
        explanation <- "Environmental_filtering"
      } else {
        explanation <- "Positive_species_interaction"
      }
    }
  } else if(association == "Segregated"){
    if(diff_dist < 0.05){
      if(diff_env < 0.05){
        explanation <- "Dispersal_limitation_or_environmental_filtering"
      } else {
        explanation <- "Dispersal_limitation"
      }
    } else {
      if(diff_env < 0.05){
        explanation <- "Environmental_filtering"
      } else {
        explanation <- "Negative_species_interaction"
      }
    }
  }
  
  res_assoc <- rbind.data.frame(res_assoc, 
                                cbind.data.frame(sp1 = sp1, sp2 = sp2,
                                                 Year = assoc.Sp2[l,"Year"],
                                                 association, explanation)
  )
}

res_assoc <- as.tbl(res_assoc)
table(res_assoc$explanation)


res_assoc <- res_assoc %>% rename(sp1_name = sp1, sp2_name = sp2)


res_assoc[res_assoc$explanation == "Negative_species_interaction",]
res_assoc[res_assoc$explanation == "Positive_species_interaction",]

############################
## change in cooccurrence ##
############################
library(igraph)

dat_used <- assoc.Sp

# optional: if filtered associations are used
dat_used <- dat_used %>% mutate(association = ifelse(explanation == "Positive_species_interaction", "Aggregated",
                                                     ifelse(explanation == "Negative_species_interaction", 
                                                            "Segregated", "Random")))

net_stat <- c()
for(l in unique(dat_used$Data)){
  dat2 <- dat_used %>% filter(Data == l) 
  for(i in unique(dat2$Year)){
    for(k in c("Aggregated", "Segregated")){
      dat1 <- dat2 %>% filter(association == k,  Year == i) %>% mutate(sp1_name = as.character(sp1_name),
                                                                       sp2_name = as.character(sp2_name)) %>%
        dplyr::select(sp1_name,sp2_name)
      
      g1 <- graph_from_data_frame(dat1, directed = F,
                                  vertices = unique(c(dat1$sp1_name, dat1$sp2_name)))
      
      net_stat <- rbind.data.frame(net_stat,
                                   cbind.data.frame(Data = l,
                                                    Year = i,
                                                    association = k,
                                                    sum = sum(degree(g1)), # sum links
                                                    Link_density = mean(degree(g1)), # links per species
                                                    Connectance = sum(degree(g1)) / (length(E(g1))^2), # connectance
                                                    Modularity = modularity(cluster_walktrap(g1)) # modularity
                                   )
      )
    }
  }
}

net_stat %>% gather(-1,-2,-3,-4, key = "type", value = "value") %>% group_by(Data, type, association) %>% 
  mutate(max = first(value)) %>%
  ungroup() %>%
  mutate(value_prop = value/max) %>% filter(association == "Aggregated") %>%
  ggplot(aes(y = value, x = Year, color = Data))  + 
  geom_abline(slope = 0, intercept = 1) +
  geom_point() + geom_path() +
  facet_wrap(Data~type, scales = "free") + 
  geom_smooth() + theme_bw()


distri_deg <- c()
for(l in unique(dat_used$Data)){
  dat2 <- dat_used %>% filter(Data == l) 
  for(i in unique(dat2$Year)){
    for(k in c("Aggregated", "Segregated")){
      dat1 <- dat2 %>% filter(association == k,  Year == i) %>% mutate(sp1_name = as.character(sp1_name),
                                                                       sp2_name = as.character(sp2_name)) %>%
        dplyr::select(sp1_name,sp2_name)
      
      g1 <- graph_from_data_frame(dat1, directed = F,
                                  vertices = unique(c(dat1$sp1_name, dat1$sp2_name)))
      
      distri_deg <- rbind.data.frame(distri_deg,
                                     cbind.data.frame(Data = l,
                                                      Year = i,
                                                      association = k,
                                                      degree_distribution = degree_distribution(g1) # modularity
                                     )
      )
    }
  }
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

distri_deg  %>% group_by(Data, Year, association) %>% 
  summarise(degree_distribution = Mode(degree_distribution)) %>%
  filter(association == "Aggregated") %>%
  ggplot(aes(y = degree_distribution, x = Year, color = Data)) + geom_point() + geom_path() +
  theme_bw()


#######################
## clean environment ##
#######################

gdata::keep(sites,
            sites_temp,
            dist_sites,
            data_counts,
            suitability,
            assoc.Sp,
            net_stat, distri_deg,
            sure=T)
