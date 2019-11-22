library(data.table)
library(tidyverse)
library(raster)
library(sp)
# library(landscapemetrics)
library(rgeos)
# library(netassoc)
library(dismo)
# library(lmerTest)
library(igraph)
# library(future.apply)
library(viridis)
library(cooccur)


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
list.spSWE <- read.csv("../Data/Birds - Sweden/Species_list.csv")

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
countsSWE <- countsSWE %>% left_join(list.spSWE %>% dplyr::select(1,3), by = c("Species" = "art")) %>% mutate(Species = latin) %>% dplyr::select(-latin)

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
  
  for(k in unique(dat_temp1$Species)){
    
    sp_obs <- dat_temp1 %>% filter(Species == k)
    sp_No_obs <- dat_temp1 %>% filter(!paste(Site, Year) %in% paste(sp_obs$Site, sp_obs$Year)) %>% 
      mutate(n = 0, Species = unique(sp_obs$Species)) %>% 
      unique()
    
    sp_dat <- bind_rows(sp_obs, sp_No_obs) %>% left_join(sites_temp)
    for(z in 6:ncol(sp_dat)){
      sp_dat[is.na(sp_dat[,z]), z] <- mean(c(sp_dat[,z])[[1]], na.rm = TRUE)
    }
    
    m <- NULL
    tryCatch({
      m <- mahal(sp_dat %>% filter(n > 0) %>% dplyr::select(6:9))
      pred <- predict(m, sp_dat %>% dplyr::select(6:9))
      pred <- 1/(-pred + 1 + 1)
      
      suitability <-  rbind.data.frame(suitability,
                                       cbind.data.frame(sp_dat[,c(1,2,3,5)], 
                                                        suitability = pred)
      )
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    if(is.null(m)){
      
      centroid = colMeans(sp_dat %>% filter(n > 0) %>% dplyr::select(6:9))
      
      pred <- c(rep(1, nrow(sp_dat %>% filter(n > 0))),
                1/(apply(sp_dat %>% filter(n == 0) %>% dplyr::select(6:9), 1, function(x)dist(rbind(x, centroid)))+1))
      suitability <-  rbind.data.frame(suitability,
                                       cbind.data.frame(sp_dat[,c(1,2,3,5)], 
                                                        suitability = pred)
      )
    }
  }
}

dist_max <- c()
for(j in unique(data_counts$Data)){
  dat_temp1 <- data_counts %>% filter(Data == j)
  
  for(k in unique(dat_temp1$Species)){
    dat_temp2 <- dat_temp1 %>% filter(Species == k)
    
    for(l in unique(dat_temp2$Year)){
      sp_obs <- dat_temp2 %>% filter(Year == l)
      
      dist_max <- rbind.data.frame(dist_max,
                                   cbind.data.frame(
                                     Data = j, Species = k, Year = l,
                                     dist = ifelse(nrow(sp_obs) > 1, 
                                                   max(apply(dist_sites[rownames(dist_sites) %in% 
                                                                          sp_obs$Site,colnames(dist_sites) %in% sp_obs$Site],
                                                             1, function(x)min(x[x != 0]))), 0)
                                   )
      )
      
    }
  }
}
dist_max_sum <- dist_max %>% group_by(Data, Species) %>% summarise(dist_max = max(dist))

### clip by distance to nearest site with species

suitability_clipped <- c()
for(j in unique(data_counts$Data)){
  dat_temp1 <- data_counts %>% filter(Data == j)
  
  for(k in unique(dat_temp1$Species)){
    dat_temp2 <- dat_temp1 %>% filter(Species == k)
    
    for(l in unique(dat_temp2$Year)){
      
      if(suitability %>% filter(Year == l, Species == k, Data == j) %>% nrow > 0){
        
        sp_obs <- data_counts %>% filter(Year == l, Species == k, Data == j)
        
        dat_temp3 <- suitability %>% filter(Year == l, Species == k, Data == j) %>% filter(Site %in% rownames(dist_sites))
        
        if(nrow(sp_obs) > 1){
          dat_temp3 <- left_join(
            dat_temp3, 
            data.frame(dist = apply(dist_sites[rownames(dist_sites) %in% dat_temp3$Site,
                                               colnames(dist_sites) %in% sp_obs$Site],
                                    1, min)) %>% rownames_to_column("Site"),
            by = "Site"
          )
          dat_temp3[dat_temp3$dist > dist_max_sum %>% filter(Data == j, Species == k) %>% pull(dist_max),"suitability"] <- 0
          
        } else if(nrow(sp_obs) == 1){
          dat_temp3 <- left_join(
            dat_temp3, 
            cbind.data.frame(Site = rownames(dist_sites)[rownames(dist_sites) %in% dat_temp3$Site],
                             dist = as.vector(dist_sites[rownames(dist_sites) %in% dat_temp3$Site,
                                                         colnames(dist_sites) %in% sp_obs$Site])),
            by = "Site"
          )
          dat_temp3[dat_temp3$dist > dist_max_sum %>% filter(Data == j, Species == k) %>% pull(dist_max),"suitability"] <- 0
          
        }else if(nrow(sp_obs) == 0){
          dat_temp3[,"suitability"] <- 0
        }
        
        suitability_clipped <- rbind.data.frame(suitability_clipped, dat_temp3[,-6])
        
      }
    }
  }
}

#################################
####        Filter data       ###
#################################

data_counts %>% group_by(Data) %>% 
  summarise(n = length(unique(Year)), min = min(Year), max = max(Year))

data_counts %>% group_by(Data, Year) %>% 
  summarise(n = length(unique(Site))) %>% 
  ggplot(aes(x = Year, y = n)) + 
  geom_point() + 
  geom_smooth() +
  facet_grid(~Data)

## keep sites surveyed more than 70% of years
sel_sites <- data_counts %>%  group_by(Data, Site) %>% summarise(n = length(unique(Year)), min = min(Year), max = max(Year)) %>%
  group_by(Data) %>%
  mutate(keep = ifelse(unique(Data) == "NL" & n > 17, "Yes", 
                       ifelse(unique(Data) == "FIN" & n > 12, "Yes", 
                              ifelse(unique(Data) == "SWE" & n > 14, "Yes","No")))) %>%
  filter(keep %in% "Yes", !Data %in% "FIN") %>% pull(Site)

data_counts_filtered <- data_counts %>% filter(Site %in% sel_sites, Site %in% sites$Site)

data_counts_filtered <- data_counts %>% filter(Year >= 2007)
data_counts_filtered %>% group_by(Data, Year) %>% 
  summarise(n = length(unique(Site))) %>% 
  ggplot(aes(x = Year, y = n)) + 
  geom_point() + 
  geom_smooth() +
  facet_grid(~Data)


# check temporal change in latitude
data_counts_filtered %>% left_join(sites) %>% 
  group_by(Data, Year) %>% summarise(Y = mean(Y, na.rm = T)) %>%
  ggplot(aes(y = Y, x = Year)) + geom_point() + facet_wrap(~Data, scale = "free")

data_counts %>% left_join(sites) %>% filter(Data != "FIN") %>%
  group_by(Data, Year) %>% summarise(Y = mean(Y, na.rm = T)) %>%
  ggplot(aes(y = Y, x = Year)) + geom_point() + facet_wrap(~Data, scale = "free")

#######################################
### compute non-random associations ###
#######################################

assoc.Sp <- vector(mode="list", 3)
aa = 0 
for(n in unique(data_counts_filtered$Data)){
  data_counts_temp <- data_counts_filtered[data_counts_filtered$Data == n,]
  suitability_temp <- suitability_clipped[suitability_clipped$Data == n,]
  
  aa = aa + 1
  bb = 0
  
  for(m in c(min(data_counts_temp$Year) : max(data_counts_temp$Year))){
    
    cat(paste("Data:", n))
    cat("\n")
    cat(paste("Year:", m))
    cat("\n")
    cat("\n")
    
    
    bb = bb + 1
    
    obs_comm <- data_counts_temp %>% ungroup() %>% filter(Year == m) %>% 
      mutate(n = ifelse(n > 0, 1, 0)) %>%
      spread(Species, n, fill = 0) %>% as.data.frame %>% 
      mutate(id = paste0(Site,"_",Year)) %>% dplyr::select(-Site,-Year,-Data) %>% column_to_rownames("id")
    
    null_comm <- suitability_temp %>% ungroup() %>% 
      filter(Year == m, Site %in% substr(rownames(obs_comm),1,nchar(rownames(obs_comm))-5),
             Species %in% colnames(obs_comm)) %>% mutate(suitability = ifelse(suitability > .8, 1, 0)) %>%
      spread(Species, suitability, fill = 0) %>% as.data.frame %>% 
      mutate(id = paste0(Site,"_",Year)) %>% dplyr::select(-Site,-Year,-Data) %>% column_to_rownames("id")
    null_comm <- null_comm[order(rownames(null_comm)),order(colnames(null_comm))]
    
    
    net_test <- cooccur(
      mat = t(obs_comm),
      site_mask = t(null_comm),
      type = "spp_site",
      spp_names = T,
      prob = "comb",
      thresh = F
    )
    
    net_test_dat <- net_test$results %>% 
      mutate(type = ifelse(p_gt < .05, "Positive", ifelse(p_lt < .05, "Negative", "Random")))
    
    net_test_neg <- net_test_dat %>% filter(type == "Negative") %>% 
      mutate(weight = abs(scale(obs_cooccur) - scale(exp_cooccur)))
    net_test_neg <- graph_from_data_frame(net_test_neg[,c(10:11,13)], directed = F, vertices = levels(net_test_dat$sp2_name))
    
    
    net_test_pos <- net_test_dat %>% filter(type == "Positive") %>% 
      mutate(weight = abs(scale(obs_cooccur) - scale(exp_cooccur)))
    net_test_pos <- graph_from_data_frame(net_test_pos[,c(10:11,13)], directed = F, vertices = levels(net_test_dat$sp2_name))
    
    res <- list(net_test_dat, net_test_pos, net_test_neg)
    names(res) <- c("Data", "Positive", "Negative")
    
    assoc.Sp[[aa]][[bb]] <- res
    
  }
  names(assoc.Sp[[aa]]) <- c(min(data_counts_temp$Year) : max(data_counts_temp$Year))
}
names(assoc.Sp) <- unique(data_counts_filtered$Data)


############################
## change in cooccurrence ##
############################

net_stat <- c()
for(l in 1:length(assoc.Sp)){
  dat2 <- assoc.Sp[[l]] 
  for(i in 1:length(dat2)){
    for(k in c("Positive", "Negative")){
      
      g1 <- dat2[[i]][[k]]
      
      net_stat <- rbind.data.frame(net_stat,
                                   cbind.data.frame(Data = names(assoc.Sp)[l],
                                                    Year = names(dat2)[i],
                                                    Association = k,
                                                    sum = sum(degree(g1)), # sum links
                                                    Link_density = mean(degree(g1)), # links per species
                                                    Connectance = sum(degree(g1)) / (length(V(g1))^2), # connectance
                                                    Modularity = modularity(cluster_fast_greedy(g1)) # modularity
                                   )
      )
    }
  }
}

distri_deg <- c()
for(l in 1:length(assoc.Sp)){
  dat2 <- assoc.Sp[[l]] 
  for(i in 1:length(dat2)){
    for(k in c("Positive", "Negative")){
      
      g1 <- dat2[[i]][[k]]
      
      distri_deg <- rbind.data.frame(distri_deg,
                                     cbind.data.frame(Data = names(assoc.Sp)[l],
                                                      Year = names(dat2)[i],
                                                      Association = k,
                                                      degree_distribution = median(degree_distribution(g1))
                                     )
      )
    }
  }
}

net_stat <- net_stat %>% mutate(Modularity = ifelse(Modularity < 0, 0 , Modularity))

net_stat %>% gather(-1,-2,-3,-4, key = "type", value = "value") %>% 
  ungroup() %>% mutate(Year = as.numeric(as.character(Year))) %>% 
  filter(type %in% c("Connectance", "Modularity"), !(Data == "FIN" & Association == "Negative")) %>% 
  ggplot(aes(y = value, x = Year, color = Data, fill = Data)) + 
  geom_point() + 
  geom_path() +
  facet_grid(type~Association, scale = "free_y") + 
  scale_color_viridis(discrete=TRUE) + scale_fill_viridis(discrete=TRUE)

net_stat %>% gather(-1,-2,-3,-4, key = "type", value = "value") %>% mutate(Year = as.numeric(as.character(Year))) %>%
  filter(type %in% c("Connectance", "Modularity"), !(Data == "FIN" & Association == "Negative")) %>% 
  group_by(Data, Association, type) %>% do(trend = lm(value ~ Year, data = .)) %>% 
  broom::tidy(trend) %>%
  filter(term == "Year")

##############################
##     species richness     ##
##############################

data_counts_filtered %>% group_by(Data, Year) %>% summarise(n = length(unique(Species))) %>%
  mutate(Year = as.numeric(as.character(Year))) %>%
  ggplot(aes(y = n, x = Year)) + 
  geom_point() + 
  geom_smooth() +
  facet_grid(Data~., scale = "free_y") + 
  theme_classic()

data_counts_filtered %>% group_by(Data, Year) %>% summarise(n = length(unique(Species))) %>%
  group_by(Data) %>% 
  do(trend = lm(n ~ Year, data = .)) %>% broom::tidy(trend)

data_counts_filtered %>% group_by(Data, Year) %>% summarise(n = length(unique(Site)))%>%
  group_by(Data) %>% 
  do(trend = lm(n ~ Year, data = .)) %>% broom::tidy(trend)

############################
##     beta diversity     ##
############################

library(betapart)

beta_div <- c()
for(n in unique(data_counts_filtered$Data)){
  data_counts_temp <- data_counts_filtered[data_counts_filtered$Data == n,]
  
  for(m in c(min(data_counts_temp$Year) : max(data_counts_temp$Year))){
    
    obs_comm <- data_counts_temp %>% ungroup() %>% filter(Year == m) %>% 
      mutate(n = ifelse(n > 0, 1, 0)) %>%
      spread(Species, n, fill = 0) %>% as.data.frame %>% 
      mutate(id = paste0(Site,"_",Year)) %>% dplyr::select(-Site,-Year,-Data) %>% column_to_rownames("id")
    
    beta_div <- rbind.data.frame(beta_div, 
                                 cbind.data.frame(Data = n, Year = m, t(unlist(beta.multi(obs_comm)))
                                 )
    )
    
  }
}

beta_div %>% 
  ggplot(aes(x = Year, y = beta.SOR)) + 
  geom_point() + 
  geom_smooth() +
  facet_wrap(~Data, scale = "free") +
  theme_classic()

beta_div %>% group_by(Data) %>% do(trend = lm(beta.SOR ~ Year, data = .)) %>% broom::tidy(trend)


beta_div %>% left_join(net_stat %>% mutate(Year = as.numeric(as.character(Year)))) %>%
  ggplot(aes(x = beta.SOR, y = Connectance, color = Association)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_grid(~Data) +
  theme_classic()


beta_div %>% left_join(net_stat %>% mutate(Year = as.numeric(as.character(Year)))) %>%
  group_by(Data, Association) %>% do(cor = lm(Connectance ~ beta.SOR, data = .)) %>% broom::tidy(cor)


data_counts_filtered %>% group_by(Data, Year) %>% summarise(n = length(unique(Site))) %>%
  ggplot(aes(x = Year, y = n)) + 
  geom_point() + 
  geom_smooth() +
  facet_grid(~Data) +
  theme_classic()

data_counts_filtered %>% group_by(Data, Year) %>% summarise(n = length(unique(Site))) %>%
  left_join(net_stat %>% mutate(Year = as.numeric(as.character(Year)))) %>%
  ggplot(aes(x = n, y = Connectance, color = Association)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_grid(~Data) +
  theme_classic()


#############################
##   results per species   ##
#############################

sp_stat <- c()
for(l in 1:length(assoc.Sp)){
  dat2 <- assoc.Sp[[l]] 
  for(i in 1:length(dat2)){
    for(k in c("Positive", "Negative")){
      
      g1 <- dat2[[i]][[k]]
      
      n.sp.poss <- as.data.frame(
        table(dat2[[i]][["Data"]]$sp1_name)+
          table(dat2[[i]][["Data"]]$sp2_name)
      ) %>% rename(Species = Var1)
      
      
      sp.as <- unlist(strsplit(attr(E(g1), "vnames"), "\\|"))
      sp.all <- attr(V(g1), "names")
      
      if(length(sp.as) > 0){
        sp_stat <- rbind.data.frame(
          sp_stat,
          cbind.data.frame(Data = names(assoc.Sp)[l],
                           Year = names(dat2)[i],
                           Association = k,
                           rbind.data.frame(as.data.frame(table(sp.as)) %>% 
                                              rename(Species = sp.as, N = Freq),
                                            cbind.data.frame(Species = setdiff(sp.all, sp.as), N = 0)),
                           n.sp = length(sp.all)
          ) %>% left_join(n.sp.poss)
        )
      } else {
        sp_stat <- rbind.data.frame(
          sp_stat,
          cbind.data.frame(Data = names(assoc.Sp)[l],
                           Year = names(dat2)[i],
                           Association = k,
                           cbind.data.frame(Species = sp.all, N = 0),
                           n.sp = length(sp.all)
          ) %>% left_join(n.sp.poss)
        )
      }
    }
  }
}


sp_btw <- c()
for(l in 1:length(assoc.Sp)){
  dat2 <- assoc.Sp[[l]] 
  for(i in 1:length(dat2)){
    for(k in c("Positive", "Negative")){
      
      g1 <- dat2[[i]][[k]]
      
      sp.as <- unique(unlist(strsplit(attr(E(g1), "vnames"), "\\|")))
      sp.all <- attr(V(g1), "names")
      
      sp_btw <- rbind.data.frame(
        sp_btw,
        cbind.data.frame(Data = names(assoc.Sp)[l],
                         Year = names(dat2)[i],
                         Association = k,
                         betweenness(g1) %>% as.data.frame %>% rownames_to_column("Species") %>% rename(bet = ".")
        )
      )
    }
  }
}

sp_stat <- sp_stat %>% group_by(Data, Species, Association) %>% mutate(n.nonr = length(N[N > 0])) %>% filter()

sp_stat2 <- sp_stat %>% 
  dplyr::select(-Freq) %>% spread(key = Association, value = N) %>%
  mutate(Total = Positive + Negative, Year = as.numeric(as.character(Year))) %>% 
  gather(-1,-2,-3,-4, key = Association, value = N)

sp_stat_trends_tot <- c()
for(i in unique(sp_stat2$Data)){
  for(j in unique(sp_stat2$Association)){
    sp_stat_trends_temp <- sp_stat2 %>% filter(Data == i, Association == j)
    for(k in unique(sp_stat_trends_temp$Species)){
      m <- NULL
      tryCatch(
        m <- glm(cbind(N, (n.sp-1)) ~ Year, family = binomial, data = sp_stat_trends_temp %>% filter(Species == k)),
        warning = function(w) {
          print("warning")
        })
      if(!is.null(m)){
        if(dim(summary(m)$coefficient)[1]>1){
            sp_stat_trends_tot <- rbind.data.frame(sp_stat_trends_tot,
                                                   cbind.data.frame(
                                                     Data = i,
                                                     Association = j,
                                                     Species = k,
                                                     trend.assoc = summary(m)$coefficient[2,1],
                                                     se.assoc = summary(m)$coefficient[2,2]
                                                   )
            )
        }
      }
    }
  }
}


sp_btw_trends_tot <- sp_btw %>% group_by(Data, Species, Association) %>% 
  mutate(Year = as.numeric(as.character(Year))) %>% 
  do(trend = lm(bet ~ Year, data = .)) %>%
  broom::tidy(trend) %>% filter(term == "Year") %>% dplyr::select(-term)

# occupancy trends
trends_per_species <- data_counts %>% group_by(Data, Year) %>% summarise(n.sites = length(unique(Site))) %>% 
  right_join(data_counts) %>% group_by(Species, Year, Data) %>% 
  summarise(occ.sites = n(), tot.sites = unique(n.sites)) %>% 
  group_by(Species, Data) %>% do(trend = glm(cbind(occ.sites,tot.sites) ~ Year, family = binomial, data = .)) %>%
  broom::tidy(trend) %>% filter(term == "Year") %>% dplyr::select(1,2,4,5) %>% rename(trend.occ = estimate,
                                                                                      se.occ = std.error) %>%
  group_by(Data) %>% mutate(se.occ = scale(se.occ))


# trend assoc ~ occupancy trend
# total
sp_stat_trends_tot %>% 
  left_join(trends_per_species) %>% filter(!(Data == "FIN" & Association == "Negative")) %>%
  ggplot(aes(x = trend.occ, y = trend.assoc)) + 
  geom_point(pch = 1) + 
  geom_smooth(method = "lm", aes(weight = 1/se.assoc)) + 
  geom_hline(yintercept = 0, linetype = 2, color = "grey") + 
  geom_vline(xintercept = 0, linetype = 2, color = "grey") + 
  facet_grid(Data~Association, scale = "free") + 
  theme_classic() +
  scale_y_continuous("Trend in associations") + 
  scale_x_continuous("Occupancy trend")

sp_stat_trends_tot %>% left_join(trends_per_species) %>% 
  filter(!(Data == "FIN" & Association == "Negative")) %>% 
  group_by(Data, Association) %>%
  do(fit = lm(trend.assoc ~ trend.occ, weight = 1/se.assoc, data = .)) %>% 
  broom::tidy(fit) %>% filter(term == "trend.occ")


trends_per_species %>% group_by(Data) %>% mutate(tot = length(unique(Species))) %>% 
  group_by(Data, trend = ifelse(trend < 0, "negative", "positive")) %>% 
  summarise(n = n(), prop = n / unique(tot)) 

sp_stat_trends_tot %>% 
  left_join(trends_per_species) %>% filter(std.error < 1) %>% 
  group_by(Data) %>% 
  ggplot(aes(x = Data, y = abs(estimate), fill = Association)) + geom_boxplot()


par(mfrow = c(3,3))
for(i in unique(sp_stat_trends_tot$Data)){
  for(j in unique(sp_stat_trends_tot$Association)){
    dat <- sp_stat_trends_tot %>% 
      left_join(trends_per_species) %>% 
      filter(Data == i, Association == j)
    m <- lm(trend.assoc ~ trend.occ, weights = 1/se.assoc, data = dat)
    # hist(dat$estimate)
    plot(m, 2, main = paste(i , j))
  }
}

# ratio
sp_stat_trends_rat %>%
  left_join(trends_per_species) %>% 
  ggplot(aes(x = trend, y = estimate)) + 
  geom_point() + 
  geom_smooth(method = "lm", aes(weight = 1/(std.error+0.0001))) + 
  facet_grid(Data~., scale = "free")

sp_stat_trends_rat %>% left_join(trends_per_species) %>% 
  group_by(Data) %>%
  do(fit = lm(estimate ~ trend, weights = 1/(std.error+0.0001), data = .)) %>% broom::tidy(fit) %>% filter(term == "trend")




#######################
## clean environment ##
#######################

gdata::keep(sites,
            sites_temp,
            dist_sites,
            data_counts,
            data_counts_filtered,
            beta_div,
            sel_sites,
            suitability,
            suitability_clipped,
            dist_max,
            dist_max_sum,
            assoc.Sp,
            net_stat, 
            distri_deg,
            sp_stat,
            sp_stat_trends_tot,
            trends_per_species,
            sure=T)
