library(tidyverse)
# library(broom)
# library(ape)
library(vegan)
library(foreach)
library(broom)
library(cowplot)

##### load data #####

data <- read_csv("../Data/butterfly_data_invasion.csv") %>% filter(!is.na(Species))


##### identify expanding species #####
sp.occupancy.trend <- data %>% group_by(Year) %>% summarise(n.sites = length(unique(Site))) %>% 
  right_join(data) %>% group_by(Species, Year, country = ifelse(grepl("_NL", Site), "NL", "FIN")) %>% summarise(Prop.sites = n()/unique(n.sites)) %>% 
  group_by(Species, country) %>% do(trend = gam(Prop.sites ~ Year, data = .)) %>% tidy(trend) %>% filter(term == "Year")

sp.occupancy.trend %>% group_by(country) %>% filter(estimate == max(estimate))


p1 <- data %>% group_by(Year) %>% summarise(n.sites = length(unique(Site))) %>% 
  right_join(data) %>% group_by(Species, Year, country = ifelse(grepl("_NL", Site), "NL", "FIN")) %>% summarise(Prop.sites = n()/unique(n.sites)) %>% filter(country == "NL", Species == "Pararge aegeria") %>% 
  ggplot(aes(x =  Year, y = Prop.sites)) + geom_point() + geom_smooth()

p2 <- data %>% group_by(Year) %>% summarise(n.sites = length(unique(Site))) %>% 
  right_join(data) %>% group_by(Species, Year, country = ifelse(grepl("_NL", Site), "NL", "FIN")) %>% summarise(Prop.sites = n()/unique(n.sites)) %>% filter(country == "FIN", Species == "Araschnia levana") %>% 
  ggplot(aes(x =  Year, y = Prop.sites)) + geom_point() + geom_smooth()

plot_grid(p1, p2, nrow = 2, labels = c("NL - Pararge aegeria", "FIN - Araschnia levana"))


decay <- function(dat){
  com <- dat %>% select(Species,Year,n) %>% spread(., Species, n, fill = 0)
  d <- vegdist(com[,-c(1:4)], method = "bray")
  dist.com <- cbind.data.frame(dist = as.matrix(d)[,1], Year = com[,1])
  m.test <- glm(dist ~ Year, data = dist.com %>% mutate(Year = Year - min(Year)), family = binomial, 
                method = "detect_separation")
  if(!m.test$separation){
    m <- glm(dist ~ Year, data = dist.com %>% mutate(Year = Year - min(Year)), family = binomial)
    return(coef(m)[[2]])
  }
}

decay_AL <- data %>% filter(Site %in% (data %>% group_by(country = ifelse(grepl("_NL", Site), "NL", "FIN")) %>% 
                             filter(Species == "Araschnia levana", country == "FIN"))$Site) %>% 
  group_by(Site) %>% 
  filter(!length(unique(Year))<2) %>% do(decay = decay(.)) %>% tidy(decay)

decay_NO.AL <- data %>% filter(Site %in% (data %>% group_by(country = ifelse(grepl("_NL", Site), "NL", "FIN")) %>% 
                                         filter(Species != "Araschnia levana", country == "FIN"))$Site) %>% 
  group_by(Site) %>% 
  filter(!length(unique(Year))<2) %>% do(decay = decay(.)) %>% tidy(decay)




ggplot(sp.occupancy.trend, aes(x =  Year, y = Prop.sites, color = Species)) + geom_point() + 
  geom_smooth(se = F, method  ="gam", formula = y ~ s(x, k = 3)) +
  theme(legend.position = "none")


##### similarity decay #####

## all data ###

res_all <- c()
for(i in unique(data$Site)){
  
  data.temp.site <- data %>% filter(Site == i, n > 0)
  
  if(length(unique(data.temp.site$Species))>4 & length(unique(data.temp.site$Year))>4){
    
    com <- data.temp.site %>% spread(., Species, n, fill = 0)
    
    d <- vegdist(com[,-c(1:4)], method = "bray")
    
    dist.com <- cbind.data.frame(dist = as.matrix(d)[,1], Year = com[,1])
    m.test <- glm(dist ~ Year, data = dist.com %>% mutate(Year = Year - min(Year)), family = binomial, 
                  method = "detect_separation")
    
    if(!m.test$separation){
      
      m <- glm(dist ~ Year, data = dist.com %>% mutate(Year = Year - min(Year)), family = binomial)
      
      res_all <- rbind.data.frame(res_all, cbind.data.frame(Site = i, n.Year.site = length(unique(data.temp.site$Year)),
                                                            coef = coef(m)[[2]]))
    }
  }
}

## species jackknife ##

res_sp <- c()
for(i in unique(data$Site)){
  
  data.temp.site <- data %>% filter(Site == i, n > 0)
  
  if(length(unique(data.temp.site$Species))>4 & length(unique(data.temp.site$Year))>4){
    
    for(j in unique(data.temp.site$Species)) {
      
      data.temp.sp <- data.temp.site %>% filter(Species != j)
      
      com <- data.temp.sp %>% spread(., Species, n, fill = 0)
      
      d <- vegdist(com[,-c(1:4)], method = "bray")
      
      dist.com <- cbind.data.frame(dist = as.matrix(d)[,1], Year = com[,1])
      m.test <- glm(dist ~ Year, data = dist.com %>% mutate(Year = Year - min(Year)), family = binomial, 
                    method = "detect_separation")
      
      if(!m.test$separation){
        
        m.test <- glm(dist ~ Year, data = dist.com %>% mutate(Year = Year - min(Year)), family = binomial)
        
        res_sp <- rbind.data.frame(res_sp, cbind.data.frame(Species = j, Site = i, 
                                                            n.Year.sp = nrow(data.temp.site) - nrow(data.temp.sp),
                                                            rel.ab = sum((data.temp.site %>% filter(Species == j))$n)/
                                                              sum((data.temp.site)$n),
                                                            coef = coef(m)[[2]]))
      }
    }
  }
}

res_sp.site <- res_sp %>% left_join(res_all, by = "Site") %>% mutate(sp_effect = (coef.y - coef.x) / coef.y) %>% 
  rename(coef_sp_excluded = coef.x, coef_all_sp = coef.y)

res_sp_mean <- res_sp.site %>% group_by(Species) %>% summarise(mean_effect = mean(sp_effect, na.rm = T),
                                                               n.site = n(),
                                                               n.Year.sp = sum(n.Year.sp),
                                                               rel.ab = mean(rel.ab)) %>%
  ungroup() %>% mutate(Species = forcats::fct_reorder(Species, mean_effect))

ggplot(res_sp_mean, aes(y = Species, x = mean_effect)) + geom_point()


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


res_sp_mean_traits <- res_sp_mean %>% left_join(traits.butterflies, by = c("Species" = "Scientific_name")) %>% 
  select(-Species.y, -FF_code)

res_sp_mean_traits %>% gather(., "Trait", "Val", -1,-2,-3, -4, -5) %>% 
  ggplot(aes(x = Val, y = mean_effect, color = n.Year.sp)) +
  geom_point() + geom_smooth(method ="lm") + 
  facet_wrap(~Trait, scale = "free")

library(mgcv)
mm <- gam(mean_effect ~ s(PC1)+s(PC2)+s(PC3)+s(PC4), weight = I(1/n.site),  data = res_sp_mean_traits)

mm <- lm(mean_effect ~ PC1 + PC2 + PC3 + PC4, weight = I(n.Year.sp), data = res_sp_mean_traits)


