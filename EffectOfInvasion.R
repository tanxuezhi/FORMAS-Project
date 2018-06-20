library(emmeans)
library(lmerTest)
library(visreg)
library(rotl)
library(FD)
library(mgcv)
library(patchwork)
library(foreach)
library(tidyverse)

##### load data #####
sites_FIN <- read_csv(file = "../Data/Butterflies - Finland/Sites_FIN_ETRS89_landcover.csv")
sites_NL <- as.tbl(left_join(rio:::import(file = "../Data/Butterflies - Netherlands/Sites_NL.xlsx"),
                             read_csv(file = "../Data/Butterflies - Netherlands/Sites_NL_ETRS89_landcover.csv")))

traits.butterflies <- as.tbl(read_csv("../Data/Butterflies - Netherlands/SpeciesTraits_WDV2014.csv"))
traits.butterflies[grep("Colias croceus", traits.butterflies$Scientific_name), "Scientific_name"] <- "Colias crocea"
traits.butterflies[grep("walbum", traits.butterflies$Scientific_name), "Scientific_name"] <- "Satyrium w-album"
traits.butterflies[grep("Neozephyrus quercus", traits.butterflies$Scientific_name), "Scientific_name"] <- "Favonius quercus"
traits.butterflies[grep("lycaon", traits.butterflies$Scientific_name), "Scientific_name"] <- "Hyponephele lycaon"
traits.butterflies[grep("Polygonia", traits.butterflies$Scientific_name), "Scientific_name"] <- "Nymphalis c-album"
traits.butterflies[grep("Inachis io", traits.butterflies$Scientific_name), "Scientific_name"] <- "Aglais io"
traits.butterflies[grep("tithonus", traits.butterflies$Scientific_name), "Scientific_name"] <- "Pyronia tithonus"
traits.butterflies[grep("alcon", traits.butterflies$Scientific_name), "Scientific_name"] <- "Phengaris alcon"

butterfly_habitat <- read_csv("../Data/Butterfly_biotopes.csv")
butterfly_habitat[grep("Neozephyrus quercus", butterfly_habitat$Species), "Species"] <- "Favonius quercus"

## Netherlands
data <- readRDS("../Data/Butterflies - Netherlands/AllSpecies_reg_gam_ind_20171206_algroutes.rds")

data %>% filter(YEAR > 1991, SPECIES == "landkaartje", regional_gam > 0) %>% mutate(SITE = as.factor(SITE)) %>%
  group_by(YEAR) %>% summarise(nSites = n()) %>%
  ggplot(aes(y = nSites, x = YEAR)) + geom_line() + geom_point() + theme(legend.position="none") + geom_smooth()

data %>% group_by(YEAR) %>% filter(YEAR > 1991) %>%
  summarise(nSites = n(), nSite_occ = length(unique(SITE[SPECIES == "landkaartje" & regional_gam > 0]))) %>%
  mutate(prop_occ = nSite_occ / nSites) %>%  
  ggplot(aes(y = prop_occ, x = YEAR)) + geom_line() + geom_point() + theme(legend.position="none") + geom_smooth()

## Finland
countsFIN1 <- read.table("../Data/Butterflies - Finland/FINLAND_Records_1999-2015.txt", sep = ";", h=T)
countsFIN2 <- read.table("../Data/Butterflies - Finland/FINLAND_Records_2016.txt", sep = ";", h=T)[,-6]
countsFIN <- rbind(countsFIN1, countsFIN2)
colnames(countsFIN) <- gsub ("Species_Faunaeur", "Species", colnames(countsFIN))

countsFIN$Species <- gsub("Polygonia c-album", "Nymphalis c-album", countsFIN$Species)
countsFIN$Species <- gsub("Cyaniris semiargus", "Polyommatus semiargus", countsFIN$Species)
countsFIN$Species <- gsub("Leptidea juvernica", "Leptidea sinapis", countsFIN$Species)

data2 <- countsFIN %>% group_by(Species, Year, Site) %>% summarise(Individuals = max(Individuals)) %>% left_join(sites_FIN)

data2 %>% group_by(Year) %>% 
  summarise(nSites = n(), nSite_occ = length(unique(Site[Species == "Araschnia levana" & Individuals > 0]))) %>%
  mutate(prop_occ = nSite_occ / nSites) %>%  
  ggplot(aes(y = prop_occ, x = Year)) + geom_line() + geom_point() + theme(legend.position="none") + geom_smooth()


# load phylogenetic data
taxa <- unique(data2$Species)
resolved_names <- tnrs_match_names(na.omit(taxa))

tree <- get_tree_ids("ott596629")

inspect(resolved_names, taxon_name = "aporia crataegi")
resolved_names <- update(resolved_names, taxon_name = "aporia crataegi",
                         new_row_number = 2)

inspect(resolved_names, taxon_name = "nymphalis c-album")
resolved_names <- update(resolved_names, taxon_name = "nymphalis c-album",
                         new_row_number = 3)

my_tree <- tol_induced_subtree(ott_ids = resolved_names$ott_id)
my_tree <- ape::compute.brlen(my_tree)

my_tree$tip.label <- gsub("_ott[0-9]*","",my_tree$tip.label)
my_tree$tip.label <- gsub("_"," ",my_tree$tip.label)
resolved_names$unique_name <- gsub(" \\(species in Holozoa\\)", "", resolved_names$unique_name)

my_tree$tip.label <- merge(data.frame(species = my_tree$tip.label), 
                           resolved_names, 
                           by.x = "species", by.y = "unique_name", sort=F)[,2]

my_tree$tip.label[my_tree$tip.label == "nymphalis c\\-album"] <- "Nymphalis c-album"


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

my_tree$tip.label <- firstup(my_tree$tip.label)
my_tree$node.label <- "NULL"

##### compute effect of A. levena invasion for all species #####

res <- c()
res$trends <- c()
res$abundance <- c()

for(i in unique(data2$Species)[!unique(data2$Species)%in% "Araschnia levana"]){
  print(i)
  
  dat.temp <- data2 %>% filter(Species == i) %>% ungroup()
  
  data2_invMap <- data2  %>% group_by(Site) %>% filter(Species == "Araschnia levana") %>% 
    select(-1) %>% summarise(inv.Year = min(Year)) %>%
    right_join(dat.temp, by = c("Site")) %>% filter(Species == i) %>%
    mutate(Invasion = ifelse(Year > inv.Year, "Yes", "No"),
           Year = Year - min(Year)) %>% 
    mutate(Invasion = factor(ifelse(is.na(Invasion), "No", Invasion), levels = c("No", "Yes")))
  
  
  if(length(unique(data2_invMap$Invasion)) > 1 & length(unique(data2_invMap$Individuals)) > 1 & 
     length(unique(data2_invMap$Invasion)) * length(unique(data2_invMap$Site)) < nrow(data2_invMap)){
    
    tryCatch({
      m.trends <- glmer(Individuals ~ Year*Invasion + scale(X)*scale(Y) + (Year|Site), family = "poisson", 
                        data = data2_invMap,
                        nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                                      optCtrl=list(maxfun=1e10),
                                                      calc.derivs = FALSE), verbose = F)
      
      trends <- pairs(emtrends(m.trends, ~ Invasion, "Year", transform = "response"))
      
      m.abundance <- glmer(Individuals ~ Invasion + scale(X)*scale(Y) + (Year|Site), family = "poisson", 
                           data = data2_invMap,
                           nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                                         optCtrl=list(maxfun=1e10),
                                                         calc.derivs = FALSE), verbose = F)
      
      abundance <- pairs(emmeans(m.abundance, ~ Invasion, transform = "response"))
      
      res$trends <- rbind.data.frame(res$trends, 
                                     cbind.data.frame(Species = i, 
                                                      nData = length(unique(data2_invMap$Invasion)) *
                                                        length(unique(data2_invMap$Site)), 
                                                      trends))
      res$abundance <- rbind.data.frame(res$abundance, 
                                        cbind.data.frame(Species = i, 
                                                         nData = length(unique(data2_invMap$Invasion)) * 
                                                           length(unique(data2_invMap$Site)), 
                                                         abundance))
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

## plot results
filter_threshold <- 0
  

bind_rows(Abundance = res$abundance %>% arrange(estimate), 
          "Long-term trend" = res$trends %>% arrange(estimate), .id = c("Type")) %>% 
  mutate(Species = factor(Species, unique(Species)), 
         significance = ifelse(p.value < .05, "yes", "no")) %>% 
  filter(nData > filter_threshold) %>%
  ggplot(aes(x = Species, y = estimate, color = significance)) + 
  facet_wrap( ~ Type, scales = "free_x") + 
  geom_point() + 
  geom_errorbar(aes(ymin = estimate - 2*SE, ymax = estimate + 2*SE), width = 0) + 
  coord_flip() +
  geom_hline(yintercept = 0, color = "dark grey", linetype = 2) +
  scale_color_manual(values = c("grey", "black")) + 
  scale_y_continuous(name = expression(paste("Effect of ", italic("Araschnia levana"), " invasion"))) + 
  theme_classic() +
  theme(legend.position="none")

plot((res$abundance  %>% filter(nData > filter_threshold))$estimate ~ (res$trends  %>% filter(nData > filter_threshold))$estimate, 
     ylab = expression(paste("Effect of ", italic("Araschnia levana"), " invasion on abundance")),
     xlab = expression(paste("Effect of ", italic("Araschnia levana"), " invasion on long-term trend")))

## merge with traits
res <- lapply(res, function(x){left_join(x,
                                         traits.butterflies %>% dplyr:::select(-Species),
                                         by = c("Species" = "Scientific_name"))})

traits_dat <- rbind(res$trends[,c(12:15)], traits.butterflies[traits.butterflies$Scientific_name == "Araschnia levana",-c(1:5)])
func_dist <- as.matrix(gowdis(traits_dat))
colnames(func_dist) <- c(res$trends[,1], "Araschnia levana"); rownames(func_dist) <- c(res$trends[,1], "Araschnia levana")
func_dist_to_inv <- func_dist[,"Araschnia levana"]

res <- lapply(res, function(x){
  cbind.data.frame(x, func_dist_to_inv = func_dist_to_inv[-length(func_dist_to_inv)])
})


foreach(j = 1:2, .combine = rbind) %:% foreach(i = 1:10, .combine = rbind) %do% {
  m.lmPoly <- lm(estimate ~ poly(func_dist_to_inv, i) ,
                 weight = nData^2, data = res[[j]] %>% filter(nData > filter_threshold, !is.na(STI)))
  cbind.data.frame(Type = names(res)[[j]], "Polynomial degree" = i, "AICc" = MuMIn::AICc(m.lmPoly))
}

par(mfrow=c(2,2))
m.ab.lmPoly.final <- lm(estimate ~ poly(func_dist_to_inv, 3) ,
                            weight = nData^2, data = res$abundance %>% filter(nData > filter_threshold, !is.na(STI)))
visreg(m.ab.lmPoly.final, scale = "response", ylim = c(-3, 3.5), rug = F, 
       ylab = expression(paste("Effect of ", italic("Araschnia levana"))),
       xlab = expression(paste("Ecological distance to ", italic("Araschnia levana"))),
       main = "Average abundance")
points(estimate ~ func_dist_to_inv, data = res$abundance %>% filter(nData > filter_threshold, !is.na(STI)))
mtext(adj = 0, paste("Polyn. LM: F = ", round(summary(m.ab.lmPoly.final)$fstatistic[1], 2), "; P = ", round(anova(m.ab.lmPoly.final)[1,5], 2), "; R2 = ", round(summary(m.ab.lmPoly.final)$adj.r.squared, 2)), cex = .6)

m.trends.lmPoly.final <- lm(estimate ~ poly(func_dist_to_inv, 10) ,
                     weight = nData^2, data = res$trends %>% filter(nData > filter_threshold, !is.na(STI)))
visreg(m.trends.lmPoly.final, scale = "response", ylim = c(-.5, 2), rug = F, 
       ylab = expression(paste("Effect of ", italic("Araschnia levana"))),
       xlab = expression(paste("Ecological distance to ", italic("Araschnia levana"))),
       main = "Long-term trend")
points(estimate ~ func_dist_to_inv, data = res$trends %>% filter(nData > filter_threshold, !is.na(STI)))
mtext(adj = 0, paste("Polyn. LM: F = ", round(summary(m.trends.lmPoly.final)$fstatistic[1], 2), "; P = ", round(anova(m.trends.lmPoly.final)[1,5], 2), "; R2 = ", round(summary(m.trends.lmPoly.final)$adj.r.squared, 2)), cex = .6)

m.ab.gam <- gam(estimate ~ s(func_dist_to_inv) ,
                        weight = nData^2, data = res$abundance %>% filter(nData > filter_threshold, !is.na(STI)))
visreg(m.ab.gam, scale = "response", ylim = c(-3, 3.5), rug = F, 
       ylab = expression(paste("Effect of ", italic("Araschnia levana"))),
       xlab = expression(paste("Ecological distance to ", italic("Araschnia levana"))))
points(estimate ~ func_dist_to_inv, data = res$abundance %>% filter(nData > filter_threshold, !is.na(STI)))
mtext(adj = 0, paste("GAM: F = ", round(summary(m.ab.gam)$s.table[1,3], 2), "; P = ", round(summary(m.ab.gam)$s.table[1,4], 2), "; R2 = ", round(summary(m.ab.gam)$r.sq, 2)), cex = .6)

m.trends.gam <- gam(estimate ~ s(func_dist_to_inv) ,
                            weight = nData^2, data = res$trends %>% filter(nData > filter_threshold, !is.na(STI)))
visreg(m.trends.gam, scale = "response", ylim = c(-.5, 2), rug = F, 
       ylab = expression(paste("Effect of ", italic("Araschnia levana"))),
       xlab = expression(paste("Ecological distance to ", italic("Araschnia levana"))))
points(estimate ~ func_dist_to_inv, data = res$trends %>% filter(nData > filter_threshold, !is.na(STI)))
mtext(adj = 0, paste("GAM: F = ", round(summary(m.trends.gam)$s.table[1,3], 2), "; P = ", round(summary(m.trends.gam)$s.table[1,4], 2), "; R2 = ", round(summary(m.trends.gam)$r.sq, 2)), cex = .6)
