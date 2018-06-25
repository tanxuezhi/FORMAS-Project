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

countsFIN <- countsFIN %>% 
  group_by(Site, Species, Year) %>% 
  summarise(Individuals = max(Individuals)) %>% 
  ungroup() %>%
  group_by(Site) %>%
  complete(Year, nesting(Site, Species)) %>%
  mutate(Individuals = ifelse(is.na(Individuals), 0, Individuals)) %>% left_join(sites_FIN)


countsFIN %>% group_by(Year) %>% 
  summarise(nSites = n(), nSite_occ = length(unique(Site[Species == "Araschnia levana" & Individuals > 0]))) %>%
  mutate(prop_occ = nSite_occ / nSites) %>%  
  ggplot(aes(y = prop_occ, x = Year)) + geom_line() + geom_point() + theme(legend.position="none") + geom_smooth()



library(raster)
FIN_map <- getData(country="Finland", level = 0)
FIN_map <- spTransform(FIN_map, CRS("+init=epsg:3035"))

r <- raster(ext = extent(FIN_map), resolution = 5000)

FIN_map_res <- rasterize(FIN_map, r)
plot(FIN_map_res)

inv.Year <- data2 %>% group_by(Site) %>% summarise(inv.year = min(Year[Species == "Araschnia levana"]), 
                                                   X= first(X), Y = first(Y)) %>%
  with(as.data.frame(.))

inv.cells <- cbind(extract(FIN_map_res, inv.Year[,3:4], cellnumbers = T), inv.Year)[,c(1,4)] %>%
  group_by(cells) %>% summarise_all(first) %>% with(as.data.frame(.))
inv.cells <- rbind(inv.cells,
                   cbind.data.frame(cells = 1:ncell(FIN_map_res), inv.year = NA) %>% filter(!cells %in% inv.cells[,1]))

FIN_map_res.inv <- setValues(FIN_map_res, inv.cells %>% arrange(cells) %>% with(as.vector(.$inv.year)))
FIN_map_res.inv <- aggregate(FIN_map_res.inv, 10, min)


FIN_map_res.sites <- aggregate(FIN_map_res, 10, mean)
FIN_map_res.sites[extract(FIN_map_res.sites, inv.Year[,3:4], cellnumbers = T)[,1]] <- 2
library(rasterVis)
levelplot(FIN_map_res.sites, margin = F) + levelplot(FIN_map_res.inv, margin = F)



# load phylogenetic data
taxa <- unique(data2$Species)
resolved_names <- tnrs_match_names(na.omit(taxa))

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

save(my_tree, file =  "../Butterflies_tree.RData")

##### compute effect of A. levena invasion for all species #####

dat.temp <- data2 %>% group_by(Site) %>% summarise(inv.year = min(Year[Species == "Araschnia levana"]), 
                                                   year.before = length(unique(Year[Year < inv.year])),
                                                   year.after = length(unique(Year[Year > inv.year]))) %>% 
  mutate(inv.year = ifelse(inv.year == Inf, NA, inv.year)) %>% filter(year.before > 2, year.after > 2) %>%
  left_join(data2) %>% mutate(BeforeAfter = as.factor(ifelse(Year <= inv.year, "Before", "After")))


dat.temp2 <- dat.temp %>% filter(Site == 56, Species == "Aglais urticae") %>% 
  mutate(growthRate = log(Individuals+1) - log(lag(Individuals)+1))


mGR <- lm(growthRate ~ BeforeAfter, 
          data = dat.temp2)
plot(growthRate ~ BeforeAfter, data = dat.temp2)

mAb <- glm(Individuals ~ BeforeAfter, family = "poisson", 
           data = dat.temp2)
plot(Individuals ~ BeforeAfter, data = dat.temp2)


res <- c()
res$trends <- c()
res$abundance <- c()

for(i in unique(data2$Species)[!unique(data2$Species) %in% "Araschnia levana"]){
  print(i)
  
  dat.temp <- data2 %>% filter(Species == i) %>% ungroup()
  
  data2_invMap <- data2 %>% group_by(Site) %>% filter(Species == "Araschnia levana") %>% 
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



####################

dat.temp <- countsFIN  %>% filter(Species == "Araschnia levana", Individuals > 0) %>% 
  group_by(Site)%>% summarise(inv.year = min(Year),inv.severity = sum(Individuals)) %>% 
  mutate(inv.year = ifelse(inv.year == Inf, NA, inv.year)) %>% 
  right_join(countsFIN) %>% mutate(BeforeAfter = as.factor(ifelse(Year < inv.year, "Before", "After")))

for(i in unique(dat.temp$Species[!unique(dat.temp$Species) %in% "Araschnia levana"])){
  
  dat.temp2 <- dat.temp %>% filter(Species == i) %>% ungroup()
  
  dat.temp2 <- dat.temp2 %>% filter(BeforeAfter == "Before") %>% group_by(Site) %>% summarise(Before = n()) %>% 
    left_join(dat.temp2 %>% filter(BeforeAfter == "After") %>% group_by(Site) %>% summarise(After = n())) %>%
    filter(Before > 1, After > 1) %>% left_join(dat.temp2)
  
  dat.temp4 <- c()
  for(j in unique(dat.temp2$Site)){
    
    impact.site <- dat.temp2 %>% filter(Site == j)
    
    dists <- spDistsN1(dat.temp %>% filter(is.na(inv.year), Species == i) %>% 
                         group_by(Site) %>% summarise(X = unique(X), Y = unique(Y),
                                                      min.year = min(Year), max.year = max(Year)) %>% 
                         filter(min.year < unique(impact.site$inv.year), 
                                max.year > unique(impact.site$inv.year))
                       %>% dplyr::select(X,Y) %>% with(as.matrix(.)),
                       impact.site %>% group_by(Site) %>% 
                         summarise_all(first) %>% dplyr::select(X,Y) %>% with(as.matrix(.)))
    
    control.site <- dat.temp %>% filter(is.na(inv.year), Species == i) %>% 
      group_by(Site) %>% summarise_all(first) %>% slice(which.min(dists)) %>% dplyr::select(Site) %>% 
      left_join(dat.temp  %>% filter(Species == i), by = "Site") %>% 
      mutate(BeforeAfter = ifelse(Year < unique(impact.site$inv.year), "Before", "After"))
    
    
    dat.temp3 <- cbind.data.frame(bind_rows(Impact = impact.site, Control = control.site, .id = "ControlImpact"), 
                                  Pair = paste0(j ,"_", unique(control.site$Site)), Distance = min(dists))
    dat.temp4 <- rbind(dat.temp4, dat.temp3)
  }
  
  dat.temp4 <- dat.temp4 %>% mutate(Year = Year - min(Year))
  
  m <- glmer(Individuals ~ BeforeAfter*ControlImpact + (1|Pair/Site) + (1|Site:Year), 
             weights = 1/(scale(Distance,center = F)),
             family = "poisson", data = dat.temp4,
             nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                           optCtrl=list(maxfun=1e10),
                                           calc.derivs = FALSE), verbose = F)
}

m1 <- glmer(Individuals ~ BeforeAfter*ControlImpact + (1|Pair/Site) + (1|Year), 
            weights = 1/(scale(Distance,center = F)),
            family = "poisson", data = dat.temp4,
            nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                          optCtrl=list(maxfun=1e10),
                                          calc.derivs = FALSE), verbose = F)

m2 <- glmer(Individuals ~ BeforeAfter*ControlImpact + (1|Pair/Site) + (1|Site:Year), 
            weights = 1/(scale(Distance,center = F)),
            family = "poisson", data = dat.temp4,
            nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                          optCtrl=list(maxfun=1e10),
                                          calc.derivs = FALSE), verbose = F)

m3 <- glmer(Individuals ~ BeforeAfter*ControlImpact + (Year|Pair/Site) + (1|Site:Year), 
            weights = 1/(scale(Distance,center = F)),
            family = "poisson", data = dat.temp4,
            nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                          optCtrl=list(maxfun=1e10),
                                          calc.derivs = FALSE), verbose = F)

m4 <- glmer(Individuals ~ BeforeAfter*ControlImpact + (Year|Pair/Site) + (1|Year), 
            weights = 1/(scale(Distance,center = F)),
            family = "poisson", data = dat.temp4,
            nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                          optCtrl=list(maxfun=1e10),
                                          calc.derivs = FALSE), verbose = F)

m5 <- glmer(Individuals ~ BeforeAfter*ControlImpact + (1|Pair/Site), 
            weights = 1/(scale(Distance,center = F)),
            family = "poisson", data = dat.temp4,
            nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                          optCtrl=list(maxfun=1e10),
                                          calc.derivs = FALSE), verbose = F)

m6 <- glmer(Individuals ~ BeforeAfter*ControlImpact + (Year|Pair/Site), 
            weights = 1/(scale(Distance,center = F)),
            family = "poisson", data = dat.temp4,
            nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                          optCtrl=list(maxfun=1e10),
                                          calc.derivs = FALSE), verbose = F)

m7 <- glm(Individuals ~ BeforeAfter*ControlImpact, 
          weights = 1/(scale(Distance,center = F)),
          family = "poisson", data = dat.temp4)

m8 <- glmer(Individuals ~ BeforeAfter*ControlImpact + (1|Pair) + (1|Year), 
            weights = 1/(scale(Distance,center = F)),
            family = "poisson", data = dat.temp4,
            nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                          optCtrl=list(maxfun=1e10),
                                          calc.derivs = FALSE), verbose = F)

m9 <- glmer(Individuals ~ BeforeAfter*ControlImpact + (Year|Pair) + (1|Year), 
            weights = 1/(scale(Distance,center = F)),
            family = "poisson", data = dat.temp4,
            nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                          optCtrl=list(maxfun=1e10),
                                          calc.derivs = FALSE), verbose = F)

m10 <- glmer(Individuals ~ BeforeAfter*ControlImpact + (Year|Pair) + (1|Pair:Year), 
             weights = 1/(scale(Distance,center = F)),
             family = "poisson", data = dat.temp4,
             nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                           optCtrl=list(maxfun=1e10),
                                           calc.derivs = FALSE), verbose = F)

m11 <- glmer(Individuals ~ BeforeAfter*ControlImpact + (1|Pair) + (1|Pair:Year), 
             weights = 1/(scale(Distance,center = F)),
             family = "poisson", data = dat.temp4,
             nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                           optCtrl=list(maxfun=1e10),
                                           calc.derivs = FALSE), verbose = F)




MuMIn::model.sel(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11)
emmeans(m2, ~BeforeAfter*ControlImpact, transform = "response") %>% with(as.data.frame(.)) %>% 
  mutate(BeforeAfter = factor(BeforeAfter, levels = c("Before", "After"))) %>%
  ggplot(aes(x = BeforeAfter, color = ControlImpact, y = rate)) + geom_point(position = position_dodge(width = .3)) + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0, position = position_dodge(width = .3)) + 
  ggtitle(expression(paste("Effect of ", italic("Araschnia levana"), " invasion on", italic(" Aglais urticae")))) + 
  scale_y_continuous("Average abundance") + 
  scale_color_discrete(name = "", labels = c("Control site", "Invaded site")) + 
  scale_x_discrete("") + 
  theme_classic()
