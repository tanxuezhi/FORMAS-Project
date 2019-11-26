library(tidyverse)
library(data.table)
library(lmerTest)
library(MuMIn)
library(doFuture)
registerDoFuture()
plan(multiprocess, workers = 7)

### habitat data
frag_data <- read_csv("../Connectivity - Fragmentation/Fragmentation/Frag_indices_Allhab.csv")
sites_FIN <- read_csv(file = "../Data/Butterflies - Finland/Sites_FIN_ETRS89_landcover.csv")
sites_NL <- read_csv(file = "../Data/Butterflies - Netherlands/Sites_NL_ETRS89_landcover.csv")

### butterfly data
traits.butterflies <- as.tbl(read_csv("../Data/Traits/SpeciesTraits_WDV2014.csv"))
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

# Finland
countsFIN1 <- fread("../Data/Butterflies - Finland/FINLAND_Records_1999-2015.txt", sep = ";")
countsFIN2 <- fread("../Data/Butterflies - Finland/FINLAND_Records_2016.txt", sep = ";")[, -6]
countsFIN <- rbind(countsFIN1, countsFIN2)
colnames(countsFIN) <- gsub("Species_Faunaeur", "Species", colnames(countsFIN))
colnames(countsFIN) <- gsub("Individuals", "n", colnames(countsFIN))

countsFIN$Species <- gsub("Polygonia c-album", "Nymphalis c-album", countsFIN$Species)
countsFIN$Species <- gsub("Cyaniris semiargus", "Polyommatus semiargus", countsFIN$Species)
countsFIN$Species <- gsub("Leptidea juvernica", "Leptidea sinapis", countsFIN$Species)

countsFIN <- countsFIN %>%
  group_by(Site, Species, Year) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  complete(nesting(Site, Year), Species) %>%
  ungroup() %>%
  mutate(
    n = ifelse(is.na(n), 0, 1),
  ) %>%
  ungroup() %>%
  mutate(Site = paste0(Site, "_FIN"))

# NL
countsNL1 <- readRDS("../Data/Butterflies - Netherlands/AllSpecies_reg_gam_ind_20171206_algroutes.rds") %>%
  as.tbl() %>%
  dplyr::select(1, 2, 3, 4) %>%
  rename(n = regional_gam)
countsNL2 <- rio::import("../Data/Butterflies - Netherlands/MissingSpecies.xlsx") %>%
  as.tbl() %>%
  dplyr::select(1, 3, 4, 5) %>%
  rename(n = Ntot, SITE = Site) %>%
  mutate(n = ifelse(n == -1, 0, n))
countsNL <- bind_rows(countsNL1, countsNL2)
countsNL$SPECIES <- gsub("steppeparelmoervlinder       ", "steppevlekvlinder", countsNL$SPECIES)
countsNL$SPECIES <- gsub("kaasjeskruiddikkopje", "kaasjeskruiddikkop", countsNL$SPECIES)

countsNL <- countsNL %>%
  complete(nesting(SITE, YEAR), SPECIES) %>%
  mutate(n = ifelse(is.na(n) | n == 0, 0, 1)) %>%
  rename(Year = YEAR, Site = SITE, Species = SPECIES) %>%
  dplyr::select(Year, Site, Species, n)

countsNL <- left_join(
  countsNL,
  traits.butterflies %>% mutate(Species = casefold(Species))
) %>%
  dplyr::mutate(Species = Scientific_name, Site = paste0(Site, "_NL")) %>%
  select(2, 1, 3, 4)

# merge
pres_abs <- bind_rows(NL = countsNL, FIN = countsFIN, .id = "Country")

### select sites with at least 9 years
counts.fil <- pres_abs %>%
  group_by(Site) %>%
  summarise(n = n_distinct(Year)) %>%
  filter(n > 8)
pres_abs <- pres_abs %>% filter(Site %in% counts.fil$Site)

res <- c()
for (k in unique(pres_abs$Country)) {
  dat.country <- pres_abs %>% filter(Country == k)
  for (j in unique(dat.country$Species)) {
    for (i in unique(dat.country$Site)) {
      dat <- dat.country %>%
        filter(Site == i, Species == j)
      dat$period <- cut_interval(dat$Year, 3, labels = F)
      dat <- dat %>%
        group_by(period) %>%
        summarise(n = max(n))
      
      if (dat[1, 2] == 0 & dat[3, 2] == 1) {
        event <- "gain"
      }
      if (dat[1, 2] == 0 & dat[3, 2] == 0) {
        event <- "absence"
      }
      if (dat[1, 2] == 1 & dat[3, 2] == 1) {
        event <- "persistence"
      }
      if (dat[1, 2] == 1 & dat[3, 2] == 0) {
        event <- "loss"
      }
      
      res <- rbind.data.frame(res, cbind.data.frame(Country = k, Species = j, Site = i, event))
    }
  }
}

table(res[, c("Country", "event")])

write_csv(res, "../Long-term-event.csv")
res <- read_csv("../Long-term-event.csv")

#######
# analyses
#######
# fragmentation data
frag_data <- read_csv("../Connectivity - Fragmentation/Fragmentation/Frag_indices_Allhab.csv")

# merge with habitat and fragmentation
data <- res %>% left_join(butterfly_habitat)

data <- bind_rows(
  left_join(
    data %>% dplyr::filter(Habitat == "Forest"), frag_data,
    by = c("Site" = "Site", "Habitat" = "Habitat", "Country" = "country")
  ),
  left_join(
    data %>% dplyr::filter(Habitat == "Grassland" | Habitat == "Wetland") %>%
      mutate(Habitat = "Open"),
    frag_data,
    by = c("Site" = "Site", "Habitat" = "Habitat", "Country" = "country")
  ),
  left_join(
    data %>% dplyr::filter(is.na(Habitat)) %>%
      mutate(Habitat = "Generalist"),
    frag_data,
    by = c("Site" = "Site", "Habitat" = "Habitat", "Country" = "country")
  )
)

# add grid cell (50x50km)
sites <- bind_rows(
  sites_FIN %>% mutate(Site = paste0(Site, "_FIN")),
  sites_NL %>% mutate(Site = paste0(Site, "_NL"))
)

library(raster)

but.pts <- SpatialPoints(sites[, c("X", "Y")])
but.r <- raster(ext = extent(but.pts) * 1.1, resolution = 50000)
values(but.r) <- c(1:ncell(but.r))
gridCell50 <- raster::extract(but.r, but.pts)
sites <- bind_cols(sites, gridCell50 = gridCell50)

data <- data %>%
  left_join(sites) %>%
  as.tbl() %>%
  filter(!is.na(Scale))
data <- data %>%
  left_join(traits.butterflies %>% dplyr::select(-1, -3), by = c("Species" = "Scientific_name")) %>%
  filter(!is.na(STI))

butterflies.cti.presence <- as.tbl(fread("../Data/cti_butterflies_data.csv")) %>% dplyr::filter(type == "Presence")
data <- left_join(data, butterflies.cti.presence %>% group_by(Site) %>% summarise(cti = mean(cti)))
data <- data %>% mutate(STI_rel = cti - STI)

## big loop
library(DHARMa)


res.mods <- foreach(k = unique(data$Scale), .combine = rbind.data.frame) %dopar% {
  
  # colonisation
  dat.Col <- data %>%
    mutate(Colonisation = ifelse(event == "gain", 1, ifelse(event == "absence", 0, NA))) %>%
    filter(!is.na(Colonisation), Scale == k) %>%
    mutate_at(.vars = c(7, 8, 9, 10, 13:18, 20), .funs = function(x) stdize(x, prefix = F))
  
  
  m.col <- glmer(Colonisation ~ STI_rel * PC3 +
                   STI_rel * PC4 +
                   STI_rel * PC1 +
                   PC1 * PLAND * CLUMPY +
                   STI_rel * PLAND * CLUMPY +
                   Habitat + X * Y +
                   (1 | gridCell50 / Site) + (1 | Species),
                 family = binomial,
                 data = dat.Col,
                 nAGQ = 0, control = glmerControl(
                   optimizer = "nloptwrap",
                   optCtrl = list(maxfun = 1e10),
                   calc.derivs = FALSE
                 ), verbose = F
  )
  
  # optional : draw qqplots of residuals
  # simulationOutput_col <- simulateResiduals(fittedModel = m.col, n = 500, use.u = T)
  # png(file = paste0("../qqplot_col_scale_", k, ".png"), width = 4.5, height = 4, units = "in", res = 300)
  # plotQQunif(simulationOutput_col, testUniformity = T, testOutliers = T)
  # dev.off()
  
  # summary(m.col)
  
  
  # extinction
  dat.Ext <- data %>%
    mutate(Extinction = ifelse(event == "loss", 1, ifelse(event == "persistence", 0, NA)), CLUMPY = log(CLUMPY + 1.0001)) %>%
    filter(!is.na(Extinction), Scale == k) %>%
    mutate_at(.vars = c(7, 8, 9, 10, 13:18, 20), .funs = function(x) stdize(x, prefix = F))
  
  
  m.ext <- glmer(Extinction ~ STI_rel * PC3 +
                   STI_rel * PC4 +
                   STI_rel * PC1 +
                   PC1 * PLAND * CLUMPY +
                   STI_rel * PLAND * CLUMPY +
                   Habitat + X * Y +
                   (1 | gridCell50 / Site) + (1 | Species),
                 family = binomial,
                 data = dat.Ext,
                 nAGQ = 0, control = glmerControl(
                   optimizer = "nloptwrap",
                   optCtrl = list(maxfun = 1e10),
                   calc.derivs = FALSE
                 ), verbose = F
  )
  
  # optional : draw qqplots of residuals
  # simulationOutput_ext <- simulateResiduals(fittedModel = m.ext, n = 500, use.u = T)
  # png(file = paste0("../qqplot_ext_scale_", k, ".png"), width = 4.5, height = 4, units = "in", res = 300)
  # plotQQunif(simulationOutput_ext, testUniformity = T, testOutliers = T)
  # dev.off()
  
  # summary(m.ext)
  
  rbind.data.frame(
    cbind.data.frame(
      Scale = k,
      Process = "Colonisation",
      estimate = fixef(m.col),
      confint(m.col, method = "boot")[-c(1:3), ]
    ) %>% rownames_to_column("Variables"),
    cbind.data.frame(
      Scale = k,
      Process = "Extinction",
      estimate = fixef(m.ext),
      confint(m.ext, method = "boot")[-c(1:3), ]
    ) %>% rownames_to_column("Variables")
  )
}

write_csv(res.mods, "../rev_res_mods_251119.csv")


###################
##     PLOTS     ##
###################

##### Coefficients #####

res <- read_csv("../rev_res_mods_251119.csv")
names(res)[5:6] <- c("lwr", "upr")
res$Scale <- as.factor(res$Scale)
levels(res$Scale) <- paste(c(1, 3, 5, 10, 20, 30, 50), "km")

res$Variables <- gsub("poly_rescale\\(poly\\(PC1, 2\\), 2\\)2", "PC1^2", res$Variables)
res$Variables <- gsub("poly_rescale\\(poly\\(PC1, 2\\), 2\\)1", "PC1", res$Variables)
res$Variables <- gsub("STI_rel", "STI", res$Variables)
res$Variables <- gsub("PC1", "Spatial use", res$Variables)
res$Variables <- gsub("PC3", "Generation time", res$Variables)
res$Variables <- gsub("PC4", "Resource specialisation", res$Variables)
res$Variables <- gsub("PLAND", "Proportion of SNH", res$Variables)
res$Variables <- gsub("CLUMPY", "Aggregation of SNH", res$Variables)
res$Variables <- gsub(":", " × ", res$Variables)
res$Variables <- gsub("Proportion of SNH × Aggregation of SNH", "Proportion × Aggregation of SNH", res$Variables)
res$Variables <- gsub("X", "Longitude", res$Variables)
res$Variables <- gsub("Y", "Latitude", res$Variables)

res$Variables <- gsub("STI × Proportion × Aggregation of SNH", "STI × Proportion ×\nAggregation of SNH", res$Variables)
res$Variables <- gsub("Spatial use × Proportion × Aggregation of SNH", "Spatial use × Proportion ×\nAggregation of SNH", res$Variables)


res <- res %>% mutate(Significancy = ifelse(lwr < 0 & upr < 0 | lwr > 0 & upr > 0, "Significant", "Not significant"))

res$Scale <- factor(res$Scale, levels = rev(paste(c(1, 3, 5, 10, 20, 30, 50), "km")))

res$Variables <- factor(res$Variables, levels = c(
  "(Intercept)",
  "STI", "Spatial use", "Generation time", "Resource specialisation",
  "STI × Spatial use", "STI × Generation time", "STI × Resource specialisation",
  "HabitatGeneralist", "HabitatOpen",
  "Proportion of SNH", "Aggregation of SNH", "Proportion × Aggregation of SNH",
  "Spatial use × Proportion of SNH", "Spatial use × Aggregation of SNH",
  "Spatial use × Proportion ×\nAggregation of SNH",
  "STI × Proportion of SNH", "STI × Aggregation of SNH",
  "STI × Proportion ×\nAggregation of SNH",
  "Longitude", "Latitude", "Longitude × Latitude"
))


res <- res %>%
  arrange(Variables, Scale, Process) %>%
  mutate(Group = rep(c(
    rep("Various", 1),
    rep("Traits", 7),
    rep("Various", 2),
    rep("Fragmentation", 3),
    rep("Fragmentation × traits", 6),
    rep("Various", 3)
  ), each = 14))

res$Group <- factor(res$Group, levels = c("Various", "Traits", "Fragmentation", "Fragmentation × traits"))

### plot ###

ggplot(
  res %>% dplyr::filter(
    grepl("STI|Spatial use|Generation time|Resource specialisation|Proportion|Aggregation", Variables)
  ) %>%
    mutate(Variables = factor(Variables, levels = rev(levels(Variables)))),
  aes(x = Variables, y = estimate, color = Scale)
) +
  geom_hline(yintercept = 0, lty = 2, color = "black") +
  geom_point(position = position_dodge(.8)) + facet_grid(Group ~ Process, scales = "free", space = "fixed", switch = "y") +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(.8)) +
  coord_flip() +
  scale_y_continuous("Coefficients +/- 95% CI") +
  scale_x_discrete("") +
  scale_color_manual(values = colorRampPalette(c("black", "grey"))(7)) +
  theme_classic() +
  guides(color = guide_legend(reverse = TRUE)) +
  theme(strip.background = element_rect(colour=NA, fill=NA), 
        strip.placement = "outside",
        axis.text = element_text(color = "black"),
        strip.text.x = element_text(size = 11, face = "bold"),
        panel.border = element_blank()
  )


###############################
### drivers of colonisation ###
###     and extinction at   ###
###      20km scale         ###
###############################


### extract data at scale = 20km ###
dat <- sp_site_occTrend(left_join(dat.occ.trend, unique(sites[, c("Site", "gridCell50")])), .2, .8)
dat <- left_join(dat, sp) %>% left_join(sites)
dat <- stdize(dat, prefix = F, omit.cols = c("Trend", "Scale", "pred.then", "pred.now", "nYear"))

dat.Col <- dat %>%
  mutate(Colonisation = ifelse(Trend == "Colonisation", 1, ifelse(Trend == "No colonisation", 0, NA))) %>%
  dplyr::filter(Scale == 20000) %>%
  dplyr::select(-nYear)
dat.Col$Habitat <- as.factor(dat.Col$Habitat)
dat.Ext <- dat %>%
  mutate(Extinction = ifelse(Trend == "Extinction", 1, ifelse(Trend == "Persistence", 0, NA))) %>%
  dplyr::filter(Scale == 20000) %>%
  dplyr::select(-nYear)
dat.Ext$Habitat <- as.factor(dat.Ext$Habitat)

scaleList <- list(
  scale = attr(dat, "scaled:scale"),
  center = attr(dat, "scaled:center")
)


##################
## colonisation ##
##################
library(DHARMa)

m.col <- glmer(Colonisation ~ STI_rel * PC3 +
                 STI_rel * PC4 +
                 STI_rel * PC1 +
                 PC1 * PLAND * CLUMPY +
                 STI_rel * PLAND * CLUMPY +
                 Habitat + X * Y +
                 (1 | gridCell50 / Site) + (1 | Species),
               family = binomial,
               data = dat.Col,
               nAGQ = 0, control = glmerControl(
                 optimizer = "nloptwrap",
                 optCtrl = list(maxfun = 1e10),
                 calc.derivs = FALSE
               ), verbose = F
)
summary(m.col)

testDispersion(simulateResiduals(m.col))

##################
### Extinction ###
##################
m.ext <- glmer(Extinction ~ STI_rel * PC3 +
                 STI_rel * PC4 +
                 STI_rel * PC1 +
                 PC1 * PLAND * CLUMPY +
                 STI_rel * PLAND * CLUMPY +
                 Habitat + X * Y +
                 (1 | gridCell50 / Site) + (1 | Species),
               family = binomial,
               data = dat.Ext,
               nAGQ = 0, control = glmerControl(
                 optimizer = "nloptwrap",
                 optCtrl = list(maxfun = 1e10),
                 calc.derivs = FALSE
               ), verbose = F
)
summary(m.ext)
