library(unmarked)
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

table(res$event)

write_csv(res, "../Long-term-event.csv")

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
res.mods <- foreach(k = unique(butterflies$Scale), .combine = rbind.data.frame) %dopar% {

  # colonisation
  dat.Col <- data %>%
    mutate(Colonisation = ifelse(event == "gain", 1, ifelse(event == "absence", 0, NA))) %>%
    filter(!is.na(Colonisation)) %>%
    mutate_at(.vars = c(7, 8, 9, 10, 13:18, 20), .funs = function(x) stdize(x, prefix = F))


  m.col <- glmer(Colonisation ~ STI_rel * PC3 +
    STI_rel * PC4 +
    STI_rel * PC1 +
    PC1 * PLAND * CLUMPY +
    STI_rel * PLAND * CLUMPY +
    Habitat + X * Y +
    (1 | gridCell50 / Site) + (1 | Species),
  family = binomial,
  data = subset(dat.Col, dat.Col$Scale == 20000),
  nAGQ = 0, control = glmerControl(
    optimizer = "nloptwrap",
    optCtrl = list(maxfun = 1e10),
    calc.derivs = FALSE
  ), verbose = F
  )

  # summary(m.col)


  # extinction
  dat.Ext <- data %>%
    mutate(Extinction = ifelse(event == "loss", 1, ifelse(event == "persistence", 0, NA))) %>%
    filter(!is.na(Extinction)) %>%
    mutate_at(.vars = c(7, 8, 9, 10, 13:18, 20), .funs = function(x) stdize(x, prefix = F))


  m.ext <- glmer(Extinction ~ STI_rel * PC3 +
    STI_rel * PC4 +
    STI_rel * PC1 +
    PC1 * PLAND * CLUMPY +
    STI_rel * PLAND * CLUMPY +
    Habitat + X * Y +
    (1 | gridCell50 / Site) + (1 | Species),
  family = binomial,
  data = subset(dat.Ext, dat.Ext$Scale == 20000),
  nAGQ = 0, control = glmerControl(
    optimizer = "nloptwrap",
    optCtrl = list(maxfun = 1e10),
    calc.derivs = FALSE
  ), verbose = F
  )

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

write_csv(res.mods, "../rev_res_mods_221119.csv")
