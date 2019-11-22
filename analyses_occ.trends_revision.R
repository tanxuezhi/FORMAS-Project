library(unmarked)
library(tidyverse)
library(data.table)

# habitat data
frag_data <- read_csv("../Connectivity - Fragmentation/Fragmentation/Frag_indices_Allhab.csv")
sites_FIN <- read_csv(file = "../Data/Butterflies - Finland/Sites_FIN_ETRS89_landcover.csv")

# butterfly data
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
  group_by(Site, Species, Year, Date) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  complete(nesting(Site, Year, Date), Species) %>%
  ungroup() %>%
  mutate(
    n = ifelse(is.na(n), 0, 1),
    Year = year(as.Date(Date, format = "%d.%m.%Y")),
    Date = yday(as.Date(Date, format = "%d.%m.%Y"))
  ) %>%
  ungroup() %>%
  mutate(Site = paste0(Site, "_FIN"))

countsFIN.fil <- countsFIN %>%
  group_by(Site) %>%
  summarise(n = n_distinct(Year)) %>%
  filter(n > 8)
countsFIN <- countsFIN %>% filter(Site %in% countsFIN.fil$Site)

res <- c()
for(j in unique(countsFIN$Species)){
  for (i in unique(countsFIN$Site)) {
    dat <- countsFIN %>%
      filter(Site == i, Species == j)
    dat$period <- cut_interval(dat$Year, 3, labels = F)
    dat <- dat %>%
      group_by(period) %>%
      summarise(n = max(n))
    
    if(dat[1,2] == 0 & dat[3,2] == 1){
      event = "gain"
    }
    if(dat[1,2] == 0 & dat[3,2] == 0){
      event = "absence"
    }
    if(dat[1,2] == 1 & dat[3,2] == 1){
      event = "persistence"
    }
    if(dat[1,2] == 1 & dat[3,2] == 0){
      event = "loss"
    }
    
    res <- rbind.data.frame(res, cbind.data.frame(Species = j, Site = i, event))
  }
}



res <- c()
for (i in unique(countsFIN$Species)) {
  dat <- countsFIN %>%
    filter(Species == i) %>%
    select(-Species) %>%
    ungroup() %>%
    left_join(sites_FIN %>% mutate(Site = paste0(Site, "_FIN"))) %>%
    left_join(frag_data) %>%
    mutate(yr = as.factor(Year), landcover = as.factor(Landcover)) %>%
    filter(Scale == 10000, Habitat == "Generalist")
  
  dat.f <- formatMult(dat[, c(2, 1, 3, 4, 2)])
  siteCovs(dat.f) <- dat %>%
    select(1, 5, 6, 7, 10, 11) %>%
    group_by(Site) %>%
    summarise_all(unique) %>%
    mutate(Landcover = as.factor(Landcover))
  yearlySiteCovs(dat.f) <- list(yr = as.factor(dat$Year))
  # summary(dat.f)
  
  dat <- dat %>% filter(Year == 2010)
  dat.o <- unmarkedFrameOccu(dat %>% select(1, 3, 4) %>% spread("Date", "n") %>% select(-1),
                             siteCovs = dat %>% select(1, 5, 6, 7, 10, 11) %>% group_by(Site) %>% summarise_all(unique) %>%
                               mutate(Landcover = as.factor(Landcover)),
                             obsCovs = list(Date = matrix(rep(unique(dat$Date), length(unique(dat$Site))),
                                                          nrow = length(unique(dat$Site)), byrow = F
                             ))
  )
  summary(dat.o)
  
  fm.o <- occu(~ Date + I(Date^2) ~ X + Y + Landcover + PLAND, data = dat.o, se = F)
  
  fm <- colext(
    psiformula = ~ PLAND + X + Y,
    gammaformula = ~PLAND,
    epsilonformula = ~PLAND,
    pformula = ~ as.factor(Year) + Date + I(Date^2),
    data = dat.f,
    se = F,
    method = "Nelder-Mead", control = list(maxit = 1e9)
  )
  
  
  pred <- cbind.data.frame(
    Site = unique(dat$Site),
    Species = i,
    init = predict(fm,
                   type = "psi",
                   newdata = data.frame(Site = unique(dat$Site))
    )[, 1],
    Col.prob = predict(fm,
                       type = "col",
                       newdata = data.frame(Site = unique(dat$Site))
    )[, 1],
    Ext.prob = predict(fm,
                       type = "ext",
                       newdata = data.frame(Site = unique(dat$Site))
    )[, 1]
  )
  
  res <- rbind.data.frame(res, pred)
}
