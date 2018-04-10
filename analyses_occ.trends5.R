library(data.table)
library(visreg)
library(MuMIn)
library(lmerTest)
source("functions.R")

butterflies <- as.tbl(fread("C:/Local Folder (c)/butterflies_occ.csv"))

dat.occ.trend <- sp_site_occupancy_trend(butterflies)
write_csv(dat.occ.trend,"../sp_site_coefs.csv")

dat <- read_csv("../sp_site_coefs.csv")

dat <- left_join(dat, butterflies %>% dplyr::select(2,3,5,6,10:20,23) %>% group_by(Species, Site, Scale) %>% summarise_all(first))

dat_col <- dat %>% mutate(colProba = ifelse(Estimate <= 0, 0, 1)) %>% 
  stdize(., prefix = F, omit.cols = c("colProba", "Scale", "nYear"))
dat_ext <- dat %>% mutate(extProba = ifelse(Estimate < 0, 1, 0)) %>% 
  stdize(., prefix = F, omit.cols = c("extProba", "Scale", "nYear"))


m.col <- glmer(colProba ~  STI_rel + PC3 + PC4 + 
                 PLAND * CLUMPY * PC1 + 
                 PLAND * CLUMPY * I(PC1^2) + 
                 Habitat + X*Y + 
                 (1|Site) + (1|Species), 
               family = binomial,
               data = dat_col %>% 
                 dplyr::filter(Scale == 20000, nYear > 10),
               nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                             optCtrl=list(maxfun=1e10),
                                             calc.derivs = FALSE), verbose = F)

summary(m.col)
bandTF = T

par(mfrow=c(2,2))
visreg(m.col, xvar = "PC1", scale = "response", rug = F, by = "CLUMPY", breaks = c(-1,1), 
       overlay = T, band=bandTF, cond = list(PLAND = -1), ylim = c(0,1))
visreg(m.col, xvar = "PC1", scale = "response", rug = F, by = "CLUMPY", breaks = c(-1,1),
       overlay = T, band=bandTF, cond = list(PLAND = 1), ylim = c(0,1))

visreg(m.col, xvar = "STI_rel", scale = "response", rug = F, band=bandTF, ylim = c(0,1))
visreg(m.col, xvar = "PC3", scale = "response", rug = F, band=bandTF, ylim = c(0,1))


m.ext <- glmer(extProba ~  STI_rel + PC3 + PC4 + 
                 PLAND * CLUMPY * PC1 + 
                 PLAND * CLUMPY * I(PC1^2) + 
                 Habitat + X*Y + 
                 (1|Site) + (1|Species), 
               family = binomial,
               data = dat_ext %>% 
                 dplyr::filter(Scale == 20000, nYear > 10),
               nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                             optCtrl=list(maxfun=1e10),
                                             calc.derivs = FALSE), verbose = F)

summary(m.ext)
bandTF = T

par(mfrow=c(2,2))
visreg(m.ext, xvar = "PC1", scale = "response", rug = F, by = "CLUMPY", breaks = c(-1,1), 
       overlay = T, band=bandTF, cond = list(PLAND = -1), ylim = c(0,1))
visreg(m.ext, xvar = "PC1", scale = "response", rug = F, by = "CLUMPY", breaks = c(-1,1), 
       overlay = T, band=bandTF, cond = list(PLAND = 1), ylim = c(0,1))

visreg(m.ext, xvar = "STI_rel", scale = "response", rug = F, band=bandTF, ylim = c(0,1))
visreg(m.ext, xvar = "PC3", scale = "response", rug = F, band=bandTF, ylim = c(0,1))
