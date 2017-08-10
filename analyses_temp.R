library(dplyr)
library(ggplot2)
library(cowplot)
library(MuMIn)
library(data.table)
library(lsmeans)
library(mgcv)
source("functions.R")

###################################
###### load and prepare data ######
###################################

### load data ###
butterflies.data <- as.tbl(fread("../Data/cti_butterflies_data.csv"))
birds.data <- as.tbl(fread("../Data/cti_birds_data.csv"))

### select scale ###
butterflies.data.scale <- butterflies.data %>% filter(Scale == 5000)
birds.data.scale <- birds.data %>% filter(Scale == 5000)

### standardize abundance and presence data ###
std.birds.data.scale_Ab <- stdize(birds.data.scale %>% 
                                    filter(type == "Abundance"), prefix = F)
std.birds.data.scale_Ab$Site <- as.factor(std.birds.data.scale_Ab$Site)

std.birds.data.scale_P <- stdize(birds.data.scale %>% 
                                    filter(type == "Presence"), prefix = F)
std.birds.data.scale_P$Site <- as.factor(std.birds.data.scale_P$Site)


std.butterflies.data.scale_Ab <- stdize(butterflies.data.scale %>% 
                                    filter(type == "Abundance"), prefix = F)
std.butterflies.data.scale_Ab$Site <- as.factor(std.butterflies.data.scale_Ab$Site)

std.butterflies.data.scale_P <- stdize(butterflies.data.scale %>% 
                                   filter(type == "Presence"), prefix = F)
std.butterflies.data.scale_P$Site <- as.factor(std.butterflies.data.scale_P$Site)

# by country
std.butterflies.data.scale_Ab_NL <- stdize(butterflies.data.scale %>% 
                                             filter(type == "Abundance" & country == "NL"), prefix = F)
std.butterflies.data.scale_Ab_NL$Site <- as.factor(std.butterflies.data.scale_Ab_NL$Site)

std.butterflies.data.scale_Ab_FIN <- stdize(butterflies.data.scale %>% 
                                              filter(type == "Abundance" & country == "FIN"), prefix = F)
std.butterflies.data.scale_Ab_FIN$Site <- as.factor(std.butterflies.data.scale_Ab_FIN$Site)

std.butterflies.data.scale_P_NL <- stdize(butterflies.data.scale %>% 
                                            filter(type == "Presence" & country == "NL"), prefix = F)
std.butterflies.data.scale_P_NL$Site <- as.factor(std.butterflies.data.scale_P_NL$Site)

std.butterflies.data.scale_P_FIN <- stdize(butterflies.data.scale %>% 
                                             filter(type == "Presence" & country == "FIN"), prefix = F)
std.butterflies.data.scale_P_FIN$Site <- as.factor(std.butterflies.data.scale_P_FIN$Site)


####################
#### by country ####
####################


## butterflies
## Abundance
# NL
m_frag_butterflies_Ab1_NL <- gamm(cti ~ CLUMPY * PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                               random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                               data = data.frame(std.butterflies.data.scale_Ab_NL))
anova(m_frag_butterflies_Ab1_NL$gam)

# FIN
m_frag_butterflies_Ab1_FIN <- gamm(cti ~ CLUMPY * PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                                  random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                  data = data.frame(std.butterflies.data.scale_Ab_FIN))
anova(m_frag_butterflies_Ab1_FIN$gam)

## Presence
# NL
m_frag_butterflies_P1_NL <- gamm(cti ~ CLUMPY * PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                                  random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                  data = data.frame(std.butterflies.data.scale_P_NL))
anova(m_frag_butterflies_P1_NL$gam)

# FIN
m_frag_butterflies_P1_FIN <- gamm(cti ~ CLUMPY * PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                                   random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                   data = data.frame(std.butterflies.data.scale_P_FIN), 
                                  control=lmeControl(opt='optim'))
anova(m_frag_butterflies_P1_FIN$gam)

### test fit ###
## butterflies ##
## both
# abundance
m_frag_butterflies_Ab1 <- uGamm(cti ~ PLAND * CLUMPY * Year + LABEL3 + s(X,Y, bs = "tp"), 
                          random = list(country = ~1, Site = ~1), correlation = corCAR1(form = ~Year|Site),
                          data = data.frame(na.omit(std.butterflies.data.scale_Ab[,c(2:7,10:11,16,18)])), method = "ML")
m_frag_butterflies_Ab2 <- uGamm(cti ~ CLUMPY * Year + LABEL3 + s(X,Y, bs = "tp"), 
                          random = list(country = ~1, Site = ~1), correlation = corCAR1(form = ~Year|Site),
                          data = data.frame(na.omit(std.butterflies.data.scale_Ab[,c(2:7,10:11,16,18)])), method = "ML")
m_frag_butterflies_Ab3 <- uGamm(cti ~ PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                          random = list(country = ~1, Site = ~1), correlation = corCAR1(form = ~Year|Site),
                          data = data.frame(na.omit(std.butterflies.data.scale_Ab[,c(2:7,10:11,16,18)])), method = "ML")
m_frag_butterflies_Ab4 <- uGamm(cti ~ Year + LABEL3 + s(X,Y, bs = "tp"), 
                          random = list(country = ~1, Site = ~1), correlation = corCAR1(form = ~Year|Site),
                          data = data.frame(na.omit(std.butterflies.data.scale_Ab[,c(2:7,10:11,16,18)])), method = "ML")

model.sel(CLUMPY_PLAND = m_frag_butterflies_Ab1,
          CLUMPY = m_frag_butterflies_Ab2,
          PLAND = m_frag_butterflies_Ab3,
          NoLand = m_frag_butterflies_Ab4)[,12:15]

# presence
m_frag_butterflies_P1 <- uGamm(cti ~ PLAND * CLUMPY * Year + LABEL3 + s(X,Y, bs = "tp"), 
                                random = list(country = ~1, Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                data = data.frame(na.omit(std.butterflies.data.scale_P[,c(2:7,10:11,16,18)])), method = "ML")
m_frag_butterflies_P2 <- uGamm(cti ~ CLUMPY * Year + LABEL3 + s(X,Y, bs = "tp"), 
                                random = list(country = ~1, Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                data = data.frame(na.omit(std.butterflies.data.scale_P[,c(2:7,10:11,16,18)])), method = "ML")
m_frag_butterflies_P3 <- uGamm(cti ~ PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                                random = list(country = ~1, Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                data = data.frame(na.omit(std.butterflies.data.scale_P[,c(2:7,10:11,16,18)])), method = "ML")
m_frag_butterflies_P4 <- uGamm(cti ~ Year + LABEL3 + s(X,Y, bs = "tp"), 
                                random = list(country = ~1, Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                data = data.frame(na.omit(std.butterflies.data.scale_P[,c(2:7,10:11,16,18)])), method = "ML")

model.sel(CLUMPY_PLAND = m_frag_butterflies_P1,
          CLUMPY = m_frag_butterflies_P2,
          PLAND = m_frag_butterflies_P3,
          NoLand = m_frag_butterflies_P4)[,12:15]

## NL
# Abundance
m_frag_butterflies_Ab1 <- uGamm(cti ~ PLAND * CLUMPY * Year + LABEL3 + s(X,Y, bs = "tp"), 
                                random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                data = data.frame(na.omit(subset(std.butterflies.data.scale_Ab, std.butterflies.data.scale_Ab$country == "NL")[,c(2:7,10:11,16,18)])), method = "ML")
m_frag_butterflies_Ab2 <- uGamm(cti ~ CLUMPY * Year + LABEL3 + s(X,Y, bs = "tp"), 
                                random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                data = data.frame(na.omit(subset(std.butterflies.data.scale_Ab, std.butterflies.data.scale_Ab$country == "NL")[,c(2:7,10:11,16,18)])), method = "ML")
m_frag_butterflies_Ab3 <- uGamm(cti ~ PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                                random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                data = data.frame(na.omit(subset(std.butterflies.data.scale_Ab, std.butterflies.data.scale_Ab$country == "NL")[,c(2:7,10:11,16,18)])), method = "ML")
m_frag_butterflies_Ab4 <- uGamm(cti ~ Year + LABEL3 + s(X,Y, bs = "tp"), 
                                random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                data = data.frame(na.omit(subset(std.butterflies.data.scale_Ab, std.butterflies.data.scale_Ab$country == "NL")[,c(2:7,10:11,16,18)])), method = "ML")

model.sel(CLUMPY_PLAND = m_frag_butterflies_Ab1,
          CLUMPY = m_frag_butterflies_Ab2,
          PLAND = m_frag_butterflies_Ab3,
          NoLand = m_frag_butterflies_Ab4)[,12:15]

# presence
m_frag_butterflies_P1 <- uGamm(cti ~ PLAND * CLUMPY * Year + LABEL3 + s(X,Y, bs = "tp"), 
                                random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                data = data.frame(na.omit(subset(std.butterflies.data.scale_P, std.butterflies.data.scale_P$country == "NL")[,c(2:7,10:11,16,18)])), method = "ML")
m_frag_butterflies_P2 <- uGamm(cti ~ CLUMPY * Year + LABEL3 + s(X,Y, bs = "tp"), 
                                random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                data = data.frame(na.omit(subset(std.butterflies.data.scale_P, std.butterflies.data.scale_P$country == "NL")[,c(2:7,10:11,16,18)])), method = "ML")
m_frag_butterflies_P3 <- uGamm(cti ~ PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                                random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                data = data.frame(na.omit(subset(std.butterflies.data.scale_P, std.butterflies.data.scale_P$country == "NL")[,c(2:7,10:11,16,18)])), method = "ML")
m_frag_butterflies_P4 <- uGamm(cti ~ Year + LABEL3 + s(X,Y, bs = "tp"), 
                                random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                data = data.frame(na.omit(subset(std.butterflies.data.scale_P, std.butterflies.data.scale_P$country == "NL")[,c(2:7,10:11,16,18)])), method = "ML")

model.sel(CLUMPY_PLAND = m_frag_butterflies_P1,
          CLUMPY = m_frag_butterflies_P2,
          PLAND = m_frag_butterflies_P3,
          NoLand = m_frag_butterflies_P4)[,12:15]

#FIN
#Abundance
m_frag_butterflies_Ab1 <- uGamm(cti ~ PLAND * CLUMPY * Year + LABEL3 + s(X,Y, bs = "tp"), 
                                random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                data = data.frame(na.omit(subset(std.butterflies.data.scale_Ab, std.butterflies.data.scale_Ab$country == "FIN")[,c(2:7,10:11,16,18)])), method = "ML")
m_frag_butterflies_Ab2 <- uGamm(cti ~ CLUMPY * Year + LABEL3 + s(X,Y, bs = "tp"), 
                                random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                data = data.frame(na.omit(subset(std.butterflies.data.scale_Ab, std.butterflies.data.scale_Ab$country == "FIN")[,c(2:7,10:11,16,18)])), method = "ML")
m_frag_butterflies_Ab3 <- uGamm(cti ~ PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                                random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                data = data.frame(na.omit(subset(std.butterflies.data.scale_Ab, std.butterflies.data.scale_Ab$country == "FIN")[,c(2:7,10:11,16,18)])), method = "ML")
m_frag_butterflies_Ab4 <- uGamm(cti ~ Year + LABEL3 + s(X,Y, bs = "tp"), 
                                random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                data = data.frame(na.omit(subset(std.butterflies.data.scale_Ab, std.butterflies.data.scale_Ab$country == "FIN")[,c(2:7,10:11,16,18)])), method = "ML")

model.sel(CLUMPY_PLAND = m_frag_butterflies_Ab1,
          CLUMPY = m_frag_butterflies_Ab2,
          PLAND = m_frag_butterflies_Ab3,
          NoLand = m_frag_butterflies_Ab4)[,12:15]

# presence
m_frag_butterflies_P1 <- uGamm(cti ~ PLAND * CLUMPY * Year + LABEL3 + s(X,Y, bs = "tp"), 
                               random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                               data = data.frame(na.omit(subset(std.butterflies.data.scale_P, std.butterflies.data.scale_P$country == "FIN")[,c(2:7,10:11,16,18)])), method = "ML", 
                               control=lmeControl(opt='optim'))
m_frag_butterflies_P2 <- uGamm(cti ~ CLUMPY * Year + LABEL3 + s(X,Y, bs = "tp"), 
                               data = data.frame(na.omit(subset(std.butterflies.data.scale_P, std.butterflies.data.scale_P$country == "FIN")[,c(2:7,10:11,16,18)])), method = "ML", 
                               control=lmeControl(opt='optim'))
m_frag_butterflies_P3 <- uGamm(cti ~ PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                               data = data.frame(na.omit(subset(std.butterflies.data.scale_P, std.butterflies.data.scale_P$country == "FIN")[,c(2:7,10:11,16,18)])), method = "ML", 
                               control=lmeControl(opt='optim'))
m_frag_butterflies_P4 <- uGamm(cti ~ Year + LABEL3 + s(X,Y, bs = "tp"), 
                               data = data.frame(na.omit(subset(std.butterflies.data.scale_P, std.butterflies.data.scale_P$country == "FIN")[,c(2:7,10:11,16,18)])), method = "ML", 
                               control=lmeControl(opt='optim'))

model.sel(CLUMPY_PLAND = m_frag_butterflies_P1,
          CLUMPY = m_frag_butterflies_P2,
          PLAND = m_frag_butterflies_P3,
          NoLand = m_frag_butterflies_P4)[,13:16]
# birds
# Abundance
m_frag_birds_Ab1 <- uGamm(cti ~ PLAND * CLUMPY * Year + LABEL3 + s(X,Y, bs = "tp"), 
                                random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                                data = data.frame(na.omit(std.birds.data.scale_Ab[,c(2:7,10:11,16)])), method = "ML")
m_frag_birds_Ab2 <- uGamm(cti ~ CLUMPY * Year + LABEL3 + s(X,Y, bs = "tp"), 
                          random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                          data = data.frame(na.omit(std.birds.data.scale_Ab[,c(2:7,10:11,16)])), method = "ML")
m_frag_birds_Ab3 <- uGamm(cti ~ PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                          random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                          data = data.frame(na.omit(std.birds.data.scale_Ab[,c(2:7,10:11,16)])), method = "ML")
m_frag_birds_Ab4 <- uGamm(cti ~ Year + LABEL3 + s(X,Y, bs = "tp"), 
                          random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                          data = data.frame(na.omit(std.birds.data.scale_Ab[,c(2:7,10:11,16)])), method = "ML")

model.sel(CLUMPY_PLAND = m_frag_birds_Ab1,
          CLUMPY = m_frag_birds_Ab2,
          PLAND = m_frag_birds_Ab3,
          NoLand = m_frag_birds_Ab4)[,12:15]

# Presence
m_frag_birds_P1 <- uGamm(cti ~ PLAND * CLUMPY * Year + LABEL3 + s(X,Y, bs = "tp"), 
                          random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                          data = data.frame(na.omit(std.birds.data.scale_P[,c(2:7,10:11,16)])), method = "ML")
m_frag_birds_P2 <- uGamm(cti ~ CLUMPY * Year + LABEL3 + s(X,Y, bs = "tp"), 
                          random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                          data = data.frame(na.omit(std.birds.data.scale_P[,c(2:7,10:11,16)])), method = "ML")
m_frag_birds_P3 <- uGamm(cti ~ PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                          random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                          data = data.frame(na.omit(std.birds.data.scale_P[,c(2:7,10:11,16)])), method = "ML")
m_frag_birds_P4 <- uGamm(cti ~ Year + LABEL3 + s(X,Y, bs = "tp"), 
                          random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                          data = data.frame(na.omit(std.birds.data.scale_P[,c(2:7,10:11,16)])), method = "ML")

model.sel(CLUMPY_PLAND = m_frag_birds_P1,
          CLUMPY = m_frag_birds_P2,
          PLAND = m_frag_birds_P3,
          NoLand = m_frag_birds_P4)[,12:15]
