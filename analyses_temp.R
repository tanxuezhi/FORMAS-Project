library(dplyr)
library(ggplot2)
library(cowplot)
library(MuMIn)
library(data.table)
library(gamm4)
library(lsmeans)
library(MuMIn)
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
m_frag_butterflies_Ab1_NL <- lme(cti ~ CLUMPY * PLAND * Year + LABEL3, 
                               random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                               data = data.frame(std.butterflies.data.scale_Ab_NL), na.action = na.exclude())
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

########################
#### Post-hoc tests ####
########################
library(lsmeans)


###############
# butterflies #
###############

m_frag_butterflies_Ab1 <- lme(cti ~ PLAND * CLUMPY * Year + LABEL3 + X*Y, 
                                random =  ~1|country/Site, correlation = corCAR1(form = ~Year|country/Site),
                                data = data.frame(na.omit(std.butterflies.data.scale_Ab[,c(2:7,10:11,16,18)])))
anova(m_frag_butterflies_Ab1)

lsm_butterflies_Ab1 <- lsmeans(m_frag_butterflies_Ab1, specs = ~ CLUMPY|PLAND, 
               at = list(PLAND = c(-1,0,1), CLUMPY = c(-1,0,1)), 
               trend = "Year")
PH_butterflies_Ab1 <- as.data.frame(summary(lsm_butterflies_Ab1))
pairs(lsm_butterflies_Ab1)

m_frag_butterflies_P1 <- lme(cti ~ PLAND * CLUMPY * Year + LABEL3 + X*Y, 
                              random =  ~1|country/Site, correlation = corCAR1(form = ~Year|country/Site),
                              data = data.frame(na.omit(std.butterflies.data.scale_P[,c(2:7,10:11,16,18)])))
anova(m_frag_butterflies_P1)

lsm_butterflies_P1 <- lsmeans(m_frag_butterflies_P1, specs = ~ CLUMPY|PLAND, 
                               at = list(PLAND = c(-1,0,1), CLUMPY = c(-1,0,1)), 
                               trend = "Year")
PH_butterflies_P1 <- as.data.frame(summary(lsm_butterflies_P1))
pairs(lsm_butterflies_P1)

PH_butterflies <- rbind.data.frame(cbind.data.frame(type ="Abundance", PH_butterflies_Ab1),
                                   cbind.data.frame(type ="Presence", PH_butterflies_P1))

colnames(PH_butterflies)[2:4] <- c("Fragmentation", "Area_of_SNH", "CTI_trend")
PH_butterflies[,2] <- as.factor(PH_butterflies[,2])
levels(PH_butterflies[,2]) <- c("High","Medium","Low")
PH_butterflies[,3] <- as.factor(PH_butterflies[,3])
levels(PH_butterflies[,3]) <- c("Little area of SNH","Medium area of SMH", "Large area of SNH")

ggplot(PH_butterflies, aes(x = Fragmentation, y = CTI_trend)) + 
  facet_grid(type ~ Area_of_SNH) + 
  geom_bar(stat="identity") +
  scale_y_continuous("CTI trend") +
  geom_errorbar(aes(ymin = CTI_trend - SE, ymax = CTI_trend + SE), width = .2)


#########
# birds #
#########


m_frag_birds_Ab1 <- lme(cti ~ PLAND * CLUMPY * Year + LABEL3 + X*Y, 
                              random =  ~1|Site, correlation = corCAR1(form = ~Year|Site),
                              data = data.frame(na.omit(std.birds.data.scale_Ab[,c(2:7,10:11,16)])))
anova(m_frag_birds_Ab1)

lsm_birds_Ab1 <- lsmeans(m_frag_birds_Ab1, specs = ~ CLUMPY|PLAND, 
                               at = list(PLAND = c(-1,0,1), CLUMPY = c(-1,0,1)), 
                               trend = "Year")
PH_birds_Ab1 <- as.data.frame(summary(lsm_birds_Ab1))
pairs(lsm_birds_Ab1)

m_frag_birds_P1 <- lme(cti ~ PLAND * CLUMPY * Year + LABEL3 + X*Y, 
                             random =  ~1|Site, correlation = corCAR1(form = ~Year|Site),
                             data = data.frame(na.omit(std.birds.data.scale_P[,c(2:7,10:11,16)])))
anova(m_frag_birds_P1)

lsm_birds_P1 <- lsmeans(m_frag_birds_P1, specs = ~ CLUMPY|PLAND, 
                              at = list(PLAND = c(-1,0,1), CLUMPY = c(-1,0,1)), 
                              trend = "Year")
PH_birds_P1 <- as.data.frame(summary(lsm_birds_P1))
pairs(lsm_birds_P1)

PH_birds <- rbind.data.frame(cbind.data.frame(type ="Abundance", PH_birds_Ab1),
                                   cbind.data.frame(type ="Presence", PH_birds_P1))

colnames(PH_birds)[2:4] <- c("Fragmentation", "Area_of_SNH", "CTI_trend")
PH_birds[,2] <- as.factor(PH_birds[,2])
levels(PH_birds[,2]) <- c("High","Medium","Low")
PH_birds[,3] <- as.factor(PH_birds[,3])
levels(PH_birds[,3]) <- c("Little area of SNH","Medium area of SMH", "Large area of SNH")

ggplot(PH_birds, aes(x = Fragmentation, y = CTI_trend)) + 
  facet_grid(type ~ Area_of_SNH) + 
  geom_bar(stat="identity") +
  scale_y_continuous("CTI trend") +
  geom_errorbar(aes(ymin = CTI_trend - SE, ymax = CTI_trend + SE), width = .2)


##################
boxplot()
aov(PLAND ~ LABEL3, data = data.frame(std.butterflies.data.scale_Ab %>% filter))
