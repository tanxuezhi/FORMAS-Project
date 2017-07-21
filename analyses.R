library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(cowplot)
library(MuMIn)
library(HH)
library(data.table)
library(lsmeans)
library(mgcv)
source("functions.R")


###############################
########## Load data ########## 
###############################

butterflies.data <- as.tbl(fread("../Data/cti_butterflies_data.csv"))
birds.data <- as.tbl(fread("../Data/cti_birds_data.csv"))


##############################
#### cti change over time ####
##############################

birds.data.sel <- birds.data %>% filter(Scale == 1000, type == "Abundance")
butterflies.data.sel <- butterflies.data %>% filter(Scale == 1000, type == "Abundance")

par(mfrow=c(1,2))

birds.data.sel$Year <- as.factor(birds.data.sel$Year)
m1 <- lmer(cti ~ Year + (1|Site), data = birds.data.sel)
ctiYear_lsmeans1 <- lsmeans(m1, specs = "Year")
ctiYear_lsmeans1 <- summary(ctiYear_lsmeans1)
ctiYear_lsmeans1$Year <- as.vector(ctiYear_lsmeans1$Year)
plot(lsmean ~ Year, ctiYear_lsmeans1, type = "o", pch = 16, main = "Birds - Sweden", ylab = "CTI")

butterflies.data.sel$Year <- as.factor(butterflies.data.sel$Year)
m2 <- lmer(cti ~ Year + (1|Site), data = subset(butterflies.data.sel, butterflies.data.sel$country == "NL"))
ctiYear_lsmeans2 <- lsmeans(m2, specs = "Year")
ctiYear_lsmeans2 <- summary(ctiYear_lsmeans2)
ctiYear_lsmeans2$Year <- as.vector(ctiYear_lsmeans2$Year)
plot(lsmean ~ Year, ctiYear_lsmeans2, type = "o", pch = 16, main = "Butterflies - NL & FIN", ylab = "CTI", ylim = c(7.9,9.25), col = "blue")

m3 <- lmer(cti ~ Year + (1|Site), data = subset(butterflies.data.sel, butterflies.data.sel$country == "FIN"))
ctiYear_lsmeans3 <- lsmeans(m3, specs = "Year")
ctiYear_lsmeans3 <- summary(ctiYear_lsmeans3)
ctiYear_lsmeans3$Year <- as.vector(ctiYear_lsmeans3$Year)
points(lsmean ~ Year, ctiYear_lsmeans3, type = "o", pch = 16, col = "red")


###################################
##### Effect of fragmentation #####
###################################

##### Test scale effect ####
## select only sites without NA
butterflies.data.sel <- butterflies.data %>% 
  filter(!Site %in% c((butterflies.data %>% filter(is.na(CLUMPY)) %>% count(Site))[,"Site"])$Site)

birds.data.sel <- birds.data %>% 
  filter(!Site %in% c((birds.data %>% filter(is.na(CLUMPY)) %>% count(Site))[,"Site"])$Site)

## run models across scales ## (export plots 800 x 600)
par(mfrow=c(2,2))
scaleFrag_butterflies_Ab <- scaleTest(butterflies.data %>% filter(type == "Abundance"), 
                                      main = "Butterflies - Abundance")
scaleFrag_birds_Ab <- scaleTest(birds.data %>% filter(type == "Abundance"), 
                                main = "Birds - Abundance")
scaleFrag_butterflies_P <- scaleTest(butterflies.data %>% filter(type == "Presence"), 
                                     main = "Butterflies - Presence")
scaleFrag_birds_P <- scaleTest(birds.data %>% filter(type == "Presence"), 
                               main = "Birds - Presence")

# merge results
scaleFrag <- rbind(cbind.data.frame(species = "Butterflies", type = "Abundance", scaleFrag_butterflies_Ab),
                   cbind.data.frame(species = "Birds", type = "Abundance", scaleFrag_birds_Ab),
                   cbind.data.frame(species = "Butterflies", type = "Presence", scaleFrag_butterflies_P),
                   cbind.data.frame(species = "Birds", type = "Presence", scaleFrag_birds_P))

# save results
save.image("scaleTests.RData")
# load results
load("scaleTests.RData")

# plot results
ggplot(data = scaleFrag[scaleFrag$Variable %in% c("PLAND:Year", "CLUMPY:Year", "CLUMPY:PLAND:Year"),],
       aes(x= Scale, y = Estimate, color = Variable)) + facet_grid(type ~ species) +
  geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE), width = 0) +
  geom_point(size = 2) + 
  geom_text(aes(x= Scale, y = (Estimate + .04 * diff(range(Estimate))),
                label = Significance))+
  geom_hline(yintercept=0, lty = 2) +
  geom_line() + 
  scale_color_manual("Interactions", 
                     values = c("#EE7600", "#698B69", "#00B2EE"), 
                     labels = c("Clumpiness x\n% SNH x Year", 
                                "Clumpiness x Year",
                                "% SNH x Year")) +
  scale_x_continuous("Spatial scale (m)")

## select data at optimal scale (not reliable) ##
# butterflies.data.scale.Ab <- butterflies.data %>% 
#   filter(Scale == scaleFrag_butterflies_Ab[order(scaleFrag_butterflies_Ab$AICc),"Scale"][1])
# butterflies.data.scale.P <- butterflies.data %>% 
#   filter(Scale == scaleFrag_butterflies_Ab[order(scaleFrag_butterflies_P$AICc),"Scale"][1])
# birds.data.scale.Ab <- birds.data %>% 
#   filter(Scale == scaleFrag_birds_Ab[order(scaleFrag_birds_Ab$AICc),"Scale"][1])
# birds.data.scale.P <- birds.data %>% 
#   filter(Scale == scaleFrag_birds_Ab[order(scaleFrag_birds_P$AICc),"Scale"][1])

## select data at user-defined scale ##
butterflies.data.scale <- butterflies.data %>% filter(Scale == 10000)
birds.data.scale <- birds.data %>% filter(Scale == 25000)


##### Run analyses ####
## standardize data (or not)
std.butterflies.data.scale <- stdize(butterflies.data.scale %>% filter(type == "Abundance"), prefix = F)
std.butterflies.data.scale$Site <- as.factor(std.butterflies.data.scale$Site)

std.birds.data.scale <- stdize(birds.data.scale.Ab %>% filter(type == "Abundance"), prefix = F)
std.birds.data.scale$Site <- as.factor(std.birds.data.scale$Site)

## lmm ##
m_frag_butterflies_Ab <- lmer(cti ~ CLUMPY * PLAND * Year + LABEL3 + X*Y + (1|country/Site), 
                              data = std.butterflies.data.scale)

anova(m_frag_butterflies_Ab)

## gamm ##
# butterflies
m_frag_butterflies_Ab1 <- gamm(cti ~ CLUMPY * PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                              random = list(country = ~1, Site = ~1), correlation = corCAR1(form = ~Year|Site),
                              data = data.frame(std.butterflies.data.scale))
anova(m_frag_butterflies_Ab1$gam)
summary(m_frag_butterflies_Ab1$gam)

# birds
m_frag_birds_Ab1 <- gamm(cti ~ CLUMPY * PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                               random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                               data = data.frame(std.birds.data.scale))
anova(m_frag_birds_Ab1$gam)
summary(m_frag_birds_Ab1$gam)


# compare with model without correlation structure
m_frag_butterflies_Ab2 <- gamm(cti ~ CLUMPY * PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                               random = list(country = ~1, Site = ~1),
                               data = data.frame(std.butterflies.data.scale))
anova(m_frag_butterflies_Ab2$gam)

model.sel(m_frag_butterflies_Ab1$lme, m_frag_butterflies_Ab2$lme)


## plot three-way interaction ##
# 2d plots
vis.2d(m_frag_butterflies_Ab1$gam, "Year", "CLUMPY", "PLAND", n = 10, origin = std.butterflies.data.scale)
vis.2d(m_frag_birds_Ab1$gam, "Year", "CLUMPY", "PLAND", n = 10, origin = std.birds.data.scale)

# plot by groups
vis.1d(m_frag_butterflies_Ab1$gam, "Year", "CLUMPY", "PLAND", origin = std.butterflies.data.scale)
vis.1d(m_frag_birds_Ab1$gam, "Year", "CLUMPY", "PLAND", origin = std.birds.data.scale)
