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
scaleFrag_butterflies_Ab <- scaleTest(butterflies.data.sel %>% filter(type == "Abundance"), 
                                      main = "Butterflies - Abundance", AIC =T)
scaleFrag_birds_Ab <- scaleTest(birds.data.sel %>% filter(type == "Abundance"), 
                                main = "Birds - Abundance", AIC =T)
scaleFrag_butterflies_P <- scaleTest(butterflies.data.sel %>% filter(type == "Presence"), 
                                     main = "Butterflies - Presence", AIC =T)
scaleFrag_birds_P <- scaleTest(birds.data.sel %>% filter(type == "Presence"), 
                               main = "Birds - Presence", AIC =T)

# merge results
scaleFrag <- rbind(cbind.data.frame(species = "Butterflies", type = "Abundance", scaleFrag_butterflies_Ab),
                   cbind.data.frame(species = "Birds", type = "Abundance", scaleFrag_birds_Ab),
                   cbind.data.frame(species = "Butterflies", type = "Presence", scaleFrag_butterflies_P),
                   cbind.data.frame(species = "Birds", type = "Presence", scaleFrag_birds_P))

# save results
save(scaleFrag, file = "scaleTests.RData")
# load results
load("scaleTests.RData")

# plot results
ggplot(data = scaleFrag[scaleFrag$Variable %in% c("PLAND:Year", "CLUMPY:Year", "CLUMPY:PLAND:Year"),],
       aes(x= Scale, y = Estimate, color = Variable)) + facet_grid(type ~ species) +
  geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE), width = 0) +
  geom_point(size = 2) + 
  geom_text(aes(x= Scale + 500, y = (Estimate + .04 * diff(range(Estimate))),
                label = Significance))+
  geom_hline(yintercept=0, lty = 2) +
  geom_line() + 
  scale_color_manual("Interactions", 
                     values = c("#EE7600", "#698B69", "#00B2EE"), 
                     labels = c("Clumpiness x\n% SNH x Year", 
                                "Clumpiness x Year",
                                "% SNH x Year")) +
  scale_x_continuous("Spatial scale (m)")


## select data at optimal scale  ##
butterflies.data.scale.Ab <- butterflies.data %>%
  filter(Scale == scaleFrag_butterflies_Ab[order(scaleFrag_butterflies_Ab$AICc),"Scale"][1])
butterflies.data.scale.P <- butterflies.data %>%
  filter(Scale == scaleFrag_butterflies_Ab[order(scaleFrag_butterflies_P$AICc),"Scale"][1])
birds.data.scale.Ab <- birds.data %>%
  filter(Scale == scaleFrag_birds_Ab[order(scaleFrag_birds_Ab$AICc),"Scale"][1])
birds.data.scale.P <- birds.data %>%
  filter(Scale == scaleFrag_birds_Ab[order(scaleFrag_birds_P$AICc),"Scale"][1])

## select data at user-defined scale ##
butterflies.data.scale <- butterflies.data %>% filter(Scale == 5000)
birds.data.scale <- birds.data %>% filter(Scale == 5000)


##### Run analyses ####
### gamm ###
## butterflies
# Abundance
std.butterflies.data.scale_Ab <- stdize(butterflies.data.scale.Ab %>% filter(type == "Abundance"), prefix = F)
std.butterflies.data.scale_Ab$Site <- as.factor(std.butterflies.data.scale_Ab$Site)

m_frag_butterflies_Ab1 <- gamm(cti ~ CLUMPY * PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                               random = list(country = ~1, Site = ~1), correlation = corCAR1(form = ~Year|Site),
                               data = data.frame(std.butterflies.data.scale_Ab))
anova(m_frag_butterflies_Ab1$gam)
# summary(m_frag_butterflies_Ab1$gam)

# Presence
std.butterflies.data.scale_P <- stdize(butterflies.data.scale.P %>% filter(type == "Presence"), prefix = F)
std.butterflies.data.scale_P$Site <- as.factor(std.butterflies.data.scale_P$Site)

m_frag_butterflies_P1 <- gamm(cti ~ CLUMPY * PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                              random = list(country = ~1, Site = ~1), correlation = corCAR1(form = ~Year|Site),
                              data = data.frame(std.butterflies.data.scale_P))
anova(m_frag_butterflies_P1$gam)
# summary(m_frag_butterflies_P1$gam)


## birds
# Abundance
std.birds.data.scale_Ab <- stdize(birds.data.scale.Ab %>% filter(type == "Abundance"), prefix = F)
std.birds.data.scale_Ab$Site <- as.factor(std.birds.data.scale_Ab$Site)

m_frag_birds_Ab1 <- gamm(cti ~ CLUMPY * PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                         random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                         data = data.frame(std.birds.data.scale_Ab))
anova(m_frag_birds_Ab1$gam)

# Presence
std.birds.data.scale_P <- stdize(birds.data.scale %>% filter(type == "Presence"), prefix = F)
std.birds.data.scale_P$Site <- as.factor(std.birds.data.scale_P$Site)

m_frag_birds_P1 <- gamm(cti ~ CLUMPY * PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                        random = list(Site = ~1), correlation = corCAR1(form = ~Year|Site),
                        data = data.frame(std.birds.data.scale_P))
anova(m_frag_birds_P1$gam)

# plot by groups (export 800*500)
plot_butterflies_Ab <- vis.1d(m_frag_butterflies_Ab1$gam, "Year", "CLUMPY", "PLAND", origin = std.butterflies.data.scale_Ab)
plot_butterflies_P <- vis.1d(m_frag_butterflies_P1$gam, "Year", "CLUMPY", "PLAND", origin = std.butterflies.data.scale_P)

plot_birds_Ab <- vis.1d(m_frag_birds_Ab1$gam, "Year", "CLUMPY", "PLAND", origin = std.birds.data.scale_Ab)
plot_birds_P <- vis.1d(m_frag_birds_P1$gam, "Year", "CLUMPY", "PLAND", origin = std.birds.data.scale_P)

plot_data <- rbind.data.frame(cbind.data.frame(type= "Abundance", species = "Butterflies", plot_butterflies_Ab),
                              cbind.data.frame(type= "Presence", species = "Butterflies", plot_butterflies_P),
                              cbind.data.frame(type= "Abundance", species = "Birds", plot_birds_Ab),
                              cbind.data.frame(type= "Presence", species = "Birds", plot_birds_P))
write.csv(plot_data, "../plot_data_selected_scale.csv", row.names = F)

plot_data2 <- c()
for(i in unique(plot_data$species)){
  for(j in unique(plot_data$type)){
    temp <- plot_data[plot_data$species == i & plot_data$type == j,]
    temp$PLAND <- as.factor(temp$PLAND)
    temp$PLAND <- reorder(temp$PLAND, order(temp$PLAND, decreasing = F))
    levels(temp$PLAND) <- c("Little area of SNH","Large area of SNH")
    temp$CLUMPY <- as.factor(temp$CLUMPY)
    temp$CLUMPY <- reorder(temp$CLUMPY, order(temp$CLUMPY, decreasing = T))
    levels(temp$CLUMPY) <- c("Highly fragmentated","Moderately fragmentated","Little fragmentated")
    plot_data2 <- rbind.data.frame(plot_data2, temp)
  }
}
plot_data <- plot_data2
rm("temp", "plot_data2")

ggplot(data = subset(plot_data, plot_data$species == "Butterflies"), aes(x = Year, y = pred, color = CLUMPY)) + geom_line() + 
  facet_grid(type ~ PLAND, scales = "free") + 
  scale_y_continuous("Community Temperature Index") +
  scale_x_continuous("Year") +
  labs(color = "Habitat fragmentation") +
  scale_color_manual(values = c("darkorange4", "darkorange3", "orange"))

ggplot(data = subset(plot_data, plot_data$species == "Birds"), aes(x = Year, y = pred, color = CLUMPY)) + geom_line() + 
  facet_grid(type ~ PLAND, scales = "free") + 
  scale_y_continuous("Community Temperature Index") +
  scale_x_continuous("Year") +
  labs(color = "Habitat fragmentation") +
  scale_color_manual(values = c("darkorange4", "darkorange3", "orange"))


# compare with model without correlation structure
m_frag_butterflies_Ab2 <- gamm(cti ~ CLUMPY * PLAND * Year + LABEL3 + s(X,Y, bs = "tp"), 
                               random = list(country = ~1, Site = ~1),
                               data = data.frame(std.butterflies.data.scale))
anova(m_frag_butterflies_Ab2$gam)

model.sel(m_frag_butterflies_Ab1$lme, m_frag_butterflies_Ab$lme)


## plot three-way interaction ##
# 2d plots
vis.2d(m_frag_butterflies_Ab1$gam, "Year", "CLUMPY", "PLAND", n = 10, origin = std.butterflies.data.scale)
vis.2d(m_frag_birds_Ab1$gam, "Year", "CLUMPY", "PLAND", n = 10, origin = std.birds.data.scale)

