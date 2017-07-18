library(rio)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(cowplot)
library(visreg)
library(MuMIn)
library(HH)
library(fields)
source("functions.R")

##### import CTI data and merge with site data #####
lc_class <- rio::import(file = "../Landcover/clc_legend.xls", which = 1L)

### butterflies ###
# sites data
sites_NL <- read.csv(file = "../Data/Butterflies - Netherlands/Sites_NL_ETRS89_landcover.csv")
sites_FIN <- read.csv(file = "../Data/Butterflies - Finland/Sites_FIN_ETRS89_landcover.csv")

# fragmentation data
frag_data_NL <- read.csv("../Connectivity/Fragmentation/NL/Frag_indices.csv")
frag_data_FIN <- read.csv("../Connectivity/Fragmentation/FIN/Frag_indices.csv")

sites_NL <- merge(sites_NL, frag_data_NL, by.x = "Site", by.y = "LID")
sites_NL <- merge(sites_NL, lc_class, by.x = "Landcover", by.y = "CLC_CODE")
sites_NL <- sites_NL[,-which(names(sites_NL) %in% "GRID_CODE")]

sites_FIN <- merge(sites_FIN, frag_data_FIN, by.x = "Site", by.y = "LID")
sites_FIN <- merge(sites_FIN, lc_class, by.x = "Landcover", by.y = "GRID_CODE")
sites_FIN <- sites_FIN[,-which(names(sites_FIN) %in% "CLC_CODE")]

# cti data
# Abundance-weighted #
cti_AB_NL <- read.csv2("../Data/Butterflies - Netherlands/cti_site_year_abundance_2017-05-19.csv", dec = ".")[,-3]
cti_AB_NL <- merge(cti_AB_NL, sites_NL, by.x = "Site", by.y = "Site")

cti_AB_FIN <- read.csv("../Data/Butterflies - Finland/CTI_Abundance_FINLAND_1999-2016.csv", dec = ".")
cti_AB_FIN <- merge(cti_AB_FIN, sites_FIN, by.x = "Site", by.y = "Site")

cti_AB_butterflies <- as.tbl(bind_rows(cbind.data.frame(cti_AB_NL, country = "NL"), cbind.data.frame(cti_AB_FIN, country = "FIN")))
cti_AB_butterflies$Site <- paste0(cti_AB_butterflies$Site, "_", cti_AB_butterflies$country)

# P/A #
cti_PA_NL <- read.csv2("../Data/Butterflies - Netherlands/cti_site_year_presence_2017-05-19.csv", dec = ".")[,-3]
cti_PA_NL <- merge(cti_PA_NL, sites_NL, by.x = "Site", by.y = "Site")

cti_PA_FIN <- read.csv("../Data/Butterflies - Finland/CTI_presence_FINLAND_1999-2016.csv", dec = ".")
cti_PA_FIN <- merge(cti_PA_FIN, sites_FIN, by.x = "Site", by.y = "Site")

cti_PA_butterflies <- as.tbl(bind_rows(cbind.data.frame(cti_PA_NL, country = "NL"), cbind.data.frame(cti_PA_FIN, country = "FIN")))
cti_PA_butterflies$Site <- paste0(cti_PA_butterflies$Site, "_", cti_PA_butterflies$country)


### birds ###
# sites data
sites_SWE <- read.csv(file = "../Data/Birds - Sweden/Sites_SWE_ETRS89_landcover.csv")

# fragmentation data
frag_data_SWE <- read.csv("../Connectivity/Fragmentation/SWE/Frag_indices.csv")

sites_SWE <- merge(sites_SWE, frag_data_SWE, by.x = "Site", by.y = "LID")
sites_SWE <- merge(sites_SWE, lc_class, by.x = "Landcover", by.y = "GRID_CODE")

# cti data
# Abundance-weighted #
cti_AB_SWE <- read.csv("../Data/Birds - Sweden/CTI_abundance_Sweden_1996-2016.csv", dec = ".")
cti_AB_SWE <- as.tbl(inner_join(cti_AB_SWE, sites_SWE, by = "Site"))

# P/A #
cti_PA_SWE <- read.csv("../Data/Birds - Sweden/CTI_presence_Sweden_1996-2016.csv", dec = ".")
cti_PA_SWE <- as.tbl(inner_join(cti_PA_SWE, sites_SWE, by = "Site"))


##### try first preliminary analyses #####

## cti change over time
par(mfrow=c(1,2))

cti_AB_SWE$Year <- as.factor(cti_AB_SWE$Year)
m1 <- lmer(cti ~ Year + (1|Site), data = cti_AB_SWE)
ctiYear_lsmeans <- lsmeansLT(m1, test.effs = "Year")
ctiYear_lsmeans$lsmeans.table$Year <- as.vector(ctiYear_lsmeans$lsmeans.table$Year)
plot(Estimate ~ Year, ctiYear_lsmeans$lsmeans.table, type = "o", pch = 16, main = "Birds - Sweden", ylab = "CTI")

cti_AB_butterflies$Year <- as.factor(cti_AB_butterflies$Year)
m2 <- lmer(cti ~ Year + (1|Site), data = subset(cti_AB_butterflies, cti_AB_butterflies$country == "NL"))
ctiYear_lsmeans <- lsmeansLT(m2, test.effs = c("Year"))
ctiYear_lsmeans$lsmeans.table$Year <- as.vector(ctiYear_lsmeans$lsmeans.table$Year)
plot(Estimate ~ Year, ctiYear_lsmeans$lsmeans.table, type = "o", pch = 16, main = "Butterflies - NL & FIN", ylab = "CTI", ylim = c(7.9,9.25), col = "blue")

m3 <- lmer(cti ~ Year + (1|Site), data = subset(cti_AB_butterflies, cti_AB_butterflies$country == "FIN"))
ctiYear_lsmeans <- lsmeansLT(m3, test.effs = c("Year"))
ctiYear_lsmeans$lsmeans.table$Year <- as.vector(ctiYear_lsmeans$lsmeans.table$Year)
points(Estimate ~ Year, ctiYear_lsmeans$lsmeans.table, type = "o", pch = 16, col = "red")


## effect of NDVI variability

# test scale
cti_AB <- subset(cti_AB, !cti_AB$Site %in% unique(subset(cti_AB, is.na(cti_AB$stSD.ndvi))$Site))

m_ltv_AB <- list()
for(i in unique(cti_AB$agreg.fac)){
  m_ltv_AB[[i]] <- lmer(cti ~ Year * sd.ndvi + LABEL1 + (1|km10square/Site), data = filter(cti_AB, agreg.fac == i), na.action = na.fail)
}
m_ltv_AB <- unlist(m_ltv_AB)
names(m_ltv_AB) <- paste("aggregation factor =", unique(cti_AB$agreg.fac))
sel.scale <- model.sel(m_ltv_AB)


# long-term variability
m_ltv_AB <- lmer(cti ~ Year * ltSD.ndvi + LABEL2 + (1|km10square/Site), data = cti_AB)
summary(m_ltv_AB)
pred.m_ltv_AB <- visreg(m_ltv_AB, xvar = "Year", ylab = "CTI", scale = "response", 
                        by = "ltSD.ndvi", gg = T, breaks = 3, plot = F)
p.ltv_AB <- ggplot(data = pred.m_ltv_AB$fit, aes(x = Year, y = visregFit, color = as.factor(ltSD.ndvi))) + 
  geom_line() + 
  scale_y_continuous(name = "CTI") +
  theme(legend.position = c(0.8,0.1)) + 
  scale_color_manual(name = "Long-term NDVI variability", labels = c("Low", "Medium", "High"), values = c("royalblue1", "palegreen3", "tomato2")) 


# short-term variability
m_stv_AB <- lmer(cti ~ Year * stSD.ndvi + LABEL2 + (1|km10square/Site), data = cti_AB)
summary(m_stv_AB)
pred.m_ltv_AB <- visreg(m_stv_AB, xvar = "Year", ylab = "CTI", scale = "response", 
                        by = "stSD.ndvi", gg = T, breaks = 3, plot = F)
p.stv_AB <- ggplot(data = pred.m_ltv_AB$fit, aes(x = Year, y = visregFit, color = as.factor(stSD.ndvi))) + 
  geom_line() + 
  scale_y_continuous(name = "CTI") +
  theme(legend.position = c(0.8,0.1)) + 
  scale_color_manual(name = "Short-term NDVI variability", labels = c("Low", "Medium", "High"), values = c("royalblue1", "palegreen3", "tomato2")) 

# plot
plot_grid(p.ltv_AB, p.stv_AB, ncol = 2)


##### Effect of fragmentation #####
# response variable  = cti
vif(data.frame(cti_AB_butterflies[,c("LSI","PLAND","Year")]))
vif(data.frame(cti_AB_SWE[,c("LSI","PLAND","Year")]))

cti_AB_butterflies$LABEL3 <- as.factor(cti_AB_butterflies$LABEL3)
cti_AB_SWE$LABEL3 <- as.factor(cti_AB_SWE$LABEL3)

cti_AB_butterflies %>% group_by(Scale) %>% summarise(n = length(Site))
cti_AB_SWE %>% group_by(Scale) %>% summarise(n = length(Site))

cti_AB_butterflies <- cti_AB_butterflies %>% filter(Scale > 2000)
cti_AB_SWE <- cti_AB_SWE %>% filter(Scale >= 3000)

par(mfrow=c(1,2))
scaleFrag_butterflies <- scaleTest(cti_AB_butterflies)
scaleFrag_birds <- scaleTest(cti_AB_SWE)

cti_AB_butterflies_sel <- cti_AB_butterflies[cti_AB_butterflies$Scale == scaleFrag_butterflies[scaleFrag_butterflies$AICc == min(scaleFrag_butterflies$AICc),1],]
cti_AB_SWE_sel <- cti_AB_SWE[cti_AB_SWE$Scale == scaleFrag_birds[scaleFrag_birds$AICc == min(scaleFrag_birds$AICc),1],]

m_frag_butterflies <- lmer(cti ~ CLUMPY * PLAND * Year + X*Y + (1|country/Site), data = cti_AB_butterflies_sel, na.action = na.fail)
m_frag_butterflies_std <- lmer(cti ~ CLUMPY * PLAND * Year + X*Y + (1|country/Site), data = stdize(cti_AB_butterflies_sel, prefix = F), na.action = na.fail)
anova(m_frag_butterflies_std)
summary(m_frag_butterflies_std)
sjp.lmer(m_frag_butterflies_std, type = "fe")

m_frag_SWE <- lmer(cti ~ LSI * PLAND * Year + LABEL3 + X*Y + (1|Site), data = cti_AB_SWE_sel, na.action = na.fail)
m_frag_SWE_std <- lmer(cti ~ LSI * PLAND * Year + X*Y + LABEL3 + (1|Site), data = stdize(cti_AB_butterflies_sel, prefix = F), na.action = na.fail)
anova(m_frag_SWE_std)
summary(m_frag_SWE_std)


ncf:::spline.correlog(cti_AB$X, cti_AB$Y, resid(m_frag_std, type  = "pearson"),
                      resamp=100, na.rm = T)

# plot three-way interaction
# 2d plots

a <- vis.2d(m_frag_butterflies_std, "Year", "CLUMPY", "PLAND", n = 50)

# grouped by classes (both factors)
PLAND_val <- c(mean(cti_AB_butterflies_sel$PLAND)-sd(cti_AB_butterflies_sel$PLAND), mean(cti_AB_butterflies_sel$PLAND), mean(cti_AB_butterflies_sel$PLAND)+sd(cti_AB_butterflies_sel$PLAND))
LSI_val <- c(mean(cti_AB_butterflies_sel$CLUMPY)-sd(cti_AB_butterflies_sel$CLUMPY), mean(cti_AB_butterflies_sel$CLUMPY), mean(cti_AB_butterflies_sel$CLUMPY)+sd(cti_AB_butterflies_sel$CLUMPY))

res_class <- c()
for(i in PLAND_val){
  for (j in LSI_val){
    p <- visreg(m_frag_butterflies, xvar = "Year", cond = list(PLAND = i, CLUMPY = j), plot = F)
    res_class <- bind_rows(res_class, p$fit)
  }
}

# version 1
PLAND_names <- c("Small","Medium","Large"); names(PLAND_names) <- PLAND_val
LSI_names <- c("Little fragmentated","Moderately fragmentated","Highly fragmentated"); names(LSI_names) <- LSI_val
names_fac <- c(PLAND_names, LSI_names)

ggplot(data = res_class, aes(x = Year, y = visregFit)) + geom_line() + 
  facet_grid(PLAND ~ CLUMPY, labeller = as_labeller(names_fac)) + 
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.2) +
  scale_y_continuous("Community temperature Index")


# version 2
col.PLAND <- c("dodgerblue", "gold2", "firebrick"); names(col.PLAND) <- PLAND_val

ggplot(data = res_class, aes(x = Year, y = visregFit, color = as.factor(PLAND))) + geom_line() + 
  facet_grid( ~ CLUMPY, labeller = as_labeller(LSI_names)) + 
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr, fill = as.factor(PLAND)), alpha=0.1, colour=NA) +
  scale_y_continuous("Community temperature Index") +
  scale_color_manual("Area of SNH", labels = names_fac, values=col.PLAND) +
  scale_fill_manual("Area of SNH", labels = names_fac, values=col.PLAND)


# grouped by classes (not fragmentation)
PLAND_val <- c(mean(cti_AB$PLAND)-sd(cti_AB$PLAND), mean(cti_AB$PLAND), mean(cti_AB$PLAND)+sd(cti_AB$PLAND))
LSI_val <- seq(from = min(cti_AB$LSI), to = max(cti_AB$LSI), length.out = 20)

res_class_PLAND <- c()
for(i in PLAND_val){
  for (j in LSI_val){
    p <- visreg(m_frag, xvar = "Year", cond = list(PLAND = i, LSI = j), plot = F)
    slope <- coefficients(lm(visregFit ~ Year, data = p$fit))[2]
    slopeLw <- coefficients(lm(visregLwr ~ Year, data = p$fit))[2]
    slopeUp <- coefficients(lm(visregUpr ~ Year, data = p$fit))[2]
    res_class_PLAND <- bind_rows(res_class_PLAND, tibble(PLAND = i, LSI = j, slope = slope, slopeLw = slopeLw, slopeUp = slopeUp))
  }
}

ggplot(data = res_class_PLAND, aes(x = LSI, y = slope, color = as.factor(PLAND))) + geom_line() + 
  geom_ribbon(aes(ymin=slopeLw, ymax=slopeUp, fill = as.factor(PLAND)), alpha=0.1, colour=NA) +
  scale_y_continuous("Temporal trend of CTI") + 
  scale_x_continuous("Fragmentation") + 
  scale_color_manual("Area of SNH", labels = PLAND_names, values=col.PLAND) +
  scale_fill_manual("Area of SNH", labels = PLAND_names, values=col.PLAND)

