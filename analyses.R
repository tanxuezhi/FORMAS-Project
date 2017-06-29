library(rio)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(cowplot)
library(visreg)
library(MuMIn)

##### import CTI data and merge with site data #####
#load sites
sites <- read.csv(file = "../Data/Butterflies - Netherlands/Sites_ETRS89_landcover.csv")
lc_class <- rio::import(file = "../Landcover/clc_legend.xls", which = 1L)
sites <- merge(sites, lc_class, by.x = "Landcover", by.y = "CLC_CODE")

# site data
ndvi_data <- read.csv("../Data/NDVI/NDVI_NL/ndvi.points.sum.scv")
frag_data <- read.csv("../Data/Butterflies - Netherlands/landFrag.csv")

cti_AB <- read.csv2("../Data/Butterflies - Netherlands/cti_site_year_abundance_2017-05-19.csv", dec = ".")
cti_AB$Site <- as.factor(cti_AB$Site)
cti_AB <- merge(cti_AB, frag_data)
cti_AB <- merge(cti_AB, sites)
cti_AB_mean <- aggregate(cti ~ Year, cti_AB, mean)


##### try first preliminary analyses #####

## cti change over time
m_AB <- lmer(cti ~ Year + (1|LABEL2/Site), data = cti_AB)
summary(m_AB)
visreg(m_AB, xvar = "Year", ylab = "CTI", scale = "response", ylim = c(8.9,9.3))
points(cti ~ Year, cti_AB_mean, pch = 16, type = "b")


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


##### Effect of connectivity #####
m_frag_AB <- lmer(cti ~ Year * patch.cohesion.index	
 + LABEL2 + (1|Site), data = stdize(cti_AB, prefix = F), na.action = na.fail)
summary(m_frag_AB)

dredge(m_frag_AB, fixed = c("Year", "LABEL3"))

# plot
pred.m_frag_AB <- visreg(m_frag_AB, xvar = "Year", ylab = "CTI", scale = "response", 
                        by = "edge.density", gg = T, breaks = 3, plot = F)
ggplot(data = pred.m_frag_AB$fit, aes(x = Year, y = visregFit, color = as.factor(edge.density))) + 
  geom_line() + 
  scale_y_continuous(name = "CTI") +
  theme(legend.position = c(0.8,0.1)) + 
  scale_color_manual(name = "edge.density", labels = c("Low", "Medium", "High"), values = c("royalblue1", "palegreen3", "tomato2")) 


