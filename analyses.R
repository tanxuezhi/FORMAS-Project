library(rio)
library(lme4)
library(lmerTest)
library(ggplot2)
library(cowplot)
library(visreg)
library(MuMIn)

##### import CTI data and merge with site data #####

sites <- rio::import(file = "../Data/Butterflies - Netherlands/Sites.xlsx", which = 1L)
ndvi_data <- read.csv(paste0("../Data/NDVI/NDVI_NL/ndvi.buffers.sum.scv"))
lc_class <- rio::import(file = "../Landcover/clc_legend.xls", which = 1L)

cti_AB <- read.csv2("../Data/Butterflies - Netherlands/cti_site_year_abundance_2017-05-19.csv", dec = ".")
cti_AB$Site <- as.factor(cti_AB$Site)
cti_AB <- merge(cti_AB, ndvi_data)
cti_AB <- merge(cti_AB, sites)
cti_AB <- merge(cti_AB, lc_class, by.x = "Landcover", by.y = "CLC_CODE")
cti_AB_mean <- aggregate(cti ~ Year, cti_AB, mean)


##### try first preliminary analyses #####

## cti change over time
m_AB <- lmer(cti ~ Year + (1|LABEL2/Site), data = cti_AB)
summary(m_AB)
visreg(m_AB, xvar = "Year", ylab = "CTI", scale = "response", ylim = c(8.9,9.3))
points(cti ~ Year, cti_AB_mean, pch = 16, type = "b")


## effect of NDVI variability
# long-term variability
m_ltv_AB <- lmer(cti ~ Year * sd.LT_ndvi + (1|LABEL2/Site), data = cti_AB)
summary(m_ltv_AB)
pred.m_ltv_AB <- visreg(m_ltv_AB, xvar = "Year", ylab = "CTI", scale = "response", 
                        by = "sd.LT_ndvi", gg = T, breaks = 3, plot = F)
p.ltv_AB <- ggplot(data = pred.m_ltv_AB$fit, aes(x = Year, y = visregFit, color = as.factor(sd.LT_ndvi))) + 
  geom_line() + 
  scale_y_continuous(name = "CTI") +
  theme(legend.position = c(0.8,0.1)) + 
  scale_color_manual(name = "NDVI variability", labels = c("Low", "Medium", "High"), values = c("royalblue1", "palegreen3", "tomato2")) 


# short-term variability
m_stv_AB <- lmer(cti ~ Year * CV.ST_ndvi + (1|LABEL2/Site), data = cti_AB)
summary(m_stv_AB)
pred.m_ltv_AB <- visreg(m_stv_AB, xvar = "Year", ylab = "CTI", scale = "response", 
                        by = "CV.ST_ndvi", gg = T, breaks = 3, plot = F)
p.stv_AB <- ggplot(data = pred.m_ltv_AB$fit, aes(x = Year, y = visregFit, color = as.factor(CV.ST_ndvi))) + 
  geom_line() + 
  scale_y_continuous(name = "CTI") +
  theme(legend.position = c(0.8,0.1)) + 
  scale_color_manual(name = "NDVI variability", labels = c("Low", "Medium", "High"), values = c("royalblue1", "palegreen3", "tomato2")) 
