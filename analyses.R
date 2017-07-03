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
frag_data <- read.csv("../Connectivity/Fragmentation/NL/Frag_indices_NL.csv")

sites <- merge(sites, frag_data, by.x = "Site", by.y = "LID")


cti_AB <- read.csv2("../Data/Butterflies - Netherlands/cti_site_year_abundance_2017-05-19.csv", dec = ".")
cti_AB$Site <- as.factor(cti_AB$Site)
cti_AB <- merge(cti_AB, sites, by.x = "Site", by.y = "Site")
cti_AB_mean <- aggregate(cti ~ Year, cti_AB, mean)


##### try first preliminary analyses #####

## cti change over time
m_AB <- lmer(cti ~ Year + (1|LABEL3/Site), data = cti_AB)
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


##### Effect of fragmentation #####
# response variable  = slope of cti ~ year
slope_year_site <- c()
for(i in unique(cti_AB$Site)){
  m_slope_year <- lm(cti ~ Year, data = cti_AB[cti_AB$Site %in% i,])
  slope <- cbind.data.frame(Site = i, slope = coefficients(m_slope_year)[2], nyears = dim(m_slope_year$model)[1])
  rownames(slope) <- NULL
  slope_year_site <- rbind.data.frame(slope_year_site, slope)
}

sites <- merge(sites, slope_year_site)

m_frag_slope <- lm(slope ~ LSI * CA + LABEL3, weight =nyears, data = na.omit(sites[sites$nyears > 10,]), na.action = na.fail)
anova(m_frag_slope)
summary(m_frag_slope)

visreg2d(m_frag_slope, xvar = "CA", yvar = "LSI", scale = "response")

# response variable  = cti
vif(cti_AB[,c("NLSI","PROX_MN","Year")])


cti_AB$LABEL3 <- as.factor(cti_AB$LABEL3)
m_frag <- lmer(cti ~ NLSI * PROX_MN * Year + LABEL3 + (1|Site), data = cti_AB, na.action = na.fail)
m_frag_std <- lmer(cti ~ NLSI * PROX_MN * Year + LABEL3 + (1|Site), data = stdize(cti_AB, prefix = F), na.action = na.fail)
anova(m_frag_std)
summary(m_frag_std)

dredge(m_frag_std)

m_frag_reduced <- lmer(cti ~ PROX_MN * Year + LABEL3 + (1|Site), data = cti_AB, na.action = na.fail)
anova(m_frag_reduced)
summary(m_frag_reduced)
visreg(m_frag_reduced, xvar = "Year", by = "PROX_MN", scale = "response")


PROX_MN_val <- seq(from = min(cti_AB$PROX_MN), to = max(cti_AB$PROX_MN), length.out = 20)
NLSI_val <- seq(from = min(cti_AB$NLSI), to = max(cti_AB$NLSI), length.out = 20)

res <- c()
n.max <- length(CA_val) * length(NLSI_val)
n <- 0
for(i in PROX_MN_val){
  for (j in NLSI_val){
    n = n+1
    cat(paste(round(n/n.max*100, 2), "%, parameters :", "PROX_MN =", round(i,2), ", LSI =", round(j,2), "\n"))
    p <- visreg(m_frag, xvar = "Year", cond = list(PROX_MN = i, NLSI = j), plot = F)
    slope <- coefficients(lm(visregFit ~ Year, data = p$fit))[2]
    res <- bind_rows(res, tibble(PROX_MN = i, NLSI = j, slope = slope))
  }
}

quilt.plot(res,  nx = length(PROX_MN_val), ny = length(NLSI_val), xlab = "Proximity index", ylab = "Normalized Landscape Shape Index")
points(NLSI ~ PROX_MN, data = cti_AB)


