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

##### import CTI data and merge with site data #####

# sites data
sites_NL <- read.csv(file = "../Data/Butterflies - Netherlands/Sites_NL_ETRS89_landcover.csv")
sites_FIN <- read.csv(file = "../Data/Butterflies - Finland/Sites_FIN_ETRS89_landcover.csv")

lc_class <- rio::import(file = "../Landcover/clc_legend.xls", which = 1L)

# ndvi_data <- read.csv("../Data/NDVI/NDVI_NL/ndvi.points.sum.scv")
frag_data_NL <- read.csv("../Connectivity/Fragmentation/NL/Frag_indices_NL.csv")
frag_data_FIN <- read.csv("../Connectivity/Fragmentation/FIN/Frag_indices_FIN.csv")

sites_NL <- merge(sites_NL, frag_data_NL, by.x = "Site", by.y = "LID")
sites_NL <- merge(sites_NL, lc_class, by.x = "Landcover", by.y = "CLC_CODE")
sites_NL <- sites_NL[,-which(names(sites_NL) %in% "GRID_CODE")]

sites_FIN <- merge(sites_FIN, frag_data_FIN, by.x = "Site", by.y = "LID")
sites_FIN <- merge(sites_FIN, lc_class, by.x = "Landcover", by.y = "GRID_CODE")
sites_FIN <- sites_FIN[,-which(names(sites_FIN) %in% "CLC_CODE")]

# cti data
cti_AB_NL <- read.csv2("../Data/Butterflies - Netherlands/cti_site_year_abundance_2017-05-19.csv", dec = ".")[,-3]
cti_AB_NL <- merge(cti_AB_NL, sites_NL, by.x = "Site", by.y = "Site")
cti_AB_NL_mean <- aggregate(cti ~ Year, cti_AB_NL, mean)

cti_AB_FIN <- read.csv("../Data/Butterflies - Finland/CTI_Abundance_FINLAND_1999-2016.csv", dec = ".")
cti_AB_FIN <- merge(cti_AB_FIN, sites_FIN, by.x = "Site", by.y = "Site")
cti_AB_FIN_mean <- aggregate(cti ~ Year, cti_AB_FIN, mean)

cti_AB <- rbind(cbind.data.frame(cti_AB_NL, country = "NL"), cbind.data.frame(cti_AB_FIN, country = "FIN"))
cti_AB$Site <- paste0(cti_AB$Site, "_", cti_AB$country)

##### try first preliminary analyses #####

## cti change over time
m_AB <- lmer(cti ~ Year*country + LABEL3 + (1|Site), data = cti_AB)
summary(m_AB)
anova(m_AB)
visreg(m_AB, xvar = "Year", by = "country", ylab = "CTI", scale = "response")

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
vif(cti_AB[,c("NLSI","CA","Year")])

cti_AB$LABEL3 <- as.factor(cti_AB$LABEL3)
m_frag <- lmer(cti ~ NLSI * CA * Year + LABEL3 + country + (1|Site), data = cti_AB, na.action = na.fail)
m_frag_std <- lmer(cti ~ NLSI * CA * Year + country + LABEL3 + (1|Site), data = stdize(cti_AB, prefix = F), na.action = na.fail)
anova(m_frag_std)
summary(m_frag_std)

ncf:::spline.correlog(cti_AB$X, cti_AB$Y, resid(m_frag_std, type  = "pearson"),
                      resamp=100, na.rm = T)

dredge(m_frag_std, fixed = c("LABEL3","country","Year"))
m_frag_reduced_std <- lmer(cti ~ NLSI * Year + CA * Year + LABEL3 + country + (1|Site), data = stdize(cti_AB, prefix = F), na.action = na.fail)
anova(m_frag_reduced_std)
summary(m_frag_reduced_std)


# plot three-way interaction
# 2d plots
CA_val <- seq(from = min(cti_AB$CA), to = max(cti_AB$CA), length.out = 20)
NLSI_val <- seq(from = min(cti_AB$NLSI), to = max(cti_AB$NLSI), length.out = 20)

res <- c()
n.max <- length(CA_val) * length(NLSI_val)
n <- 0
for(i in CA_val){
  for (j in NLSI_val){
    n = n+1
    cat(paste(round(n/n.max*100, 2), "%, parameters :", "CA =", round(i,2), ", nLSI =", round(j,2), "\n"))
    p <- visreg(m_frag, xvar = "Year", cond = list(CA = i, NLSI = j), plot = F)
    slope <- coefficients(lm(visregFit ~ Year, data = p$fit))[2]
    res <- bind_rows(res, tibble(CA = i, NLSI = j, slope = slope))
  }
}

quilt.plot(res,  nx = length(CA_val), ny = length(NLSI_val), 
           xlab = "SNH Area", ylab = "Landscape Shape Index",
           main = "Temporal trend of CTI")
points(NLSI ~ CA, data = cti_AB, col = country, pch = 16)

# grouped by classes (both factors)
CA_val <- c(mean(cti_AB$CA)-sd(cti_AB$CA), mean(cti_AB$CA), mean(cti_AB$CA)+sd(cti_AB$CA))
NLSI_val <- c(mean(cti_AB$NLSI)-sd(cti_AB$NLSI), mean(cti_AB$NLSI), mean(cti_AB$NLSI)+sd(cti_AB$NLSI))

res_class <- c()
for(i in CA_val){
  for (j in NLSI_val){
    p <- visreg(m_frag, xvar = "Year", cond = list(CA = i, NLSI = j), plot = F)
    res_class <- bind_rows(res_class, p$fit)
  }
}

# version 1
CA_names <- c("Small","Medium","Large"); names(CA_names) <- CA_val
NLSI_names <- c("Little fragmentated","Moderately fragmentated","Highly fragmentated"); names(NLSI_names) <- NLSI_val
names_fac <- c(CA_names, NLSI_names)

ggplot(data = res_class, aes(x = Year, y = visregFit)) + geom_line() + 
  facet_grid(CA ~ NLSI, labeller = as_labeller(names_fac)) + 
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.2) +
  scale_y_continuous("Community temperature Index")


# version 2
col.CA <- c("dodgerblue", "gold2", "firebrick"); names(col.CA) <- CA_val

ggplot(data = res_class, aes(x = Year, y = visregFit, color = as.factor(CA))) + geom_line() + 
  facet_grid( ~ NLSI, labeller = as_labeller(NLSI_names)) + 
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr, fill = as.factor(CA)), alpha=0.1, colour=NA) +
  scale_y_continuous("Community temperature Index") +
  scale_color_manual("Area of SNH", labels = names_fac, values=col.CA) +
  scale_fill_manual("Area of SNH", labels = names_fac, values=col.CA)



# grouped by classes (not fragmentation)
CA_val <- c(mean(cti_AB$CA)-sd(cti_AB$CA), mean(cti_AB$CA), mean(cti_AB$CA)+sd(cti_AB$CA))
NLSI_val <- seq(from = min(cti_AB$NLSI), to = max(cti_AB$NLSI), length.out = 20)

res_class_CA <- c()
for(i in CA_val){
  for (j in NLSI_val){
    p <- visreg(m_frag, xvar = "Year", cond = list(CA = i, NLSI = j), plot = F)
    slope <- coefficients(lm(visregFit ~ Year, data = p$fit))[2]
    slopeLw <- coefficients(lm(visregLwr ~ Year, data = p$fit))[2]
    slopeUp <- coefficients(lm(visregUpr ~ Year, data = p$fit))[2]
    res_class_CA <- bind_rows(res_class_CA, tibble(CA = i, NLSI = j, slope = slope, slopeLw = slopeLw, slopeUp = slopeUp))
  }
}

ggplot(data = res_class_CA, aes(x = NLSI, y = slope, color = as.factor(CA))) + geom_line() + 
  geom_ribbon(aes(ymin=slopeLw, ymax=slopeUp, fill = as.factor(CA)), alpha=0.1, colour=NA) +
  scale_y_continuous("Temporal trend of CTI") + 
  scale_x_continuous("Fragmentation") + 
  scale_color_manual("Area of SNH", labels = CA_names, values=col.CA)+
  scale_fill_manual("Area of SNH", labels = CA_names, values=col.CA)

