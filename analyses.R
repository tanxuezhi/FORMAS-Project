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
frag_data_NL <- read.csv("../Connectivity/Fragmentation/NL/Frag_indices.csv")
frag_data_FIN <- read.csv("../Connectivity/Fragmentation/FIN/Frag_indices.csv")

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
# response variable  = cti
vif(cti_AB[,c("NLSI","PLAND","Year")])
cti_AB$LABEL3 <- as.factor(cti_AB$LABEL3)
cti_AB[is.na(cti_AB$NLSI),"NLSI"] <- 0

# cti_AB <- subset(cti_AB, cti_AB$Scale > 2000)

scaleTest <- c()
for(i in unique(cti_AB$Scale)){
  m_frag_std <- lmer(cti ~ NLSI * CA * Year + country + LABEL3 + (1|Site), data = stdize(subset(cti_AB, cti_AB$Scale == i), prefix = F))
  scaleTest <- rbind.data.frame(scaleTest, cbind.data.frame(Scale = i, AICc = AICc(m_frag_std), coefInter = fixef(m_frag_std)["NLSI:CA:Year"]))
}
plot(scaleTest[order(scaleTest$Scale),], xlab = "Spatial scale (m)", ylab = "AICc", type = "o", pch = 16)
plot(scaleTest[order(scaleTest$Scale),-2], xlab = "Spatial scale (m)", ylab = "coef", type = "o", pch = 16)


cti_AB_sel <- cti_AB[cti_AB$Scale == scaleTest[order(scaleTest$AICc), "Scale"][1],]

m_frag <- lmer(cti ~ NLSI * PLAND * Year + LABEL3 + country + (1|Site), data = cti_AB_sel, na.action = na.fail)
m_frag_std <- lmer(cti ~ NLSI * CA * Year + country + LABEL3 + (1|Site), data = stdize(cti_AB_sel, prefix = F), na.action = na.fail)
anova(m_frag_std)
summary(m_frag_std)

ncf:::spline.correlog(cti_AB$X, cti_AB$Y, resid(m_frag_std, type  = "pearson"),
                      resamp=100, na.rm = T)

# plot three-way interaction
# 2d plots
PLAND_val <- seq(from = min(cti_AB$PLAND), to = max(cti_AB$PLAND), length.out = 20)
NLSI_val <- seq(from = min(cti_AB$NLSI), to = max(cti_AB$NLSI), length.out = 20)

res <- c()
n.max <- length(CA_val) * length(NLSI_val)
n <- 0
for(i in PLAND_val){
  for (j in NLSI_val){
    n = n+1
    cat(paste(round(n/n.max*100, 2), "%, parameters :", "PLAND =", round(i,2), ", nLSI =", round(j,2), "\n"))
    p <- visreg(m_frag, xvar = "Year", cond = list(PLAND = i, NLSI = j), plot = F)
    slope <- coefficients(lm(visregFit ~ Year, data = p$fit))[2]
    res <- bind_rows(res, tibble(PLAND = i/100, NLSI = j, slope = slope))
  }
}

quilt.plot(res,  nx = length(PLAND_val), ny = length(NLSI_val), 
           xlab = "SNH Area", ylab = "Landscape Shape Index",
           main = "Temporal trend of CTI")
points(NLSI ~ PLAND, data = cti_AB, pch = 16, col = country)

# grouped by classes (both factors)
PLAND_val <- c(mean(cti_AB$PLAND)-sd(cti_AB$PLAND), mean(cti_AB$PLAND), mean(cti_AB$PLAND)+sd(cti_AB$PLAND))
NLSI_val <- c(mean(cti_AB$NLSI)-sd(cti_AB$NLSI), mean(cti_AB$NLSI), mean(cti_AB$NLSI)+sd(cti_AB$NLSI))

res_class <- c()
for(i in PLAND_val){
  for (j in NLSI_val){
    p <- visreg(m_frag, xvar = "Year", cond = list(PLAND = i, NLSI = j), plot = F)
    res_class <- bind_rows(res_class, p$fit)
  }
}

# version 1
PLAND_names <- c("Small","Medium","Large"); names(PLAND_names) <- PLAND_val
NLSI_names <- c("Little fragmentated","Moderately fragmentated","Highly fragmentated"); names(NLSI_names) <- NLSI_val
names_fac <- c(PLAND_names, NLSI_names)

ggplot(data = res_class, aes(x = Year, y = visregFit)) + geom_line() + 
  facet_grid(PLAND ~ NLSI, labeller = as_labeller(names_fac)) + 
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.2) +
  scale_y_continuous("Community temperature Index")


# version 2
col.PLAND <- c("dodgerblue", "gold2", "firebrick"); names(col.PLAND) <- PLAND_val

ggplot(data = res_class, aes(x = Year, y = visregFit, color = as.factor(PLAND))) + geom_line() + 
  facet_grid( ~ NLSI, labeller = as_labeller(NLSI_names)) + 
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr, fill = as.factor(PLAND)), alpha=0.1, colour=NA) +
  scale_y_continuous("Community temperature Index") +
  scale_color_manual("Area of SNH", labels = names_fac, values=col.PLAND) +
  scale_fill_manual("Area of SNH", labels = names_fac, values=col.PLAND)



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

