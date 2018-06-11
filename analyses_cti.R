library(data.table)
library(lmerTest)
library(cowplot)
library(MuMIn)
library(visreg)
library(raster)
source("functions.R")
library(broom)
library(foreach)
library(doFuture)
library(nlme)
library(tidyverse)

###############################
########## Load data ########## 
###############################

butterflies.data.presence <- as.tbl(fread("../Data/cti_butterflies_data.csv")) %>% dplyr::filter(type == "Presence")
butterflies.data.presence <- butterflies.data.presence %>% 
  mutate(gridCell50 = ifelse(country == "NL", paste0(gridCell50, "_NL"), paste0(gridCell50, "_FIN")))

butterflies.data.abundance <- as.tbl(fread("../Data/cti_butterflies_data.csv")) %>% dplyr::filter(type == "Abundance")
butterflies.data.abundance <- butterflies.data.abundance %>% 
  mutate(gridCell50 = ifelse(country == "NL", paste0(gridCell50, "_NL"), paste0(gridCell50, "_FIN")))


#########################################
## prepare data for landscape analysis ##
#########################################

## select only sites with .05 < %SNH < .95
# butterflies.data.presence <- butterflies.data.presence %>% dplyr:::filter(PLAND <= .95 & PLAND >= .05)
# butterflies.data.abundance <- butterflies.data.abundance %>% dplyr:::filter(PLAND <= .95 & PLAND >= .05)

## standardize ##
butterflies.data.presence <- stdize(butterflies.data.presence, prefix = F, omit.cols= c("cti", "gridCell50", "Scale"))
butterflies.data.abundance <- stdize(butterflies.data.abundance, prefix = F, omit.cols= c("cti", "gridCell50", "Scale"))
scaleList <- list(scale = attr(butterflies.data.presence, "scaled:scale")[c("Year", "CLUMPY", "PLAND")],
                  center = attr(butterflies.data.presence, "scaled:center")[c("Year", "CLUMPY", "PLAND")])

## select data 
dat <- butterflies.data.presence

##############
## Analyses ##
##############

## run models for all scales and plot coefficients
mods <- foreach(i = unique(dat$Scale), .combine = rbind) %do% {
  cat(paste("Scale = ", i/1000, "km\n"))
  m <- lmer(cti ~ CLUMPY * PLAND * Year + Habitat * Year + X*Y + (1|gridCell50/Site), 
           data = dat %>% dplyr:::filter(Scale == i))
  cbind.data.frame(Scale = i, est. = fixef(m), confint(m, method = "Wald")[-c(1:3),]) %>% 
    rownames_to_column("term")
}
colnames(mods)[4:5] <- c("lower", "upper")

mods <- mods %>% dplyr::filter(term %in% c("CLUMPY:Year", "PLAND:Year", "CLUMPY:PLAND:Year"))
mods$term <- gsub("CLUMPY:PLAND:Year","Proportion × Aggregation of SNH", mods$term)
mods$term <- gsub("PLAND:Year","Proportion of SNH", mods$term)
mods$term <- gsub("CLUMPY:Year","Aggregation of SNH", mods$term)
mods$term <- factor(mods$term, levels = c("Proportion of SNH", "Aggregation of SNH","Proportion × Aggregation of SNH"))
mods$Scale <- factor(mods$Scale, levels = c(1,3,5,10,20,30,50)*1000)
levels(mods$Scale) <- paste(c(1,3,5,10,20,30,50), "km")

ggplot(mods,aes(x = Scale, y = est.)) + 
  geom_point() + 
  facet_grid( ~term) + 
  geom_errorbar(aes(x = Scale, ymin = lower, ymax = upper), width = 0) + 
  geom_hline(yintercept = 0, linetype = 2) + 
  scale_y_continuous("Coefficient +/- 95%CI") +
  theme(axis.text.x=element_text(angle=90,hjust=1), strip.background = element_blank())
  

ggsave("../Coef_cti.svg", height = 3, width = 8, dpi = 1200)

######################
#### plot results ####
######################
registerDoFuture()
options(future.globals.maxSize = Inf)
plan(multiprocess, workers = 2)

## run models and list them

mods.out <- foreach(i = unique(dat$Scale)) %do% {
  cat(paste("Scale = ", i/1000, "km\n"))
  lmer(cti ~ CLUMPY * PLAND * Year + Habitat * Year + X*Y + (1|gridCell50/Site), 
            data = dat %>% dplyr:::filter(Scale == i))
}


## make prediction dataframe

pred <- foreach(i = 1:length(mods.out), j = unique(dat$Scale), .combine = rbind) %dopar% {
  pred.temp <- predict_raster(mods.out[[i]], scaleList, 100)
  pred.temp <- cbind.data.frame(Scale = j, pred.temp)
}

## and plot it...
ggplot(pred, aes(x= x, y = y, fill = estimate)) + geom_raster() + facet_wrap(~Scale, labeller = label_wrap) +
  scale_fill_gradientn(name = "CTI trend\n", 
                       breaks =  c(min(pred$estimate),0,max(pred$estimate)),
                       labels = c(round(min(pred$estimate), 3),"0.000",round(max(pred$estimate), 3)),
                       limits = c(min(pred$estimate),max(pred$estimate)),
                       colours=c("light blue", "yellow","red")) +
  scale_x_continuous("Proportion of semi-natural habitat") + scale_y_continuous("Habitat clumpiness")


###########################
## remove unused objects ##
###########################


gdata:::keep(butterflies.data.presence, butterflies.data.abundance,
             mods, mods.out, 
             pred,
             sure = T)


