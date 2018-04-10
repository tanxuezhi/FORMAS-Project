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

###############################
########## Load data ########## 
###############################

butterflies.data.presence <- as.tbl(fread("../Data/cti_butterflies_data.csv")) %>% dplyr::filter(type == "Presence")
butterflies.data.presence <- butterflies.data.presence %>% 
  mutate(gridCell50 = ifelse(country == "NL", paste0(gridCell50, "_NL"), paste0(gridCell50, "_FIN")))

butterflies.data.abundance <- as.tbl(fread("../Data/cti_butterflies_data.csv")) %>% dplyr::filter(type == "Abundance")
butterflies.data.abundance <- butterflies.data.abundance %>% 
  mutate(gridCell50 = ifelse(country == "NL", paste0(gridCell50, "_NL"), paste0(gridCell50, "_FIN")))

ggplot(butterflies.data.presence %>% dplyr:::filter(PLAND < .95 & PLAND > .05) %>% group_by(country, Scale, Habitat) %>%
         distinct(), aes(y = CLUMPY, x = PLAND, color =country)) + 
  geom_point(alpha = .3) + facet_wrap( ~ Scale, labeller = label_wrap) + 
  scale_color_manual("Country", values = c("#00008B", "#BDB76B"), labels = c("Finland", "Netherlands")) + 
  scale_y_continuous("Habitat clumpiness") + scale_x_continuous("% Semi-natural habitat")


#########################################
## prepare data for landscape analysis ##
#########################################

## select only sites with .05 < %SNH < .95
butterflies.data.presence <- butterflies.data.presence %>% dplyr:::filter(PLAND <= .95 & PLAND >= .05)
butterflies.data.abundance <- butterflies.data.abundance %>% dplyr:::filter(PLAND <= .95 & PLAND >= .05)

## standardize ##
butterflies.data.presence <- stdize(butterflies.data.presence, prefix = F, omit.cols= c("cti", "gridCell50", "Scale"))
butterflies.data.abundance <- stdize(butterflies.data.abundance, prefix = F, omit.cols= c("cti", "gridCell50", "Scale"))
scaleList <- list(scale = attr(butterflies.data.abundance, "scaled:scale")[c("Year", "CLUMPY", "PLAND")],
                  center = attr(butterflies.data.abundance, "scaled:center")[c("Year", "CLUMPY", "PLAND")])

## select data 
dat <- butterflies.data.presence

##############
## Analyses ##
##############

## run models for all scales and plot coefficients
mods <- foreach(i = unique(dat$Scale), .combine = rbind) %do% {
  cat(paste("Scale = ", i/1000, "km\n"))
  m <- lme(cti ~ CLUMPY * PLAND * Year + Habitat * Year, random = list(gridCell50 = ~ 1, Site = ~ 1), 
           correlation = corExp(form = ~ X + Y), 
           data = dat %>% dplyr:::filter(Scale == i) %>% 
             mutate(X = jitter(X), Y = jitter(Y)))
  cbind.data.frame(Scale = i, intervals(m, which = "fixed")$fixed) %>% 
    rownames_to_column("term")
}


mods <- mods %>% dplyr::filter(term %in% c("CLUMPY:Year", "PLAND:Year", "CLUMPY:PLAND:Year"))
mods$term <- gsub("CLUMPY:PLAND:Year","% Semi-natural habitat\n x Habitat clumpiness", mods$term)
mods$term <- gsub("PLAND:Year","% Semi-natural habitat", mods$term)
mods$term <- gsub("CLUMPY:Year","Habitat clumpiness", mods$term)
mods$term <- factor(mods$term, levels = c("% Semi-natural habitat", "Habitat clumpiness","% Semi-natural habitat\n x Habitat clumpiness"))
mods$Scale <- factor(mods$Scale, levels = c(1,3,5,10,20,30,50)*1000)
levels(mods$Scale) <- paste(c(1,3,5,10,20,30,50), "km")

ggplot(mods,aes(x = Scale, y = est.)) + 
  geom_point() + 
  facet_grid(~ term) + 
  geom_errorbar(aes(x = Scale, ymin = lower, ymax = upper), width = 0) + 
  geom_hline(yintercept = 0, linetype = 2) + 
  scale_y_continuous("Coefficient +/- 95%CI")


######################
#### plot results ####
######################
registerDoFuture()
options(future.globals.maxSize = Inf)
plan(multiprocess, workers = 2)

## run models and list them

mods.out <- foreach(i = unique(dat$Scale)) %do% {
  cat(paste("Scale = ", i/1000, "km\n"))
  lme(cti ~ CLUMPY * PLAND * Year + Habitat * Year, random = list(gridCell50 = ~ 1, Site = ~ 1), 
      correlation = corExp(form = ~ X + Y), 
      data = dat %>% dplyr:::filter(Scale == i) %>% 
        mutate(X = jitter(X), Y = jitter(Y)))
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
  scale_x_continuous("% Semi-natural habitat") + scale_y_continuous("Habitat clumpiness")


###########################
## remove unused objects ##
###########################


gdata:::keep(butterflies.data.presence, butterflies.data.abundance,
             mods, mods.out, 
             pred,
             sure = T)


