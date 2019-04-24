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
library(emmeans)

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
butterflies.data.presence <- stdize(butterflies.data.presence, prefix = F, 
                                    omit.cols= c("cti", "Year", "gridCell50", "Scale", "temp"))
butterflies.data.abundance <- stdize(butterflies.data.abundance, prefix = F, 
                                     omit.cols= c("cti", "Year", "gridCell50", "Scale", "temp"))
scaleList <- list(scale = attr(butterflies.data.presence, "scaled:scale")[c("CLUMPY", "PLAND")],
                  center = attr(butterflies.data.presence, "scaled:center")[c("CLUMPY", "PLAND")])

## select data 
dat <- butterflies.data.presence %>%
  mutate(Year = Year - 1990)


##############
## Analyses ##
##############

## run models for all scales and plot coefficients
mods <- foreach(i = unique(dat$Scale), .combine = rbind) %do% {
  cat(paste("Scale = ", i/1000, "km\n"))
  m <- lmer(cti ~ CLUMPY * PLAND * Year + Habitat * Year + X*Y + (Year|gridCell50/Site), 
            data = dat %>% dplyr:::filter(Scale == i))
  cbind.data.frame(Scale = i, est. = fixef(m), confint(m, method = "Wald")[-c(1:7),]) %>% 
    rownames_to_column("term")
}
colnames(mods)[4:5] <- c("lower", "upper")

mods <- mods %>% dplyr::filter(term %in% c("CLUMPY:Year", "PLAND:Year", "CLUMPY:PLAND:Year"))
mods$term <- gsub("CLUMPY:PLAND:Year","Proportion × Clumpiness of SNH", mods$term)
mods$term <- gsub("PLAND:Year","Proportion of SNH", mods$term)
mods$term <- gsub("CLUMPY:Year","Clumpiness of SNH", mods$term)
mods$term <- factor(mods$term, levels = c("Proportion of SNH", "Clumpiness of SNH","Proportion × Clumpiness of SNH"))
mods$Scale <- factor(mods$Scale, levels = c(1,3,5,10,20,30,50)*1000)
levels(mods$Scale) <- paste(c(1,3,5,10,20,30,50), "km")

ggplot(mods,aes(x = Scale, y = est.)) + 
  geom_point() + 
  facet_grid( ~term) + 
  geom_errorbar(aes(x = Scale, ymin = lower, ymax = upper), width = 0) + 
  geom_hline(yintercept = 0, linetype = 2) + 
  scale_y_continuous("Coefficient +/- 95%CI") +
  theme(axis.text.x=element_text(angle=90,hjust=1), strip.background = element_blank())

ggplot(mods %>% mutate(Significancy = ifelse(lower < 0 & upper < 0 | lower > 0 & upper > 0, "Significant", "Not significant"),
                       term = factor(mods$term, levels = rev(levels(mods$term)))), 
       aes(x = term, y = est., color = Scale)) + 
  geom_point(position = position_dodge(.8)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(.8)) + 
  geom_hline(yintercept = 0, lty = 2, color = "grey") +
  coord_flip() + 
  scale_color_manual(values = colorRampPalette(c("grey", "black"))(7)) + 
  theme_bw()


ggsave("../Coef_cti.svg", height = 3, width = 8, dpi = 1200)

### alt ###
# compare model coefficients

mods_all <- rbind.data.frame(
  cbind.data.frame(Landcover_classification = "Site-specific", 
                   rbind.data.frame(cbind.data.frame(Country = "Finland", type = "Presence", mods_sp_presence_FIN),
                                    cbind.data.frame(Country = "Finland", type = "Abundance", mods_sp_abundance_FIN),
                                    cbind.data.frame(Country = "Netherlands", type = "Presence", mods_sp_presence_NL),
                                    cbind.data.frame(Country = "Netherlands", type = "Abundance", mods_sp_abundance_NL),
                                    cbind.data.frame(Country = "All", type = "Presence", mods_sp_presence),
                                    cbind.data.frame(Country = "All", type = "Abundance", mods_sp_abundance)
                   )
  ),
  cbind.data.frame(Landcover_classification = "Generalist", 
                   rbind.data.frame(cbind.data.frame(Country = "Finland", type = "Presence", mods_gen_presence_FIN),
                                    cbind.data.frame(Country = "Finland", type = "Abundance", mods_gen_abundance_FIN),
                                    cbind.data.frame(Country = "Netherlands", type = "Presence", mods_gen_presence_NL),
                                    cbind.data.frame(Country = "Netherlands", type = "Abundance", mods_gen_abundance_NL),
                                    cbind.data.frame(Country = "All", type = "Presence", mods_gen_presence),
                                    cbind.data.frame(Country = "All", type = "Abundance", mods_gen_abundance)
                   )
  )
)

write_csv(mods_all, "../mods_cti_comp.csv")

ggplot(mods_all %>% mutate(term = factor(term, levels = rev(levels(term)))) %>%
         filter(Landcover_classification == "Generalist"), 
       aes(x = term, y = est., color = Scale)) + 
  geom_point(position = position_dodge(.8)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = position_dodge(.8)) + 
  geom_hline(yintercept = 0, lty = 2, color = "grey") +
  coord_flip() + 
  facet_grid(Country ~ type) + 
  scale_color_manual(values = colorRampPalette(c("grey", "black"))(7)) + 
  theme_bw()


######################
#### plot results ####
######################
registerDoFuture()
options(future.globals.maxSize = Inf)
plan(multiprocess, workers = 2)

## run models and list them

mods.out <- foreach(i = unique(dat$Scale)) %dopar% {
  cat(paste("Scale = ", i/1000, "km\n"))
  lmer(cti ~ CLUMPY * PLAND * Year + Habitat * Year + X*Y + (Year|gridCell50/Site), 
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


############
# alt plot #
############

mod20 <-   lmer(cti ~ CLUMPY * PLAND * Year + Habitat * Year + X*Y + (Year|gridCell50/Site), 
                data = dat %>% dplyr:::filter(Scale == 20000)) # do not standardize !!


percentiles_percountry <- dat %>% dplyr:::filter(Scale == 20000) %>% group_by(Site) %>% 
  summarize_all(first) %>% dplyr::select(1,6,7,13,15,16,17) %>% 
  mutate(PLAND = PLAND * scaleList$scale["PLAND"] + scaleList$center["PLAND"],
         CLUMPY = CLUMPY * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"]) %>% 
  gather(-1,-2,-3,-4,-5, key = "Variable", value = "Value") %>% group_by(Variable, country) %>% 
  summarise(quant_std = list(enframe(quantile(Value, probs=c(0.1,0.25,0.5,0.75,0.9))))) %>% 
  unnest

percentiles <- dat %>% dplyr:::filter(Scale == 20000) %>% group_by(Site) %>% 
  summarize_all(first) %>% dplyr::select(1,6,7,13,15,16,17) %>% 
  mutate(PLAND = PLAND * scaleList$scale["PLAND"] + scaleList$center["PLAND"],
         CLUMPY = CLUMPY * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"]) %>% 
  gather(-1,-2,-3,-4,-5, key = "Variable", value = "Value") %>% group_by(Variable) %>% 
  summarise(quant_std = list(enframe(quantile(Value, probs=c(0.1,0.25,0.5,0.75,0.9))))) %>% 
  unnest


CLUMPY_brks <- (c(0.80,0.91) - scaleList$center["CLUMPY"])/ scaleList$scale["CLUMPY"]
PLAND_brks <- (c(0.06,0.54) - scaleList$center["PLAND"])/ scaleList$scale["PLAND"]


emtrends(mod20, ~ CLUMPY | PLAND, var = "Year", 
         at = list(CLUMPY = (percentiles[percentiles$Variable == "CLUMPY", "value"][[1]] -
                               scaleList$center["CLUMPY"])/ scaleList$scale["CLUMPY"],
                   PLAND = (percentiles[percentiles$Variable == "PLAND", "value"][[1]] -
                              scaleList$center["PLAND"])/ scaleList$scale["PLAND"]))

emtrends(mod20, ~ CLUMPY | PLAND, var = "Year", at = list(CLUMPY = CLUMPY_brks,
                                                          PLAND = PLAND_brks))

# pred.test <- predict(mod20, newdata = data.frame(CLUMPY = (0.881 - scaleList$center["CLUMPY"])/scaleList$scale["CLUMPY"],
#                                                  PLAND = (0.253 - scaleList$center["PLAND"])/scaleList$scale["PLAND"],
#                                                  Year = 1995:2016,
#                                                  Habitat = "Open",
#                                                  X = median(dat$X),
#                                                  Y = median(dat$Y),
#                                                  Site = "398_NL",
#                                                  gridCell50 = "1126_NL"), re.form = NA)
# 
# lm(pred.test ~ c(1995:2016))

cti_spat_trend <- lmer(cti ~ Y + (Y|Year), data = dat %>% filter(Scale == 20000))

temp_temp_trend <- lmer(temp ~ Year + (1|Site), data = dat %>% filter(Scale == 20000))
temp_spat_trend <- lmer(temp ~ Y + (1|Year), data = dat %>% filter(Scale == 20000))

m <- lmer(temp ~ Year + X*Y + (1|Site), data = temp_summer_2)



# median frag:
0.002064154   / (-fixef(cti_spat_trend)[2] * 1000)

# low frag:
0.008174396 / (-fixef(cti_spat_trend)[2] * 1000)

# climate across all sites 1991-2016:
0.027318398 / (-(-2.953940e-06) * 1000)

visreg(mod20, "Year", by = "PLAND", breaks = PLAND_brks, rug = F, scale = "response", overlay = T, band = F, 
       cond = list(CLUMPY = CLUMPY_brks[1]), type = "contrast",
       xtrans = function(x){x + 1990},
       ylim = c(-.1, .1), axes = F, 
       line.par = list(lty = 2, col = c("#104E8B", "#CD2626")), ylab = "Community temperature index",
       legend = F)
par(new=T)
visreg(mod20, "Year", by = "PLAND", breaks = PLAND_brks, rug = F, scale = "response", overlay = T, band = F, 
       cond = list(CLUMPY = CLUMPY_brks[2]), type = "contrast",
       xtrans = function(x){x + 1990},
       ylim = c(-.1, .1), ann = F, axes = F, line.par = list(lty = 1, col = c("#104E8B", "#CD2626")), legend = F)
axis(1)
axis(2)



###########################
## remove unused objects ##
###########################


gdata:::keep(butterflies.data.presence, butterflies.data.abundance,
             mods, mods.out, 
             pred,
             sure = T)


