library(MuMIn)
library(lmerTest)
library(data.table)
library(broom)
library(cowplot)
library(doFuture)
source("functions.R")

## load data ##
butterflies.cti.presence <- as.tbl(fread("../Data/cti_butterflies_data.csv")) %>% dplyr::filter(type == "Presence")
butterflies <- bind_rows(read_csv("../Data/pres_abs_FIN_data.csv"), read_csv("../Data/pres_abs_NL_data.csv"))

butterflies <- left_join(butterflies %>% group_by(coords = paste(X,Y)),
                         butterflies.cti.presence %>% group_by(coords = paste(X,Y)) %>% summarise(cti = mean(cti)),
          by = "coords")

butterflies <- butterflies %>% mutate(STI_rel = cti - STI) %>%
  dplyr::filter(!is.na(STI)) %>%
  mutate(gridCell50 = ifelse(country == "NL", paste0(gridCell50, "_NL"), paste0(gridCell50, "_FIN")))
butterflies <- butterflies %>% as.data.frame %>% polypoly::poly_add_columns(.col = PC1, degree = 2) %>% as.tbl

write_csv(butterflies, "C:/Local Folder (c)/butterflies_occ.csv")

butterflies <- butterflies %>% stdize(., prefix = F, omit.cols = c("n", "gridCell50", "Scale"))
scaleList <- list(scale = attr(butterflies, "scaled:scale")[c("Year", "CLUMPY", "PLAND")],
                  center = attr(butterflies, "scaled:center")[c("Year", "CLUMPY", "PLAND")])

write_csv(butterflies, "C:/Local Folder (c)/butterflies_occ_std.csv")

butterflies <- as.tbl(fread("C:/Local Folder (c)/butterflies_occ_std.csv"))
scaleList <- structure(list(scale = structure(c(7.76768830494924, 0.054126736142142,0.213524890144323),
                                              .Names = c("Year", "CLUMPY", "PLAND")),
                            center = structure(c(2003.20121386844,0.867284438162297, 0.396710732822166),
                                               .Names = c("Year", "CLUMPY","PLAND"))),
                       .Names = c("scale", "center"))

##############################
#### models at all scales ####
##############################

registerDoFuture()
options(future.globals.maxSize = Inf)
plan(multiprocess, workers = 2)

foreach(i = unique(butterflies$Scale), .combine = rbind.data.frame) %dopar% {
  # model with STI and quadratic dispersal
  m <- glmer(n ~ Year * PLAND * CLUMPY * STI_rel * PC11 +
               Year * PLAND * CLUMPY * STI_rel * PC12 +
               X * Y +
               Year * Habitat +
               (1 | gridCell50/Site) + (Year | Species),
             family = binomial,
             data = butterflies %>% dplyr::filter(Scale == i) %>% 
               mutate(n = ifelse(n < .5, 0, 1)),
             nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                           optCtrl=list(maxfun=1e10),
                                           calc.derivs = FALSE), verbose = F)
  saveRDS(m, file = paste0("C:/Local Folder (c)/Models occurrence trend/model_occ_", i/1000, "km.RDS"))
  rm(m)
  gc()
}


######################
#### plot results ####
######################

mods.out <- list.files("C:/Local Folder (c)/Models occurrence trend/", full.names = T)

pred <- foreach(i = 1:length(mods.out), .combine = rbind) %do% {
  m <- readRDS(mods.out[[i]])
  pred_temp <- predict_raster2(m, scaleList, 10)
  pred_temp <- cbind.data.frame(Scale = paste(regmatches(mods.out,regexec("/model_occ_(.*?)km.RDS",mods.out))[[i]][2], "km"),
                                pred_temp)
  rm(m)
  gc()
  return(pred_temp)
}

# levels(pred$Scale) <- factor(pred$Scale, levels =  paste(c(1,3,5,10,20,30,50), "km"))


ggplot(pred, aes(x= x, y = y, fill = estimate)) + 
  geom_raster() + 
  facet_grid(STI_rel + PC1 ~ Scale) +
  scale_fill_gradientn(name = "Occupancy trend\n", 
                       breaks =  c(min(pred$estimate),0,max(pred$estimate)),
                       labels = c(round(min(pred$estimate), 3),"0.000",round(max(pred$estimate), 3)),
                       limits = c(min(pred$estimate),max(pred$estimate)),
                       colours=c("light blue", "yellow","red")) +
  scale_x_continuous("% Semi-natural habitat") + scale_y_continuous("Habitat clumpiness")



###########################
## remove unused objects ##
###########################


gdata:::keep(m_NL_3000, m_NL_30000, m_FIN_3000, m_FIN_30000, 
             sure = T)
