library(data.table)
library(lme4)
library(lmerTest)
library(MuMIn)
library(visreg)
library(cowplot)
library(nloptr)
library(rotl)
library(caper)
library(sticky)
source("functions.R")

###############################
########## Load data ########## 
###############################
#### prepare butterfly dataset ####
# butterflies.FIN <- as.tbl(fread("../Data/pres_abs_FIN_data.csv"))
# butterflies.NL <- as.tbl(fread("../Data/pres_abs_NL_data.csv"))
# 
# butterflies <- bind_rows(butterflies.FIN %>% mutate(country = "FIN"),
#                          butterflies.NL %>% mutate(country = "NL"))
# 
# traits.butterflies <- as.tbl(read_csv("../Data/Butterflytraits_Bink.csv"))
# traits.butterflies <- traits.butterflies %>% dplyr:::select(1,4,20,21)
# 
# traits.butterflies <- as.tbl(read_csv("../Data/Butterflies - Netherlands/SpeciesTraits_WDV2014.csv"))
# 
# sp.coef %>% dplyr:::filter(!Species %in% traits.butterflies$Scientific_name)
# 
# traits.butterflies[grep("-album", traits.butterflies$Scientific_name), "Scientific_name"]
# butterflies[grep("-album", butterflies$Species), "Species"]
# 
# 
# butterflies[grep("Polygonia c-album", butterflies$Species), "Species"] <- "Nymphalis c-album"
# 
# traits.butterflies[grep("Colias croceus", traits.butterflies$Scientific_name), "Scientific_name"] <- "Colias crocea"
# traits.butterflies[grep("walbum", traits.butterflies$Scientific_name), "Scientific_name"] <- "Satyrium w-album"
# traits.butterflies[grep("Neozephyrus quercus", traits.butterflies$Scientific_name), "Scientific_name"] <- "Favonius quercus"
# traits.butterflies[grep("lycaon", traits.butterflies$Scientific_name), "Scientific_name"] <- "Hyponephele lycaon"
# traits.butterflies[grep("Polygonia", traits.butterflies$Scientific_name), "Scientific_name"] <- "Nymphalis c-album"
# traits.butterflies[grep("Inachis io", traits.butterflies$Scientific_name), "Scientific_name"] <- "Aglais io"
# traits.butterflies[grep("tithonus", traits.butterflies$Scientific_name), "Scientific_name"] <- "Pyronia tithonus"
# traits.butterflies[grep("alcon", traits.butterflies$Scientific_name), "Scientific_name"] <- "Phengaris alcon"
# 
# butterflies <- as.tbl(merge(butterflies, 
#                             traits.butterflies %>% dplyr:::select(-Species, -STI), 
#                             by.x = "Species", by.y = "Scientific_name"))
# 
# write_csv(butterflies, "C:/Local Folder (c)/temp_butterflies_pres-abs_data.csv")

# load already prepared data

butterflies <- bind_rows(read_csv("../Data/pres_abs_FIN_data.csv"), read_csv("../Data/pres_abs_NL_data.csv"))

nbDat <- butterflies %>% dplyr:::filter(Scale == 3000) %>% group_by(country, Species) %>% summarise(tot = n())
# 
# # load cti data
# butterflies.data <- as.tbl(fread("../Data/cti_butterflies_data.csv"))
# butterflies.data <- butterflies.data %>% dplyr:::filter(type == "Presence")
# butterflies.data <- 
#   left_join(
#     butterflies.data %>%
#       group_by(coords = paste(X,Y)) %>% summarise(Site = first(Site)) %>%
#       ungroup(),
#     butterflies.data %>% mutate(coords = paste(X,Y)) %>% dplyr:::select(-2),
#     by = c("coords" = "coords")
#   ) %>% dplyr:::select(-1) %>%
#   group_by(Site, Year, Scale) %>% summarise_all(function(x){ifelse(is.numeric(x),mean(x, na.rm = TRUE), x)}) %>%
#   ungroup()
# 
# butterflies.data <- butterflies.data %>% dplyr:::filter(PLAND < 90) %>%
#   dplyr:::filter(!Site %in% subset(., is.na(CLUMPY))$Site)
# 
# butterflies.data <- stdize(butterflies.data, prefix = F, omit.cols= c("gridCell50", "Scale"))
# scaleList <- list(scale = attr(butterflies.data, "scaled:scale")[c("cti", "Year", "CLUMPY", "PLAND")],
#                   center = attr(butterflies.data, "scaled:center")[c("cti", "Year", "CLUMPY", "PLAND")])

# load phylogenetic data
taxa <- unique(butterflies$Species)
resolved_names <- tnrs_match_names(na.omit(taxa))

inspect(resolved_names, taxon_name = "aporia crataegi")
resolved_names <- update(resolved_names, taxon_name = "aporia crataegi",
                         new_row_number = 2)

inspect(resolved_names, taxon_name = "nymphalis c-album")
resolved_names <- update(resolved_names, taxon_name = "nymphalis c-album",
                         new_row_number = 3)

my_tree <- tol_induced_subtree(ott_ids = resolved_names$ott_id)
my_tree <- ape::compute.brlen(my_tree)

my_tree$tip.label <- gsub("_ott[0-9]*","",my_tree$tip.label)
my_tree$tip.label <- gsub("_"," ",my_tree$tip.label)
resolved_names$unique_name <- gsub(" \\(species in Holozoa\\)", "", resolved_names$unique_name)

my_tree$tip.label <- merge(data.frame(species = my_tree$tip.label), 
                           resolved_names, 
                           by.x = "species", by.y = "unique_name", sort=F)[,2]

my_tree$tip.label[my_tree$tip.label == "nymphalis c\\-album"] <- "Nymphalis c-album"


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

my_tree$tip.label <- firstup(my_tree$tip.label)
my_tree$node.label <- "NULL"

# plot(my_tree, no.margin=TRUE)

##################
## prepare data ##
##################

## select only sites without NA
butterflies <- butterflies  %>% dplyr:::filter(PLAND < .99 & PLAND > .01) %>%
  dplyr:::filter(!Site %in% subset(., is.na(CLUMPY))$Site)

## standardize ##
butterflies <- stdize(butterflies, prefix = F, omit.cols= c("gridCell50", "Scale", "n"))
scaleList2 <- list(scale = attr(butterflies, "scaled:scale")[c("Year", "CLUMPY", "PLAND")],
                   center = attr(butterflies, "scaled:center")[c("Year", "CLUMPY", "PLAND")])

## select scale and country ##
data.sel.3000.FIN <- butterflies %>% dplyr:::filter(country == "FIN", Scale == "3km")
data.sel.30000.FIN <- butterflies %>% dplyr:::filter(country == "FIN", Scale == "30km")
data.sel.3000.NL <- butterflies %>% dplyr:::filter(country == "NL", Scale == "3km")
data.sel.30000.NL <- butterflies %>% dplyr:::filter(country == "NL", Scale == "30km")

##########################
## extract coefficients ##
##########################

#########
## FIN ##
#########

## 3000 ##

coefs_fin.3000.Y <- sp_occupancy_trend_yearOnly(data.sel.3000.FIN)
coefs_fin.3000.Y$term <- gsub("Year", "Trend", coefs_fin.3000.Y$term)

coefs_fin.3000 <- sp_occupancy_trend_landscape(data.sel.3000.FIN)
coefs_fin.3000$term <- gsub("Year:PLAND", "Effect of %SNH", coefs_fin.3000$term)
coefs_fin.3000$term <- gsub("Year:CLUMPY", "Effect of clumping", coefs_fin.3000$term)

coefs_fin.3000 <- bind_rows(coefs_fin.3000.Y, coefs_fin.3000)

coefs_fin.3000$term <- as.factor(coefs_fin.3000$term); coefs_fin.3000$term <- relevel(coefs_fin.3000$term, "Trend")

coefs_fin.3000 <- coefs_fin.3000 %>% group_by(term) %>% arrange(estimate) %>%
  mutate(Species = factor(Species, Species))

ggplot(coefs_fin.3000, aes(x = Species, y = estimate)) + 
  geom_errorbar(aes(x = Species, ymin = estimate - std.error, ymax = estimate + std.error), width = 0, color = "grey") + 
  geom_point() + 
  facet_grid(~term, scale = "free") + 
  coord_flip() + 
  geom_hline(yintercept = 0)


## 30000 ##
coefs_fin.30000 <- sp_occupancy_trend_landscape(data.sel.30000.FIN)
coefs_fin.30000$term <- gsub("Year:PLAND", "Effect of %SNH", coefs_fin.30000$term)
coefs_fin.30000$term <- gsub("Year:CLUMPY", "Effect of clumping", coefs_fin.30000$term)

#########
## NL ##
#########

## 3000 ##

coefs_NL.3000.Y <- sp_occupancy_trend_yearOnly(data.sel.3000.NL)
coefs_NL.3000.Y$term <- gsub("Year", "Trend", coefs_NL.3000.Y$term)

coefs_NL.3000 <- sp_occupancy_trend_landscape(data.sel.3000.NL)
coefs_NL.3000$term <- gsub("Year:PLAND", "Effect of %SNH", coefs_NL.3000$term)
coefs_NL.3000$term <- gsub("Year:CLUMPY", "Effect of clumping", coefs_NL.3000$term)

coefs_NL.3000 <- bind_rows(coefs_NL.3000.Y, coefs_NL.3000)

coefs_NL.3000$term <- as.factor(coefs_NL.3000$term); coefs_NL.3000$term <- relevel(coefs_NL.3000$term, "Trend")

coefs_NL.3000 <- coefs_NL.3000 %>% group_by(term) %>% arrange(estimate) %>%
  mutate(Species = factor(Species, Species))

ggplot(coefs_NL.3000, aes(x = Species, y = estimate)) + 
  geom_errorbar(aes(x = Species, ymin = estimate - std.error, ymax = estimate + std.error), width = 0, color = "grey") + 
  geom_point() + 
  facet_grid(~term, scale = "free") + 
  coord_flip() + 
  geom_hline(yintercept = 0)


## 30000 ##
coefs_NL.30000 <- sp_occupancy_trend_landscape(data.sel.30000.NL)
coefs_NL.30000$term <- gsub("Year:PLAND", "Effect of %SNH", coefs_NL.30000$term)
coefs_NL.30000$term <- gsub("Year:CLUMPY", "Effect of clumping", coefs_NL.30000$term)


##############
## analyses ##
##############

coefs_fin.3000 %>% ungroup() %>% group_by(term, direction = ifelse(estimate < 0 & p.value < .05, "Negative", 
                                                                   ifelse(estimate > 0 & p.value < .05, "Positive", "Nonsignificant"))) %>% 
  summarise(n = n())
coefs_NL.3000 %>% ungroup() %>% group_by(term, direction = ifelse(estimate < 0 & p.value < .05, "Negative", 
                                                                   ifelse(estimate > 0 & p.value < .05, "Positive", "Nonsignificant"))) %>% 
  summarise(n = n())

coefs_fin.30000 %>% ungroup() %>% group_by(term, direction = ifelse(estimate < 0 & p.value < .05, "Negative", 
                                                                   ifelse(estimate > 0 & p.value < .05, "Positive", "Nonsignificant"))) %>% 
  summarise(n = n())
coefs_NL.30000 %>% ungroup() %>% group_by(term, direction = ifelse(estimate < 0 & p.value < .05, "Negative", 
                                                                  ifelse(estimate > 0 & p.value < .05, "Positive", "Nonsignificant"))) %>% 
  summarise(n = n())


plot(estimate.x ~ estimate.y, data = inner_join(coefs_fin.3000 %>% ungroup() %>% dplyr::filter(term == "Trend") %>% 
                                                  dplyr::select(Species, estimate), 
                                                coefs_NL.3000 %>% ungroup() %>% dplyr::filter(term == "Trend") %>% 
                                                  dplyr::select(Species, estimate), by = "Species"))

plot(estimate.x ~ estimate.y, data = inner_join(coefs_fin.3000 %>% ungroup() %>% dplyr::filter(term == "Effect of %SNH") %>% 
                                                  dplyr::select(Species, estimate), 
                                                coefs_NL.3000 %>% ungroup() %>% dplyr::filter(term == "Effect of %SNH") %>%
                                                  dplyr::select(Species, estimate), by = "Species"))

plot(estimate.x ~ estimate.y, data = inner_join(coefs_fin.3000 %>% ungroup() %>% dplyr::filter(term == "Effect of clumping") %>%
                                                  dplyr::select(Species, estimate), 
                                                coefs_NL.3000 %>% ungroup() %>% dplyr::filter(term == "Effect of clumping") %>%
                                                  dplyr::select(Species, estimate), by = "Species"))

#########
## FIN ##
#########

coefs_fin.3000.Trend <- as.data.frame(coefs_fin.3000 %>% dplyr:::filter(term == "Trend")) %>% column_to_rownames("Species") %>% 
  do(na.omit(.)) %>% do(transPower(.))

m.traits.fin.3000.Trend <- lm(estimate ~ (STI_sd + PC1 + Habitat + PC3)*STI, weights = 1/std.error,
                                data = coefs_fin.3000.Trend, na.action = na.fail)
car::qqPlot(resid(m.traits.fin.3000.Trend))
summary(m.traits.fin.3000.Trend)
summary(stats::step(m.traits.fin.3000.Trend))

m.traits.fin.3000.Trend.avg <- model.avg(dredge(m.traits.fin.3000.Trend), subset = cumsum(weight) <= 95)
summary(m.traits.fin.3000.Trend.avg)

# 3000 #
coefs_fin.3000.SNH <- as.data.frame(coefs_fin.3000 %>% dplyr:::filter(term == "Effect of %SNH")) %>% column_to_rownames("Species")  %>% 
  do(na.omit(.)) %>% do(transPower(.))

m.traits.fin.3000.SNH <- lm(estimate ~ (STI_sd + PC1 + Habitat + PC3)*STI, weights = 1/std.error,
                              data = coefs_fin.3000.SNH, na.action = na.fail)
car::qqPlot(resid(m.traits.fin.3000.SNH))
summary(m.traits.fin.3000.SNH)
summary(stats::step(m.traits.fin.3000.SNH))

m.traits.fin.3000.SNH.avg <- model.avg(dredge(m.traits.fin.3000.SNH), subset = cumsum(weight) <= 95)
summary(m.traits.fin.3000.SNH.avg)


coefs_fin.3000.Clumping <- as.data.frame(coefs_fin.3000 %>% dplyr:::filter(term == "Effect of clumping")) %>% column_to_rownames("Species") %>% 
  do(na.omit(.)) %>% do(transPower(.))

m.traits.fin.3000.Clumping <- lm(estimate ~ (STI_sd + PC1 + Habitat + PC3)*STI, weights = 1/std.error,
                                 data = coefs_fin.3000.Clumping, na.action = na.fail)

car::qqPlot(resid(m.traits.fin.3000.Clumping))
summary(m.traits.fin.3000.Clumping)
summary(stats::step(m.traits.fin.3000.Clumping))

m.traits.fin.3000.Clumping.avg <- model.avg(dredge(m.traits.fin.3000.Clumping), subset = cumsum(weight) <= 95)
summary(m.traits.fin.3000.Clumping.avg)

# 30000 #
coefs_fin.30000.SNH <- as.data.frame(coefs_fin.30000 %>% dplyr:::filter(term == "Effect of %SNH")) %>% column_to_rownames("Species")  %>% 
  do(na.omit(.)) %>% do(transPower(.))

m.traits.fin.30000.SNH <- lm(estimate ~ (STI_sd + PC1 + Habitat + PC3)*STI, weights = 1/std.error, 
                              data = coefs_fin.30000.SNH, na.action = na.fail)
car::qqPlot(resid(m.traits.fin.30000.SNH))
summary(m.traits.fin.30000.SNH)
summary(stats::step(m.traits.fin.30000.SNH))

m.traits.fin.30000.SNH.avg <- model.avg(dredge(m.traits.fin.30000.SNH), subset = cumsum(weight) <= 95)
summary(m.traits.fin.30000.SNH.avg)


coefs_fin.30000.Clumping <- as.data.frame(coefs_fin.30000 %>% dplyr:::filter(term == "Effect of clumping")) %>% column_to_rownames("Species")  %>% 
  do(na.omit(.)) %>% do(transPower(.))

m.traits.fin.30000.Clumping <- lm(estimate ~ (STI_sd + PC1 + Habitat + PC3)*STI, weights = 1/std.error,
                                   data = coefs_fin.30000.Clumping, na.action = na.fail)
car::qqPlot(resid(m.traits.fin.30000.Clumping))
summary(m.traits.fin.30000.Clumping)
summary(stats::step(m.traits.fin.30000.Clumping))

m.traits.fin.30000.Clumping.avg <- model.avg(dredge(m.traits.fin.30000.Clumping), subset = cumsum(weight) <= 95)
summary(m.traits.fin.30000.Clumping.avg)


########
## NL ##
########
coefs_NL.3000.Trend <- as.data.frame(coefs_NL.3000 %>% dplyr:::filter(term == "Trend")) %>% column_to_rownames("Species") %>% 
  do(na.omit(.)) %>% do(transPower(.))

m.traits.NL.3000.Trend <- lm(estimate ~ (STI_sd + PC1 + Habitat + PC3)*STI, weights = 1/std.error,
                                data = coefs_NL.3000.Trend, na.action = na.fail)
car::qqPlot(resid(m.traits.NL.3000.Trend))
summary(m.traits.NL.3000.Trend)
summary(stats::step(m.traits.NL.3000.Trend))

m.traits.NL.3000.Trend.avg <- model.avg(dredge(m.traits.NL.3000.Trend), subset = cumsum(weight) <= 95)
summary(m.traits.NL.3000.Trend.avg)

# 3000 #
coefs_NL.3000.SNH <- as.data.frame(coefs_NL.3000 %>% dplyr:::filter(term == "Effect of %SNH")) %>% column_to_rownames("Species") %>% 
  do(na.omit(.)) %>% do(transPower(.))

m.traits.NL.3000.SNH <- lm(estimate ~ (STI_sd + PC1 + Habitat + PC3)*STI, weights = 1/std.error,
                              data = coefs_NL.3000.SNH, na.action = na.fail)
car::qqPlot(resid(m.traits.NL.3000.SNH))
summary(m.traits.NL.3000.SNH)
summary(stats::step(m.traits.NL.3000.SNH))

m.traits.NL.3000.SNH.avg <- model.avg(dredge(m.traits.NL.3000.SNH), subset = cumsum(weight) <= 95)
summary(m.traits.NL.3000.SNH.avg)


coefs_NL.3000.Clumping <- as.data.frame(coefs_NL.3000 %>% dplyr:::filter(term == "Effect of clumping")) %>% column_to_rownames("Species") %>% 
  do(na.omit(.)) %>% do(transPower(.))

m.traits.NL.3000.Clumping <- lm(estimate ~ (STI_sd + PC1 + Habitat + PC3)*STI, weights = 1/std.error,
                                data = coefs_NL.3000.Clumping, na.action = na.fail)
car::qqPlot(resid(m.traits.NL.3000.Clumping))
summary(m.traits.NL.3000.Clumping)
summary(stats::step(m.traits.NL.3000.Clumping))

m.traits.NL.3000.Clumping.avg <- model.avg(dredge(m.traits.NL.3000.Clumping), subset = cumsum(weight) <= 95)
summary(m.traits.NL.3000.Clumping.avg)

# 30000 #
coefs_NL.30000.SNH <- as.data.frame(coefs_NL.30000 %>% dplyr:::filter(term == "Effect of %SNH")) %>% column_to_rownames("Species") %>% 
  do(na.omit(.)) %>% do(transPower(.))

m.traits.NL.30000.SNH <- lm(estimate ~ (STI_sd + PC1 + Habitat + PC3)*STI, weights = 1/std.error, 
                             data = coefs_NL.30000.SNH, na.action = na.fail)
car::qqPlot(resid(m.traits.NL.30000.SNH))
summary(m.traits.NL.30000.SNH)
summary(stats::step(m.traits.NL.30000.SNH))

m.traits.NL.30000.SNH.avg <- model.avg(dredge(m.traits.NL.30000.SNH), subset = cumsum(weight) <= 95)
summary(m.traits.NL.30000.SNH.avg)


coefs_NL.30000.Clumping <- as.data.frame(coefs_NL.30000 %>% dplyr:::filter(term == "Effect of clumping")) %>% column_to_rownames("Species") %>% 
  do(na.omit(.)) %>% do(transPower(.))

m.traits.NL.30000.Clumping <- lm(estimate ~ (STI_sd + PC1 + Habitat + PC3)*STI, weights = 1/std.error, 
                                  data = coefs_NL.30000.Clumping, na.action = na.fail)
car::qqPlot(resid(m.traits.NL.30000.Clumping))
summary(m.traits.NL.30000.Clumping)
summary(stats::step(m.traits.NL.30000.Clumping))

m.traits.NL.30000.Clumping.avg <- model.avg(dredge(m.traits.NL.30000.Clumping), subset = cumsum(weight) <= 95)
summary(m.traits.NL.30000.Clumping.avg)


###########
## plots ##
###########

coefs.m.traits <- rbind.data.frame(cbind.data.frame(country = "Finland", Scale = "None", Variable = "Occupancy Trend", 
                                                    importance = as.vector(importance(m.traits.fin.3000.Trend.avg)),
                                                    estimate = coef(m.traits.fin.3000.Trend.avg, full =T)[-1], 
                                                    confint(m.traits.fin.3000.Trend.avg, full =T)[-1,]) %>% rownames_to_column("Traits"),
                                   cbind.data.frame(country = "Finland", Scale = "3 Km", Variable = "%SNH", 
                                                    importance = as.vector(importance(m.traits.fin.3000.SNH.avg)),
                                                    estimate = coef(m.traits.fin.3000.SNH.avg, full =T)[-1], 
                                                    confint(m.traits.fin.3000.SNH.avg, full =T)[-1,]) %>% rownames_to_column("Traits"),
                                   cbind.data.frame(country = "Finland", Scale = "3 Km", Variable = "Clumping",
                                                    importance = as.vector(importance(m.traits.fin.3000.Clumping.avg)),
                                                    estimate = coef(m.traits.fin.3000.Clumping.avg, full =T)[-1], 
                                                    confint(m.traits.fin.3000.Clumping.avg, full =T)[-1,]) %>% rownames_to_column("Traits"),
                                   cbind.data.frame(country = "Finland", Scale = "30 Km", Variable = "%SNH", 
                                                    importance = as.vector(importance(m.traits.fin.30000.SNH.avg)),
                                                    estimate = coef(m.traits.fin.30000.SNH.avg, full =T)[-1], 
                                                    confint(m.traits.fin.30000.SNH.avg, full =T)[-1,]) %>% rownames_to_column("Traits"),
                                   cbind.data.frame(country = "Finland", Scale = "30 Km", Variable = "Clumping",
                                                    importance = as.vector(importance(m.traits.fin.30000.Clumping.avg)),
                                                    estimate = coef(m.traits.fin.30000.Clumping.avg, full =T)[-1], 
                                                    confint(m.traits.fin.30000.Clumping.avg, full =T)[-1,]) %>% rownames_to_column("Traits"),
                                   cbind.data.frame(country = "Netherlands", Scale = "None", Variable = "Occupancy Trend", 
                                                    importance = as.vector(importance(m.traits.NL.3000.Trend.avg)),
                                                    estimate = coef(m.traits.NL.3000.Trend.avg, full =T)[-1], 
                                                    confint(m.traits.NL.3000.Trend.avg, full =T)[-1,]) %>% rownames_to_column("Traits"),
                                   cbind.data.frame(country = "Netherlands", Scale = "3 Km", Variable = "%SNH", 
                                                    importance = as.vector(importance(m.traits.NL.3000.SNH.avg)),
                                                    estimate = coef(m.traits.NL.3000.SNH.avg, full =T)[-1], 
                                                    confint(m.traits.NL.3000.SNH.avg, full =T)[-1,]) %>% rownames_to_column("Traits"),
                                   cbind.data.frame(country = "Netherlands", Scale = "3 Km", Variable = "Clumping",
                                                    importance = as.vector(importance(m.traits.NL.3000.Clumping.avg)),
                                                    estimate = coef(m.traits.NL.3000.Clumping.avg, full =T)[-1], 
                                                    confint(m.traits.NL.3000.Clumping.avg, full =T)[-1,]) %>% rownames_to_column("Traits"),
                                   cbind.data.frame(country = "Netherlands", Scale = "30 Km", Variable = "%SNH", 
                                                    importance = as.vector(importance(m.traits.NL.30000.SNH.avg)),
                                                    estimate = coef(m.traits.NL.30000.SNH.avg, full =T)[-1], 
                                                    confint(m.traits.NL.30000.SNH.avg, full =T)[-1,]) %>% rownames_to_column("Traits"),
                                   cbind.data.frame(country = "Netherlands", Scale = "30 Km", Variable = "Clumping",
                                                    importance = as.vector(importance(m.traits.NL.30000.Clumping.avg)),
                                                    estimate = coef(m.traits.NL.30000.Clumping.avg, full =T)[-1], 
                                                    confint(m.traits.NL.30000.Clumping.avg, full =T)[-1,]) %>% rownames_to_column("Traits"))

coefs.m.traits <- coefs.m.traits %>% complete(Traits, nesting(country, Scale, Variable))
coefs.m.traits <- coefs.m.traits %>% mutate(importance = ifelse(is.na(importance), 0, importance))
coefs.m.traits$Traits <- factor(coefs.m.traits$Traits, 
                                levels = c("STI", "STI_sd", "PC1", "PC3",
                                           "STI:STI_sd", "PC1:STI", "PC3:STI"))
coefs.m.traits$Sig <- ifelse(coefs.m.traits$`97.5 %` * coefs.m.traits$`2.5 %` > 0 , "significant", "nonsignificant")
coefs.m.traits$Sig[is.na(coefs.m.traits$Sig)] <- "nonsignificant"


ggplot(coefs.m.traits %>% dplyr::filter(Variable != "Occupancy trend"), aes(x = Traits, y = estimate)) + 
  geom_bar(aes(x = Traits, y = importance/3, fill = Scale), position = position_dodge(.7), width = .7, stat = 'identity') +
  geom_errorbar(aes(ymin = `2.5 %`, ymax = `97.5 %`, color = Scale, linetype = Sig), width = .2, position = position_dodge(.7)) + 
  geom_point(aes(color = Scale, shape = Sig, size = Sig), position = position_dodge(.7)) + 
  geom_hline(yintercept = 0) +
  facet_grid(country ~ Variable) +
  scale_fill_manual(values = c("grey", "#1C86F230", "#F0656560")) + 
  scale_color_manual(values = c("grey3","#1C86EE", "#EE6363")) +
  scale_y_continuous("Coefficient estimates\nvariable importance / 3") + 
  scale_shape_manual(values = c(16, 15), guide = "none") + 
  scale_size_manual(values = c(1.8, 2.5), guide = "none") + 
  scale_linetype_manual(values = c(2,1), guide = "none") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("../plot.png", width = 10, height = 7.5)


###########################
## remove unused objects ##
###########################


gdata:::keep(coefs_NL.30000, coefs_NL.3000, coefs_fin.30000, coefs_fin.3000, my_tree, 
             sure = T)
