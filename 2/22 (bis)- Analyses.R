library(tidyverse)
library(doFuture)
library(future.apply)
library(lmerTest)
library(visreg)
library(abind)
library(ape)
library(cowplot)

#### load
data <- read_csv("../SWE_PA_data.csv") %>% filter(latin != "Branta leucopsis")
traits <- read_csv("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Birds trait database/traits_birds_corrected.csv")
phyl_dist <- read.csv("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Birds phylogeny/phylo_dist.csv")
sti <- read_csv("../Data/Birds Atlas/EBCC1/EBCC1_sti.csv")

#### calculate turnover per site based on two time periods
turnover.per.site <- data %>% filter(abundance > 0) %>% group_by(karta) %>% 
  mutate(yr = yr - min(yr), time = as.numeric(cut(yr, 2, labels = c(1,2)))) %>% 
  group_by(karta, time) %>% nest(latin, .key = "Species") %>% 
  mutate(Species = lapply(Species, unique)) %>% group_by(karta) %>% arrange(karta, time) %>%
  mutate(previous = lag(Species)) %>% filter(time == 2) %>% rowwise() %>%
  mutate(gained = list(setdiff(Species, previous)), 
         lost = list(setdiff(previous, Species)),
         n.gained = nrow(setdiff(Species, previous)), 
         n.lost = nrow(setdiff(previous, Species)), 
         turnover = n.lost + n.gained,
         prop.gained = n.gained/nrow(previous),
         prop.lost = n.lost/nrow(previous),
         prop.turnover = turnover / nrow(previous)) %>%
  left_join(data %>% group_by(karta) %>% summarise(x = unique(x), y = unique(y)), by = "karta") %>%
  rename(sp.first = previous, sp.last = Species)

# extract gained and lost species between 1st end 2nd time periods 
gained.sp <- table(unlist(turnover.per.site$gained)); gained.sp <- gained.sp[order(gained.sp, decreasing = T)]
lost.sp <- table(unlist(turnover.per.site$lost)); lost.sp <- lost.sp[order(lost.sp, decreasing = T)]

## exploratory plots
plot(prop.gained ~ n.first, data = turnover.per.site %>% mutate(n.first = nrow(sp.first)))
abline(lm(prop.gained ~ n.first, data = turnover.per.site %>% mutate(n.first = nrow(sp.first))))

plot(prop.lost ~ n.first, data = turnover.per.site %>% mutate(n.first = nrow(sp.first)))
abline(lm(prop.lost ~ n.first, data = turnover.per.site %>% mutate(n.first = nrow(sp.first))))

plot(I(n.lost/n.gained) ~ n.first, data = turnover.per.site %>% mutate(n.first = nrow(sp.first)))
abline(lm(I(n.lost/n.gained) ~ n.first, data = turnover.per.site %>% mutate(n.first = nrow(sp.first))))

plot(prop.turnover ~ n.first, data = turnover.per.site %>% mutate(n.first = nrow(sp.first)))
abline(lm(prop.turnover ~ n.first, data = turnover.per.site %>% mutate(n.first = nrow(sp.first), turnover = n.lost + n.gained)))

### compute association between colonisation and extinctions
## empirical
registerDoFuture()
plan(multiprocess)

Assoc <- foreach(i = unique(names(gained.sp)), .combine = rbind.data.frame) %dopar%{
  turnover.per.site[unlist(lapply(turnover.per.site$gained, 
                                  function(x){any(x$latin %in% i)})),] %>% 
    filter(n.lost >0) %>%
    rowwise() %>% do(l = cbind(Site = .$karta, Gained = i, .$lost)) %>% 
    unnest() %>% rename(Lost = latin)
}
Assoc <- Assoc %>% mutate(Gained = as.character(Gained))


Assoc <- merge(Assoc, data.frame(Gained = unique(do.call(rbind, turnover.per.site$sp.last))[[1]]),
               all = T, by = "Gained")
Assoc <- merge(Assoc, data.frame(Lost = unique(do.call(rbind, turnover.per.site$sp.first))[[1]]),
               all = T, by = "Lost")

# compute table of number of associations 
assoc.sum <- table(Assoc[c("Gained", "Lost")], useNA = "no")
assoc.sum <- assoc.sum[sort(rownames(assoc.sum)), sort(colnames(assoc.sum))]

## potential
ref_dat <- data %>% filter(abundance > 0) %>% group_by(karta) %>% 
  mutate(yr = yr - min(yr), time = as.numeric(cut(yr, 2, labels = c(1,2)))) %>% 
  filter(time == 1) %>% select(1,3) %>% unique() %>% rename(Lost = latin) %>%
  left_join(data %>% filter(abundance > 0) %>% group_by(karta) %>% 
              mutate(yr = yr - min(yr), time = as.numeric(cut(yr, 2, labels = c(1,2)))) %>% 
              filter(time == 2) %>% select(1,3) %>% unique()) %>% rename(Gained = latin) %>% filter(Gained != Lost)

ref_dat <- table(ref_dat[c("Gained", "Lost")], useNA = "no")
ref_dat <- ref_dat[sort(rownames(ref_dat)), sort(colnames(ref_dat))]

### compute no. of observed "interactions" / no. of potential interactions
perm.test <- as.data.frame(as.table(ref_dat)) %>% rename(Potential_lost = Freq) %>%
  left_join(as.data.frame(as.table(assoc.sum)) %>% rename(Observed_lost = Freq)) %>% 
  mutate_at(c(1,2),.funs = as.character) %>% filter(Gained != Lost) %>% as.tbl


# # permutation test
# 
# ab.per.site <- data %>% group_by(latin,karta) %>% summarise(n = mean(abundance, na.rm = T))
# 
# n.rep <- 1000
# 
# avg.inter <- future_replicate(n.rep, {
#   perm.dat <- turnover.per.site %>% select(1,4) %>% unnest %>% left_join(turnover.per.site %>% select(1,5)) %>% unnest %>%
#     rename(Site = karta, Lost = latin, Gained = latin1) %>%
#     filter(Lost %in% names(lost.sp)) %>% group_by(Site, Gained) %>%
#     sample_n(1) %>% ungroup
#   perm.dat <- perm.dat %>% mutate(Gained = as.character(Gained), Lost = as.character(Lost))
#   
#   perm.dat <- merge(perm.dat, data.frame(Gained = unique(do.call(rbind, turnover.per.site$sp.last))[[1]]), 
#                     all = T, by = "Gained")
#   perm.dat <- merge(perm.dat, data.frame(Lost = unique(do.call(rbind, turnover.per.site$sp.first))[[1]]), 
#                     all = T, by = "Lost")
#   
#   perm.dat.sum <- table(perm.dat[,c("Gained", "Lost")], useNA = "no")
# }
# )
# 
# # extract p-values for pairwise associations
# plan(multiprocess, workers = 3)
# 
# p.values <- future_apply(abind(avg.inter, assoc.sum), c(1,2), 
#                          function(x){(sum(x[1:n.rep] >= x[n.rep+1]) + 1) / (n.rep + 1)})
# 
# obs.quantile <- future_apply(abind(avg.inter, assoc.sum), c(1,2), 
#                              function(x){quantile(x[1:n.rep], probs = c(.05,.5,.95))})
# 
# 
# perm.test <- as.data.frame(as.table(p.values)) %>% rename(Gained = Var1, Lost = Var2, p.value = Freq) %>%
#   left_join(as.data.frame(as.table(obs.quantile)) %>% rename(Gained = Var2, Lost = Var3) %>%
#               spread(key = Var1, value = Freq)) %>% rename(Perm.median = `50%`, LowCI = `5%`, UprCI = `95%`) %>%
#   left_join(as.data.frame(as.table(assoc.sum)) %>% rename(Observed = Freq)) %>% 
#   as.tbl() %>% select(1,2,7,5,4,6,3) %>% mutate(adj.p.value = p.adjust(p.value, method = "fdr"))
# 
# 
# write_csv(perm.test, "../perm.test.csv")
# perm.test <- read_csv("../perm.test.csv")
# 
# perm.test %>% filter(adj.p.value < 0.05)
# perm.test %>% filter(p.value < .05)
# 
# 
# # example plot
# i = "Accipiter nisus"          
# j = "Accipiter gentilis"
# p.value = (sum(avg.inter[i,j,] >= assoc.sum[i,j]) + 1) / (length(avg.inter[i,j,]) + 1)
# hist(avg.inter[i,j,], xlim = c(0, max(assoc.sum[i,j], max(avg.inter[i,j,]))+1), main = p.value)
# abline(v=assoc.sum[i,j], col="red", lwd=2)
# 
# 
# perm.test %>% filter(p.value < .05) %>% group_by(Gained) %>% summarise(n = n()) %>% mutate(Gained = fct_reorder(Gained, n)) %>%
#   ggplot(aes(x = n, y = Gained)) + geom_point()
# 
# perm.test %>% filter(p.value < .05) %>% group_by(Lost) %>% summarise(n = n()) %>% mutate(Lost = fct_reorder(Lost, n)) %>%
#   ggplot(aes(x = n, y = Lost)) + geom_point()
# 
# perm.test %>% mutate(diff = ifelse(Observed - Perm.median < 0, 0, Observed - Perm.median )) %>% do(hist((.$diff)))

## add traits
funct_dist <- as.matrix(FD::gowdis(as.data.frame(traits[,c(5,6,9,12,15,18,21,24,25,28:85)])))
rownames(funct_dist) <- traits$Species; colnames(funct_dist) <- traits$Species
funct_dist <- subset(reshape2::melt(funct_dist), value!=0) %>% rename(Gained = Var1, Lost = Var2, funct_dist = value) %>% as.tbl
funct_dist <- bind_rows(funct_dist, funct_dist %>% rename(Gained = Lost, Lost = Gained))


perm.test.traits <- perm.test %>% filter(as.character(Gained) != as.character(Lost)) %>% 
  left_join(traits %>% select(2:4, 53:67), by = c("Gained" = "Species")) %>% 
  left_join(traits %>% select(2:4, 53:67), by = c("Lost" = "Species")) %>% 
  mutate(Deciduous.forest = Deciduous.forest.x + Deciduous.forest.y,
         Coniferous.forest = Coniferous.forest.x + Coniferous.forest.y,
         Woodland = Woodland.x + Woodland.y,
         Shrub = Shrub.x + Shrub.y,
         Savanna = Savanna.x + Savanna.y,
         Tundra = Tundra.x + Tundra.y,
         Grassland = Grassland.x + Grassland.y,
         Mountain.meadows = Mountain.meadows.x + Mountain.meadows.y,
         Reed = Reed.x + Reed.y,
         Swamps = Swamps.x + Swamps.y,
         Desert = Desert.x + Desert.y,
         Freshwater = Freshwater.x + Freshwater.y,
         Marine = Marine.x + Marine.y,
         Rocks = Rocks.x + Rocks.y,
         Human.settlements = Human.settlements.x + Human.settlements.y, 
         Same.order = ifelse(Order.x == Order.y, "Yes", "No"), Same.family = ifelse(Family.x == Family.y, "Yes", "No")) %>%
  mutate(Habitat_width.Gained = rowSums(.[c(7:21)]), Habitat_width.Lost = rowSums(.[c(24:38)]),
         Diff_habitat_width = Habitat_width.Gained - Habitat_width.Lost) %>% select(-c(5:38)) %>% 
  mutate_at(.funs = function(x){ifelse(x == 2, 1, 0)}, .vars = 5:19) %>% 
  mutate(Common_habitats = rowSums(.[c(5:19)])) %>% select(-c(5:19)) %>%
  left_join(traits, by = c("Gained" = "Species")) %>% 
  left_join(traits, by = c("Lost" = "Species")) %>%
  mutate(Folivore_B = Folivore_B.x + Folivore_B.y,
         Frugivore_B = Frugivore_B.x + Frugivore_B.y,
         Granivore_B = Granivore_B.x + Granivore_B.y,
         Arthropods_B = Arthropods_B.x + Arthropods_B.y,
         Other.invertebrates_B = Other.invertebrates_B.x + Other.invertebrates_B.y,
         Fish_B = Fish_B.x + Fish_B.y,
         Other.vertebrates_B = Other.vertebrates_B.x + Other.vertebrates_B.y,
         Carrion_B = Carrion_B.x + Carrion_B.y,
         Omnivore_B = Omnivore_B.x + Omnivore_B.y,
         Migrant.Gained = Long.distance.migrant.x + Short.distance.migrant.x,
         Migrant.Lost = Long.distance.migrant.y + Short.distance.migrant.y) %>%
  mutate(Diet_width.Gained = rowSums(.[c(86:94)]), Diet_width.Lost = rowSums(.[c(171:179)]),
         Diff_diet_width = Diet_width.Gained - Diet_width.Lost, 
         Diff_weight = WeightU_MEAN.x - WeightU_MEAN.y,
         weight.Gained = WeightU_MEAN.x, weight.Lost = WeightU_MEAN.y,
         diff_weight = WeightU_MEAN.x - WeightU_MEAN.y,
         repro.Gained = Clutch_MEAN.x * Broods.per.year.x, repro.Lost = Clutch_MEAN.y * Broods.per.year.y,
         diff_repro = repro.Gained - repro.Lost) %>% 
  mutate_at(.funs = function(x){ifelse(x == 2, 1, 0)}, .vars = 181:189) %>% 
  mutate(Common_diet = rowSums(.[c(181:189)])) %>% select(-c(11:189)) %>%
  mutate(Migrant.Gained = ifelse(Migrant.Gained > 0, "Migrant", "Resident"), 
         Migrant.Lost = ifelse(Migrant.Lost > 0, "Migrant", "Resident")) %>%
  left_join(funct_dist) %>% unique()

perm.test.traits$Migrant.Gained <- relevel(as.factor(perm.test.traits$Migrant.Gained), "Resident")
perm.test.traits$Migrant.Lost <- relevel(as.factor(perm.test.traits$Migrant.Lost), "Resident")

## add phylogeny
perm.test.traits <- perm.test.traits %>% left_join(phyl_dist)

# write data
write_csv(perm.test.traits, "../perm.test.traits.csv")
perm.test.traits <- read_csv("../perm.test.traits.csv")



#### Analyses

# dat.fit <- perm.test.traits %>% mutate(Obs_class2 = ifelse(adj.p.value < 0.05, 1, 0), 
#                                        pair = paste(Gained,"-", Lost),
#                                        diff = ifelse(Observed - floor(Perm.median) < 0, 0, Observed - floor(Perm.median) )) %>% 
#   filter(Gained != Lost) %>% na.omit(.)
dat.fit <- MuMIn::stdize(na.omit(perm.test.traits), prefix = F, 
                         omit.cols = c("Potential_lost" ,"Observed_lost", "Migrant.Gained", "Migrant.Lost")) %>% unique()
scaleList <- list(scale = attr(dat.fit, "scaled:scale"),
                  center = attr(dat.fit, "scaled:center"))


m1 <- glmer(cbind(Observed_lost, Potential_lost) ~ funct_dist + phyl_dist + (1|Gained) + (1|Lost),
            family = binomial, 
            data = dat.fit,
            nAGQ=0,control = glmerControl(optimizer = "nloptwrap",optCtrl=list(maxfun=1e10),calc.derivs = FALSE), 
            verbose = F, na.action = na.fail)

car::vif(m1)
car::Anova(m1)
summary(m1)

visreg(m1, xvar = "funct_dist", scale = "response")

pm1 <- sjPlot::plot_model(m1, ci.lvl = .95, axis.lim = c(0.8,1.2), vline.color = "grey", type = "std2",
                          show.value = T, value.offset = .1, title = "",
                          order.terms = c(1,2),
                          axis.labels = rev(c("Functional distance",
                                              "Phylogenetic distance"))) +
  theme_classic()

#dredge(m1)
# (Int)  fnc_dst phy_dst df    logLik    AICc delta weight
# -4.614 0.04198          4 -22649.51 45307.0  0.00  0.387
# -4.611 0.02822 0.02596  5 -22648.72 45307.4  0.43  0.312
# -4.611         0.03954  4 -22649.85 45307.7  0.70  0.273
# -4.610                  3 -22653.12 45312.2  5.23  0.028

m <- glmer(cbind(Observed_lost, Potential_lost) ~ Common_habitats + Common_diet +
             Habitat_width.Gained*Habitat_width.Lost + 
             weight.Gained*weight.Lost + 
             repro.Gained*repro.Lost + 
             Diet_width.Gained*Diet_width.Lost + 
             Migrant.Gained*Migrant.Lost + 
             (1|Gained) + (1|Lost),
           family = binomial, 
           data = dat.fit,
           nAGQ=0,control = glmerControl(optimizer = "nloptwrap",optCtrl=list(maxfun=1e10),calc.derivs = FALSE), 
           verbose = F, na.action = na.fail)

blmeco::dispersion_glmer(m)
car::vif(m)
car::Anova(m)


pm <- sjPlot::plot_model(m, ci.lvl = .95, axis.lim = c(0.05,12), vline.color = "grey", type = "std2",
                         show.value = T, value.offset = .4, title = "",
                         order.terms = c(1:17),
                         axis.labels = rev(c("No. habitats in common",
                                             "No. diet elements in common",
                                             "Habitat niche width of colonising species",
                                             "Habitat niche width of extinct species",
                                             "Body mass of colonising species",
                                             "Body mass of exinct species",
                                             "Reproduction of colonising species",
                                             "Reproduction of extinct species",
                                             "Diet niche width of colonising species",
                                             "Diet niche width of extinct species",
                                             "Migratory behaviour of colonising species",
                                             "Migratory behaviour of extinct species",
                                             "Habitat niche width of \ncolonising x extinct species",
                                             "Body mass of colonising x extinct species",
                                             "Reproduction of colonising x extinct species",
                                             "Diet niche width of \ncolonising x extinct species",
                                             "Migratory behaviour of colonising x extinct species"))) +
  theme_classic()

plot_grid(pm1,pm)


par(mfrow=c(2,2))

visreg(m, xvar = "Common_habitats", scale = "response", rug = F, band = F, ylim = c(0.02,.1),
       ylab = "Probability of negative interaction", xlab = "No. habitat in common",
       xtrans = function(x){x * scaleList$scale["Common_habitats"] + scaleList$center["Common_habitats"]})

visreg(m, xvar = "Common_diet", scale = "response", rug = F, band = F, ylim = c(0.02,.1), 
       ylab = "Probability of negative interaction", xlab = "No. diet items in common",
       xtrans = function(x){x * scaleList$scale["Common_diet"] + scaleList$center["Common_diet"]})

visreg(m, xvar = "Habitat_width.Gained", by = "Habitat_width.Lost", scale = "response", 
       rug = F, breaks = c(-1,0,1), overlay = T, band = F, ylim = c(0.02,.1),
       ylab = "Probability of negative interaction", xlab = "Habitat width of colonising species",
       xtrans = function(x){x * scaleList$scale["Habitat_width.Gained"] + scaleList$center["Habitat_width.Gained"]}, legend = F)
legend("topright", lty = 1, col = c("red", "green4", "blue"), lwd = 2,
       title = "Habitat width of extinct species", legend = c("Low", "Medium", "High"), bty = "n")

visreg(m, xvar = "Migrant.Gained", by = "Migrant.Lost", scale = "response", 
       rug = F, overlay = T, band = F, ylim = c(0.02,.1),
       ylab = "Probability of negative interaction", xlab = "Migratory behaviour of colonising species", legend = F)
legend("bottomright", lty = 1, col = c("red", "blue"), lwd = 2,
       title = "Migratory behaviour of extinct species", legend = c("Resident", "Migrant"), bty = "n")


d <- dredge(m)

###########
form <- as.formula(cbind(Observed_lost, Potential_lost) ~ Common_habitats + Common_diet +
                     Habitat_width.Gained*Habitat_width.Lost + 
                     weight.Gained*weight.Lost + 
                     repro.Gained*repro.Lost + 
                     Diet_width.Gained*Diet_width.Lost)

mods_per_sp <- dat.fit %>% group_by(Gained) %>% do(m = glm(form, data = ., family = binomial)) %>% broom::tidy(m)
