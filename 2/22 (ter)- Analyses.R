library(tidyverse)
library(broom)
library(foreach)
library(vegan)
library(betapart)
library(lmerTest)
library(MuMIn)
library(cowplot)
library(visreg)

data <- read.csv("../SWE_PA_data.csv") %>% filter(latin != "Branta leucopsis")

trends_per_species <- as.tbl(read.csv("../trends_per_species.csv"))
decay.jackknife <- as.tbl(read.csv("../decay.jackknife.csv"))
std.decay <- as.tbl(read.csv("../std.decay.csv"))
decay <- as.tbl(read.csv("../decay.csv"))

# prepare traits
traits <- read_csv("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Birds trait database/traits_birds_corrected.csv")
sti <- read_csv("../Data/Birds Atlas/EBCC1/EBCC1_sti.csv")

traits2 <- traits %>% mutate(Habitat_specialisation =  rowSums(.[c(53:67)]),
                             Diet_specialisation =  rowSums(.[c(77:85)]),
                             Reproduction = Clutch_MEAN * Broods.per.year,
                             Migrant  = Long.distance.migrant + Short.distance.migrant) %>% 
  select(-c(5:86), WingU_MEAN, WeightU_MEAN) %>% left_join(sti, by = c("Species" = "species")) %>% 
  right_join(data %>% select(latin) %>% unique(), by = c("Species" = "latin")) %>% 
  mutate(Migrant = ifelse(Migrant == 0, "Resident", ifelse(Migrant >= 1, "Migrant", NA)))


# cti trend
library(optimx)
cti.dat <- data %>% left_join(sti, by = c("latin" = "species")) %>% filter(yr > 1996) %>%
  group_by(karta, yr) %>% summarise(cti.ab = weighted.mean(sti, abundance, na.rm = T), 
                                    cti.pa = mean(sti, na.rm = T),
                                    x = unique(x), y = unique(y)) %>% ungroup() %>% mutate(yr = yr-1996, 
                                                                                           x = scale(x), y = scale(y))

cti.trend.ab <- lmer(cti.ab ~ yr + x*y + (1|karta), data = cti.dat)
cti.trend.pa <- lmer(cti.pa ~ yr + x*y + (1|karta), data = cti.dat)

summary(cti.trend.ab)
summary(cti.trend.pa)

# survey effort per year
survey_effort <- data %>% group_by(yr) %>% summarise(nSiteSurveyed = length(unique(karta)))

# species trends per species
trends_per_species <- foreach(i = unique(data$latin), .combine = rbind.data.frame) %do% {
  dat.temp <- data %>% filter(latin == i) %>% left_join(survey_effort)  %>% mutate(yr = yr - 1996)
  if(length(unique(dat.temp$karta)) > 1){
    m.ab <- glmer(abundance ~ yr + (1|karta), data = dat.temp, family = poisson)
  } else {
    m.ab <- glm(abundance ~ yr, data = dat.temp, family = poisson)
  }
  
  m.occ <- glm(cbind(nSite, nSiteSurveyed - nSite) ~ yr, data = dat.temp %>% filter(abundance > 0) %>% group_by(yr) %>% 
                 summarise(nSite = n(),
                           nSiteSurveyed = unique(nSiteSurveyed)), family = binomial)
  
  cbind.data.frame(Species = i, rbind.data.frame(cbind.data.frame(Type = "Abundance", broom::tidy(m.ab))[,1:6],
                                                 cbind.data.frame(Type = "Occupancy", broom::tidy(m.occ))))
}

write_csv(trends_per_species, "../trends_per_species.csv")

# overall species trends per site
trends_per_site <- data %>% group_by(karta) %>% 
  mutate(yr = yr - min(yr)) %>% do(trend = glmer(abundance ~ yr + (1|latin), data = ., family = poisson)) %>% broom::tidy(trend)

trends_per_site %>% left_join(data %>% group_by(karta) %>% summarise(x = unique(x), y = unique(y))) %>%
  filter(term == "yr") %>%
  ggplot(aes(x = x, y = y, fill = estimate)) + geom_tile()

# overall species richness trends per site
trends_SR_per_site <- data %>% filter(abundance > 0) %>% group_by(karta, yr) %>% 
  summarise(SR = length(unique(latin))) %>% do(trend = glm(SR ~ yr, data = ., family = poisson)) %>% broom::tidy(trend)

data %>% group_by(karta) %>% summarise(x = unique(x), y = unique(y)) %>% 
  left_join(trends_SR_per_site %>% filter(term == "yr")) %>%
  ggplot(aes(x = x, y = y, fill = estimate)) + geom_tile()

data %>% group_by(karta) %>% summarise(x = unique(x), y = unique(y)) %>% 
  left_join(trends_SR_per_site %>% filter(term == "yr")) %>%
  ggplot(aes(y = estimate, x = y)) + geom_point() + geom_smooth() + 
  scale_y_continuous("Trend in species richness") + scale_x_continuous("Latitude")


# temporal decay of similarity
decay <- foreach(i = unique(data$karta), .combine = rbind.data.frame) %do%{
  mod <- decay.model(beta.pair.abund(data %>% filter(karta == i) %>% select(2,3,5) %>% 
                                       spread(., latin, abundance, fill = 0))$beta.bray,
                     dist(unique(data %>% filter(karta == i) %>% select(2))$yr), 
                     y.type = "dissimilarities", perm=100)
  cbind.data.frame(Site = i, slope.decay = mod$b.slope, p.value = mod$p.value, rsq = mod$pseudo.r.squared)
}

decay <- decay %>% as.tbl %>% 
  left_join(data %>% group_by(karta) %>% summarise(x = unique(x), y = unique(y)), by = c("Site"  = "karta"))

write_csv(decay, "../decay.csv")

decay %>% ggplot(aes(x = x, y = y, fill = slope.decay)) + geom_tile()
decay %>% mutate(sig = ifelse(p.value < 0.05, "yes", "no")) %>% ggplot(aes(x = x, y = y, fill = p.value)) + geom_tile()
decay %>% mutate(sig = ifelse(p.value < 0.05, "yes", "no")) %>% ggplot(aes(x = sig, y = y)) + geom_boxplot()

data %>% group_by(karta) %>% summarise(x = unique(x), y = unique(y), SR = length(unique(latin))) %>% 
  left_join(decay) %>% 
  ggplot(aes(y = slope.decay, x = y)) + geom_point() + geom_smooth(method = "lm") +
  scale_y_continuous("Temporal decay of similarity") + scale_x_continuous("Latitude")

data %>% group_by(karta) %>% summarise(x = unique(x), y = unique(y), SR = length(unique(latin))) %>% 
  left_join(decay) %>% 
  ggplot(aes(y = slope.decay, x = SR)) + geom_point() + geom_smooth(method = "lm") +
  scale_y_continuous("Temporal decay of similarity") + scale_x_continuous("Species richness")

# standardize by species richness
m <- lm(slope.decay ~ SR, 
        data = data %>% group_by(karta) %>% summarise(x = unique(x), y = unique(y), SR = length(unique(latin))) %>% 
          left_join(decay))

std.decay <- bind_cols(data %>% group_by(karta) %>% 
                         summarise(x = unique(x), y = unique(y), SR = length(unique(latin))) %>% 
                         rename(Site = karta), std.decay = resid(m), decay = decay$slope.decay)
ggplot(std.decay, aes(x = x, y = y, fill = std.decay)) + geom_tile()
ggplot(std.decay, aes(x = y, y = std.decay)) + geom_point() + geom_smooth(method = "lm") +
  scale_y_continuous("Temporal decay of similarity") + scale_x_continuous("Latitude")

summary(lm(std.decay ~ SR + x*y, data = std.decay))

write_csv(std.decay, "../std.decay.csv")


# jackknife species
decay.jackknife <- foreach(j = unique(data$latin), .combine = rbind.data.frame) %:% 
  foreach(i = unique(data$karta), .combine = rbind.data.frame) %do%{
    temp.sp <- data %>% filter(karta == i) %>% select(2,3,5)
    if(j %in% temp.sp$latin){
      mod <- decay.model(beta.pair.abund(temp.sp %>% filter(latin != j) %>% 
                                           spread(., latin, abundance, fill = 0))$beta.bray,
                         dist(unique(data %>% filter(karta == i) %>% select(2))$yr), 
                         y.type = "dissimilarities", perm=100)
      cbind.data.frame(Species = j, Site = i, slope.decay = mod$b.slope, p.value = mod$p.value, rsq = mod$pseudo.r.squared)
    }
  }

write_csv(decay.jackknife, "../decay.jackknife.csv")

m2 <- lm(slope.decay ~ SR, 
         data = decay.jackknife %>% left_join(data %>% group_by(karta) %>% summarise(SR = length(unique(latin))), 
                                              by = c("Site" = "karta")))
std.decay.jackknife <- bind_cols(data %>% group_by(karta) %>% summarise(x = unique(x),
                                                                        y = unique(y),
                                                                        SR = length(unique(latin))) %>% 
                                   left_join(decay.jackknife, by = c("karta" = "Site")), std.decay = resid(m2))

sp_effect <- std.decay.jackknife %>% rename(std.decay.jackknife = std.decay) %>% left_join(std.decay) %>% 
  mutate(effect = (std.decay - std.decay.jackknife)) %>% as.tbl()

sp_effect.mean <- sp_effect %>% group_by(Species) %>% 
  summarise(effect = mean(effect), nSites = length(unique(karta)))

sp_effect.mean %>% ungroup() %>% mutate(Species = fct_reorder(Species, effect)) %>% 
  ggplot(aes(x = effect, y = Species)) + geom_point() + geom_vline(xintercept = 0)


#############################
## analyse decay by region ##
#############################

sum.std.decay <- std.decay %>% group_by(Region = cut_number(y, n = 3)) %>% 
  summarise(Richness = mean(SR), se.richness = plotrix::std.error(SR),
            decay = mean(std.decay), se.decay = plotrix::std.error(std.decay),
            n = n())
levels(sum.std.decay$Region) <- c("South", "Middle", "North")


std.decay <- std.decay %>% mutate(Region = cut_number(y, n = 3))
levels(std.decay$Region) <- c("South", "Middle", "North")

ggplot(std.decay, aes(x = Region, y = SR)) + 
  geom_boxplot() + geom_signif(comparisons = list(c("South", "North"), c("South", "Middle"), c("Middle", "North")))

ggplot(std.decay, aes(x = Region, y = std.decay)) + 
  geom_boxplot() + geom_signif(comparisons = list(c("South", "North"), c("South", "Middle"), c("Middle", "North")))


#####################################
## analyse species effect by trait ##
#####################################

sp_effect.mean <- sp_effect %>% group_by(Species, Region) %>% 
  summarise(effect = mean(effect), nSites = length(unique(karta)))

sp_effect.mean.traits <- sp_effect.mean %>% left_join(traits2) %>% na.omit

m.sp <- lm(effect ~ poly(Habitat_specialisation, 2) + poly(Diet_specialisation, 2) + 
             poly(Reproduction, 2) + Migrant +
             poly(WingU_MEAN, 2) + poly(WeightU_MEAN, 2) + poly(sti, 2), 
           weight = nSites, data = sp_effect.mean.traits, na.action = na.fail)
d.sp <- dredge(m.sp)
summary(model.avg(d.sp, subset = delta < 2))


visreg(lm(effect ~ poly(WingU_MEAN, 2) +  poly(sti, 2), 
          weight = nSites, data = sp_effect.mean.traits, na.action = na.fail), scale = "response")


###############################################
## analyse species effect by trait by site ##
###############################################

site.traitMean <- sp_effect %>% left_join(traits2) %>% group_by(karta) %>% 
  summarise(Habitat_specialisation.mean = mean(Habitat_specialisation, na.rm = T),
            Diet_specialisation.mean = mean(Diet_specialisation, na.rm = T),
            Reproduction.mean = mean(Reproduction, na.rm = T),
            WingU_MEAN.mean = mean(WingU_MEAN, na.rm = T),
            WeightU_MEAN.mean = mean(WeightU_MEAN, na.rm = T),
            sti.mean = mean(sti, na.rm = T))

sp_effect.site.traits <- sp_effect %>% left_join(traits2) %>% left_join(site.traitMean, by = "karta") %>%
  group_by(Species, karta) %>% mutate(Habitat_specialisation = Habitat_specialisation - Habitat_specialisation.mean,
                                      Diet_specialisation = Diet_specialisation - Diet_specialisation.mean,
                                      Reproduction = Reproduction - Reproduction.mean,
                                      WingU_MEAN = WingU_MEAN - WingU_MEAN.mean,
                                      sti = sti - sti.mean) %>% na.omit

m.sp.site <- lmer(effect ~ poly(Habitat_specialisation, 2) + poly(Diet_specialisation, 2) + 
                    poly(Reproduction, 2) + Migrant +
                    poly(WingU_MEAN, 2) + poly(WeightU_MEAN, 2) + poly(sti, 2) + (1|karta) + (1|Species),
                  weights = SR, 
                  data = sp_effect.site.traits, na.action = na.fail)

summary(m.sp.site)

d.sp.reg <- dredge(m.sp.reg)
avg.m.sp.reg  <- model.avg(d.sp.reg , subset = delta < 2)
summary(avg.m.sp.reg )

