library(codyn)
library(data.table)
library(lme4)
library(lmerTest)
library(MuMIn)
source("functions.R")
library(visreg)
library(lattice)
source("similarity_decay.R")

butterflies <- bind_rows(read_csv("../Data/pres_abs_FIN_data.csv"), read_csv("../Data/pres_abs_NL_data.csv"))
site_data <- butterflies %>% group_by(Site, Scale) %>% summarise(X = unique(X),
                                                                 Y = unique(Y),
                                                                 PLAND = unique(PLAND), 
                                                                 CLUMPY = unique(CLUMPY),
                                                                 Habitat = unique(Habitat),
                                                                 gridCell50 = unique(gridCell50))

## Finland ##

dat.FIN <- butterflies %>% dplyr::filter(country == "FIN", Scale == "3km", n == 1)

# extract turnover rates

# FIN_turnover <- similarity_decay(dat.FIN)

FIN_turnover <- turnover2(dat.FIN)

FIN_turnover <- left_join(FIN_turnover, site_data, by = "Site")
FIN_turnover <- stdize(FIN_turnover, prefix = F, 
                       omit.cols = c("total.turnover", "appearance", "disappearance", 
                                     "gridCell50", "total.species", "Scale", 
                                     "beta.sor", "beta.sne", "beta.sim"))

FIN_turnover.mean <- FIN_turnover %>% group_by(Site) %>% summarise_all(mean)

boxplot(beta.sor ~ Habitat, data = FIN_turnover, scales = list(relation = "free"))

summary(glm(total.turnover ~ PLAND * Habitat, weights = total.species, family = binomial, data = FIN_turnover))
summary(lm(disappearance ~ PLAND * Habitat, data = FIN_turnover.mean))
summary(lm(appearance ~ PLAND * Habitat, data = FIN_turnover.mean))



# 3000
# total turnover
m_FIN_total.turnover_3000 <- glmer(total.turnover ~ PLAND * CLUMPY + X * Y + (1|gridCell50/Site) + (1|Year), 
                                  family = binomial, weight = total.species,
                                  data = FIN_turnover %>% dplyr::filter(Scale == "3km"),
                                  nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                                                optCtrl=list(maxfun=1e10),
                                                                calc.derivs = FALSE))

summary(m_FIN_total.turnover_3000)

# appearance
m_FIN_appearance_3000 <- glmer(appearance ~ PLAND * CLUMPY + X * Y + (1|gridCell50/Site) + (1|Year), family = binomial, weight = total.species,
                              data = FIN_turnover %>% dplyr::filter(Scale == "3km"),
                              nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                                            optCtrl=list(maxfun=1e10),
                                                            calc.derivs = FALSE))

summary(m_FIN_appearance_3000)

# disappearance
m_FIN_disappearance_3000 <- glmer(disappearance ~ PLAND * CLUMPY + X * Y + (1|gridCell50/Site) + (1|Year), 
                                 family = binomial, weight = total.species,
                                 data = FIN_turnover %>% dplyr::filter(Scale == "3km"),
                                 nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                                               optCtrl=list(maxfun=1e10),
                                                               calc.derivs = FALSE))

summary(m_FIN_disappearance_3000)

# 30000
# total turnover
m_FIN_total.turnover_30000 <- glmer(total.turnover ~ PLAND * CLUMPY + X * Y + (1|gridCell50/Site) + (1|Year), 
                                   family = binomial, weight = total.species,
                                   data = FIN_turnover %>% dplyr::filter(Scale == "30km"),
                                   nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                                                 optCtrl=list(maxfun=1e10),
                                                                 calc.derivs = FALSE))

summary(m_FIN_total.turnover_30000)

# appearance
m_FIN_appearance_30000 <- glmer(appearance ~ PLAND * CLUMPY + X * Y + (1|gridCell50/Site) + (1|Year),
                               family = binomial, weight = total.species,
                               data = FIN_turnover %>% dplyr::filter(Scale == "30km"),
                               nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                                             optCtrl=list(maxfun=1e10),
                                                             calc.derivs = FALSE))

summary(m_FIN_appearance_30000)


# disappearance
m_FIN_disappearance_30000 <- glmer(disappearance ~ PLAND * CLUMPY + X * Y + (1|gridCell50/Site) + (1|Year), 
                                  family = binomial, weight = total.species,
                                  data = FIN_turnover %>% dplyr::filter(Scale == "30km"),
                                  nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                                                optCtrl=list(maxfun=1e10),
                                                                calc.derivs = FALSE))

summary(m_FIN_disappearance_30000)

visreg(m_FIN_appearance_30000, xvar = "PLAND", scale = "response")
visreg(m_FIN_appearance_3000, xvar = "PLAND", scale = "response")


## Netherlands ##

dat.NL <- butterflies %>% dplyr::filter(country == "NL", Scale == "3km", n == 1)

# extract turnover rates
NL_turnover <- turnover2(dat.NL)

NL_turnover <- left_join(NL_turnover, site_data, by = "Site")
NL_turnover <- stdize(NL_turnover, prefix = F, 
                       omit.cols = c("total.turnover", "appearance", "disappearance", "gridCell50", "total.species", "Scale"))

# 3000
# total turnover
m_NL_total.turnover_3000 <- glmer(total.turnover ~ PLAND * CLUMPY + X * Y + (1|gridCell50/Site) + (1|Year), 
                                  family = binomial, weight = total.species,
                                  data = NL_turnover %>% dplyr::filter(Scale == "3km"),
                                  nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                                                optCtrl=list(maxfun=1e10),
                                                                calc.derivs = FALSE))

summary(m_NL_total.turnover_3000)

# appearance
m_NL_appearance_3000 <- glmer(appearance ~ PLAND * CLUMPY + X * Y + (1|gridCell50/Site) + (1|Year), 
                              family = binomial, weight = total.species,
                              data = NL_turnover %>% dplyr::filter(Scale == "3km"),
                              nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                                            optCtrl=list(maxfun=1e10),
                                                            calc.derivs = FALSE))

summary(m_NL_appearance_3000)

# disappearance
m_NL_disappearance_3000 <- glmer(disappearance ~ PLAND * CLUMPY + X * Y + (1|gridCell50/Site) + (1|Year), 
                                 family = binomial, weight = total.species,
                                 data = NL_turnover %>% dplyr::filter(Scale == "3km"),
                                 nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                                               optCtrl=list(maxfun=1e10),
                                                               calc.derivs = FALSE))

summary(m_NL_disappearance_3000)

visreg(m_NL_appearance_3000, xvar = "PLAND", scale = "response")
visreg(m_NL_disappearance_3000, xvar = "PLAND", scale = "response")


# 30000
# total turnover
m_NL_total.turnover_30000 <- glmer(total.turnover ~ PLAND * CLUMPY + X * Y + (1|gridCell50/Site) + (1|Year), 
                                   family = binomial, weight = total.species,
                                   data = NL_turnover %>% dplyr::filter(Scale == "30km"),
                                   nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                                                 optCtrl=list(maxfun=1e10),
                                                                 calc.derivs = FALSE))

summary(m_NL_total.turnover_30000)

# appearance
m_NL_appearance_30000 <- glmer(appearance ~ PLAND * CLUMPY + X * Y + (1|gridCell50/Site) + (1|Year),
                               family = binomial, weight = total.species,
                               data = NL_turnover %>% dplyr::filter(Scale == "30km"),
                               nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                                             optCtrl=list(maxfun=1e10),
                                                             calc.derivs = FALSE))

summary(m_NL_appearance_30000)


# disappearance
m_NL_disappearance_30000 <- glmer(disappearance ~ PLAND * CLUMPY + X * Y + (1|gridCell50/Site) + (1|Year), 
                                  family = binomial, weight = total.species,
                                  data = NL_turnover %>% dplyr::filter(Scale == "30km"),
                                  nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                                                optCtrl=list(maxfun=1e10),
                                                                calc.derivs = FALSE))

summary(m_NL_disappearance_30000)

visreg(m_NL_disappearance_30000, xvar = "CLUMPY", by = "PLAND", scale = "response", breaks = 2)


###########################
## remove unused objects ##
###########################


gdata:::keep(FIN_turnover, NL_turnover,
             sure = T)
