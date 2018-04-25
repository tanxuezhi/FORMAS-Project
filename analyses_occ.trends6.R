source("functions.R")

dat <- sp_site_occTrend(dat.occ.trend, .2, .8)
dat <- left_join(dat, butterflies %>% dplyr::select(2,3,5,6,8,10:20,23) %>% group_by(Species, Site, Scale) %>% summarise_all(first))
dat <- stdize(dat, prefix = F, omit.cols = c("Trend", "Scale", "pred.then", "pred.now", "nYear"))
scaleList <- list(scale = attr(dat, "scaled:scale"),
                         center = attr(dat, "scaled:center"))

table(subset(dat, dat$Scale == 30000)$Trend)


dat_exp <- dat %>% dplyr::filter(Scale == 50000) %>% group_by(Species, Trend) %>% summarise(n = n()) %>% spread(Trend, n) %>%
  left_join(butterflies %>% dplyr::select(3,10:16,23) %>% group_by(Species) %>% summarise_all(first)) %>%
  dplyr::filter(!is.na(STI)) %>% mutate(Colonisation = ifelse(!is.na(Colonisation), Colonisation, 0),
                                        Extinction = ifelse(!is.na(Extinction), Extinction, 0),
                                        Persistence = ifelse(!is.na(Persistence), Persistence, 0))

dat_exp %>% ungroup() %>% arrange(desc(Colonisation)) %>% mutate(Species  = factor(Species, Species)) %>%
  gather(key = Event, value = value, Colonisation, Extinction) %>% 
  ggplot(aes(y = value, x = Species, fill = Event, group = Event)) + geom_bar(stat = "identity", position = "dodge") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(expand = c(0, 0))



dat_exp %>% ungroup() %>% arrange(desc(Colonisation)) %>% mutate(Species  = factor(Species, Species)) %>% ggplot(aes(y = Colonisation, x = Species, fill = Colonisation)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# landscape only #

dat.col <- dat %>% mutate(Colonisation = ifelse(Trend == "Colonisation", 1, ifelse(Trend == "Absence", 0, NA)))

m.col.land <- glmer(Colonisation ~ PLAND * CLUMPY + 
                 Habitat + X*Y + 
                 (1|gridCell50/Site) + (1|Species), 
               family = binomial,
               data = na.omit(subset(dat.col, dat.col$Scale == 10000)),
               nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                             optCtrl=list(maxfun=1e10),
                                             calc.derivs = FALSE), verbose = F, na.action = na.fail)
summary(m.col.land)

visreg(m.col.land, xvar = "CLUMPY", by = "PLAND", scale = "response", rug = F, breaks = c(-1,1),ylim = c(0,.5), overlay = T,
       xtrans =  function(x) x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"])

dat.ext <- dat %>% mutate(Extinction = ifelse(Trend == "Extinction", 1, ifelse(Trend == "Persistence", 0, NA)))

m.ext.land <- glmer(Extinction ~ PLAND * CLUMPY + 
                 Habitat + X*Y + 
                 (1|gridCell50/Site) + (1|Species), 
               family = binomial,
               data = na.omit(subset(dat.ext, dat.ext$Scale == 30000)),
               nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                             optCtrl=list(maxfun=1e10),
                                             calc.derivs = FALSE), verbose = F, na.action = na.fail)
summary(m.ext.land)

visreg(m.ext.land, xvar = "CLUMPY", by = "PLAND", scale = "response", rug = F, breaks = c(-1,1), ylim = c(0,.5),overlay = T,
       xtrans =  function(x) x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"])



# traits only #
dd <- na.omit(subset(dat.col, dat.col$Scale == 50000))
m.col.traits <- glmer(Colonisation ~ STI_rel * PC1 + 
                        STI_rel * PC3 + 
                        STI_rel * PC4 +
                      Habitat + X*Y + 
                      (1|gridCell50/Site) + (1|Species), 
                    family = binomial,
                    data = dd,
                    nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                                  optCtrl=list(maxfun=1e10),
                                                  calc.derivs = FALSE), verbose = F, na.action = na.fail)
summary(m.col.traits)

d.col.traits <- dredge(m.col.traits, fixed = c("X", "Y", "X:Y", "Habitat"))
summary(model.avg(d.col.traits, subset = delta < 4))


visreg(m.col.traits, xvar = "PC1", by = "STI_rel", scale = "response", rug = F, breaks = c(-1,0,1), overlay = T, 
       xtrans =  function(x) x * scaleList$scale["STI_rel"] + scaleList$center["STI_rel"])

visreg(m.col.traits, xvar = "STI_rel", by = "PC4", scale = "response", rug = F, breaks = c(-1,0,1), overlay = T, ,
       xtrans =  function(x) x * scaleList$scale["STI_rel"] + scaleList$center["STI_rel"])


m.ext.traits <- glmer(Extinction ~ STI_rel * PC1 + 
                        STI_rel * PC3 + 
                        STI_rel * PC4 + 
                      Habitat + X*Y + 
                      (1|gridCell50/Site) + (1|Species), 
                    family = binomial,
                    data = na.omit(subset(dat.ext, dat.ext$Scale == 30000)),
                    nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                                  optCtrl=list(maxfun=1e10),
                                                  calc.derivs = FALSE), verbose = F, na.action = na.fail)
summary(m.ext.traits)

visreg(m.ext.traits, xvar = "STI_rel", by = "PC1", scale = "response", rug = F, breaks = c(-1,0,1), overlay = T, 
       xtrans =  function(x) x * scaleList$scale["STI_rel"] + scaleList$center["STI_rel"])

visreg(m.ext.traits, xvar = "STI_rel", by = "PC3", scale = "response", rug = F, breaks = c(-1,0,1), overlay = T, 
       xtrans =  function(x) x * scaleList$scale["STI_rel"] + scaleList$center["STI_rel"])

visreg(m.ext.traits, xvar = "STI_rel", by = "PC4", scale = "response", rug = F, breaks = c(-1,0,1), overlay = T, 
       xtrans =  function(x) x * scaleList$scale["STI_rel"] + scaleList$center["STI_rel"])


# traits x landscape #

m.col.traitsLand <- glmer(Colonisation ~ STI_rel * PC1 + 
                            STI_rel * PC3 + 
                            STI_rel * PC4 +
                            PC1 * PLAND * CLUMPY + 
                        Habitat + X*Y + 
                        (1|gridCell50/Site) + (1|Species), 
                      family = binomial,
                      data = na.omit(subset(dat.col, dat.col$Scale == 30000)),
                      nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                                    optCtrl=list(maxfun=1e10),
                                                    calc.derivs = FALSE), verbose = F, na.action = na.fail)

summary(m.col.traitsLand)

visreg(m.col.traitsLand, xvar = "PC1", by = "PLAND", scale = "response", breaks = c(-1,0,1), rug =F, overlay = T, ylim = c(0,.5),
       xtrans =  function(x) x * scaleList$scale["PC1"] + scaleList$center["PC1"])



par(mfrow=c(2,2))

visreg(m.col.traitsLand, xvar = "CLUMPY", by = "PLAND", scale = "response", rug = F, breaks = c(-1,1), overlay = T,ylim = c(0,.25),
       xtrans =  function(x) x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"], cond = list(PC1 = -1, STI_rel = -1))

visreg(m.col.traitsLand, xvar = "CLUMPY", by = "PLAND", scale = "response", rug = F, breaks = c(-1,1),overlay = T,ylim = c(0,.25),
       xtrans =  function(x) x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"], cond = list(PC1 = -1, STI_rel = 1))

visreg(m.col.traitsLand, xvar = "CLUMPY", by = "PLAND", scale = "response", rug = F, breaks = c(-1,1),overlay = T,ylim = c(0,.25),
       xtrans =  function(x) x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"], cond = list(PC1 = 1, STI_rel = -1))

visreg(m.col.traitsLand, xvar = "CLUMPY", by = "PLAND", scale = "response", rug = F, breaks = c(-1,1),overlay = T,ylim = c(0,.25),
       xtrans =  function(x) x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"], cond = list(PC1 = 1, STI_rel = 1))



m.ext.traitsLand <- glmer(Extinction ~ STI_rel * PC1 + 
                            STI_rel * PC3 + 
                            STI_rel * PC4 +
                            PC1 * PLAND * CLUMPY + 
                            I(PC1^2) * PLAND * CLUMPY + 
                            Habitat + X*Y + 
                            (1|gridCell50/Site) + (1|Species), 
                          family = binomial,
                          data = na.omit(subset(dat.ext, dat.ext$Scale == 30000)),
                          nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                                        optCtrl=list(maxfun=1e10),
                                                        calc.derivs = FALSE), verbose = F, na.action = na.fail)

summary(m.ext.traitsLand)

par(mfrow=c(2,2))

visreg(m.ext.traitsLand, xvar = "CLUMPY", by = "PLAND", scale = "response", rug = F, breaks = c(-1,1), overlay = T,ylim = c(0,.5),
       xtrans =  function(x) x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"], cond = list(PC1 = -1, STI_rel = -1))

visreg(m.ext.traitsLand, xvar = "CLUMPY", by = "PLAND", scale = "response", rug = F, breaks = c(-1,1),overlay = T,ylim = c(0,.5),
       xtrans =  function(x) x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"], cond = list(PC1 = -1, STI_rel = 1))

visreg(m.ext.traitsLand, xvar = "CLUMPY", by = "PLAND", scale = "response", rug = F, breaks = c(-1,1),overlay = T,ylim = c(0,.5),
       xtrans =  function(x) x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"], cond = list(PC1 = 1, STI_rel = -1))

visreg(m.ext.traitsLand, xvar = "CLUMPY", by = "PLAND", scale = "response", rug = F, breaks = c(-1,1),overlay = T,ylim = c(0,.5),
       xtrans =  function(x) x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"], cond = list(PC1 = 1, STI_rel = 1))



###########################
## remove unused objects ##
###########################


gdata:::keep(butterflies, dat.occ.trend,
             sure = T)





