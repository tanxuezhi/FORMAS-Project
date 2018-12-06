library(tidyverse)
library(doFuture)
library(future.apply)
library(lmerTest)


## load data ##
# monitoring data
data <- read_csv("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Temporal_climatic_debt/SWE_PA_data.csv")
data <- data %>% filter(latin != "Branta leucopsis")
#traits data
traits <- as.tbl(read.table("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Birds trait database/Life-history characteristics of European birds.txt", sep = "\t", h = T))

data %>% filter(!latin %in% traits$Species) %>% distinct(latin) %>% print(n = 30)

traits$Species <- gsub("Spatula clypeata", "Anas clypeata", traits$Species)
traits$Species <- gsub("Mareca penelope", "Anas penelope", traits$Species)
traits$Species <- gsub("Spatula querquedula", "Anas querquedula", traits$Species)
traits$Species <- gsub("Mareca strepera", "Anas strepera", traits$Species)
traits$Species <- gsub("Linaria cannabina", "Carduelis cannabina", traits$Species)
traits$Species <- gsub("Acanthis flammea hornemanni", "Carduelis hornemanni", traits$Species)
traits$Species <- gsub("Acanthis flammea", "Carduelis flammea", traits$Species)
traits$Species <- gsub("Linaria flavirostris", "Carduelis flavirostris", traits$Species)
traits$Species <- gsub("Spinus spinus", "Carduelis spinus", traits$Species)
traits$Species <- gsub("Eudromias morinellus", "Charadrius morinellus", traits$Species)
traits$Species <- gsub("Corvus corone", "Corvus corone cornix", traits$Species)
traits$Species <- gsub("Dryobates minor", "Dendrocopos minor", traits$Species)
traits$Species <- gsub("Lagopus muta", "Lagopus mutus", traits$Species)
traits$Species <- gsub("Hydrocoloeus minutus", "Larus minutus", traits$Species)
traits$Species <- gsub("Calidris falcinellus", "Limicola falcinellus", traits$Species)
traits$Species <- gsub("Cyanecula svecica", "Luscinia svecica", traits$Species)
traits$Species <- gsub("Mergellus albellus", "Mergus albellus", traits$Species)
traits$Species <- gsub("Poecile cinctus", "Parus cinctus", traits$Species)
traits$Species <- gsub("Poecile montanus", "Parus montanus", traits$Species)
traits$Species <- gsub("Poecile palustris", "Parus palustris", traits$Species)
traits$Species <- gsub("Calidris pugnax", "Philomachus pugnax", traits$Species)
traits$Species <- gsub("Regulus ignicapilla", "Regulus ignicapillus", traits$Species)
traits$Species <- gsub("Sternula albifrons", "Sterna albifrons", traits$Species)
traits$Species <- gsub("Hydroprogne caspia", "Sterna caspia", traits$Species)
traits$Species <- gsub("Thalasseus sandvicensis", "Sterna sandvicensis", traits$Species)
traits$Species <- gsub("Morus bassanus", "Sula bassana", traits$Species)
traits$Species <- gsub("Lyrurus tetrix", "Tetrao tetrix", traits$Species)

traits <- bind_rows(traits, traits %>% filter(grepl("Loxia", traits$Species)) %>% 
                      summarise_all(mean) %>% mutate(Species = "Loxia species"))

data %>% filter(!latin %in% traits$Species) %>% distinct(latin) %>% print(n = 30)



## calculate turnover per site based on two time periods
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
  rename(sp.first = Species, sp.last = previous) %>% select(-2)

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

## compute association between colonisation and extinctions
registerDoFuture()
plan(multiprocess)

Assoc <- foreach(i = unique(names(gained.sp)), .combine = rbind.data.frame) %dopar%{
  turnover.per.site[unlist(lapply(turnover.per.site$gained, 
                                  function(x){any(x$latin %in% i)})),] %>% 
    filter(n.lost >0) %>%
    rowwise() %>% do(l = cbind(Site = .$karta, Gained = i, .$lost)) %>% 
    unnest() %>% rename(Lost = latin)
}

Assoc <- rbind.data.frame(Assoc, cbind.data.frame(Site = NA, Gained = NA, 
                                                  Lost = setdiff(unique(Assoc$Gained), unique(Assoc$Lost))))
Assoc <- rbind.data.frame(Assoc, cbind.data.frame(Site = NA, Lost = NA, 
                                                  Gained = setdiff(unique(Assoc$Lost), unique(Assoc$Gained))))

# compute table of number of associations 
assoc.sum <- table(Assoc[2:3])

# permutation test
avg.inter <- future_replicate(1000, {
  perm.dat <- bind_cols(Assoc[,1:2], Lost = sample(Assoc$Lost))
  perm.dat.sum <- table(perm.dat[2:3])
}
)

# extract p-values for pairwise associations
plan(multiprocess)

p.values <- future_apply(abind(avg.inter, assoc.sum), c(1,2), 
                         function(x){(sum(x[1:1000] >= x[1001]) + 1) / (length(x[1:1000]) + 1)})

obs.quantile = future_apply(abind(avg.inter, assoc.sum), c(1,2), 
                            function(x){quantile(x[1:1000], probs = c(.05,.5,.95))})


perm.test <- as.data.frame(as.table(p.values)) %>% rename(Gained = Var1, Lost = Var2, p.value = Freq) %>%
  left_join(as.data.frame(as.table(obs.quantile)) %>% rename(Gained = Var2, Lost = Var3) %>%
              spread(key = Var1, value = Freq)) %>% rename(Perm.median = `50%`, LowCI = `5%`, UprCI = `95%`) %>%
  left_join(as.data.frame(as.table(assoc.sum)) %>% rename(Observed = Freq)) %>% 
  as.tbl() %>% select(1,2,7,5,4,6,3) %>% mutate(adj.p.value = p.adjust(p.value, method = "fdr"))


write_csv(perm.test, "../perm.test.csv")
perm.test <- read_csv("../perm.test.csv")

perm.test %>% filter(adj.p.value < .05)
perm.test %>% filter(p.value < .05)


# example plot
i = "Jynx torquilla"          
j = "Chloris chloris"
p.value = (sum(avg.inter[i,j,] >= assoc.sum[i,j]) + 1) / (length(avg.inter[i,j,]) + 1)
hist(avg.inter[i,j,], xlim = c(0, max(assoc.sum[i,j], max(avg.inter[i,j,]))+1), main = p.value)
abline(v=assoc.sum[i,j], col="red", lwd=2)

# add traits
registerDoFuture()
plan(multiprocess)


add_traits <- function(x){
  sp1 <- as.character(x["Gained"][[1]])
  sp2 <- as.character(x["Lost"][[1]])
  
  return(
    cbind.data.frame(diff_weight = c(traits[traits$Species == sp1,"WeightU_MEAN"] -
                                       traits[traits$Species == sp2,"WeightU_MEAN"])[[1]],
                     diff_repro = ((traits[traits$Species == sp1,"Clutch_MEAN"] *
                                      traits[traits$Species == sp1,"Broods.per.year"]) -
                                     (traits[traits$Species == sp2,"Clutch_MEAN"] *
                                        traits[traits$Species == sp2,"Broods.per.year"]))[[1]],
                     diff_habitat = dist(traits[traits$Species %in% c(sp1, sp2),53:67])[1],
                     diff_habitat_width = sum(traits[traits$Species == sp1,53:67]) -
                       sum(traits[traits$Species == sp2,53:67]),
                     diff_diet = dist(traits[traits$Species %in% c(sp1, sp2),68:85])[1],
                     diff_diet_width = sum(traits[traits$Species == sp1,68:85]) -
                       sum(traits[traits$Species == sp2,68:85])
    )
  )
}

perm.test.traits <- as.tbl(cbind.data.frame(perm.test, 
                                            data.table::rbindlist(future_apply(perm.test, 1, add_traits))))

write_csv(perm.test.traits, "../perm.test.traits.csv")
perm.test.traits <- read_csv("../perm.test.csv")


perm.test.traits %>% mutate(Observed = ifelse(Observed == 0, "no", "yes")) %>% 
  ggplot(aes(y = diff_weight, x = Observed, fill = Observed)) + geom_violin(draw_quantiles = c(.05,.5,.95)) + theme_classic()



dat.fit <- perm.test.traits %>% mutate(Obs_class = ifelse(Observed == 0, 0, 1),
                                       Obs_class2 = ifelse(p.value < 0.05, 1, 0), 
                                       pair = paste(as.numeric(Gained), as.numeric(Lost)))
dat.fit <- MuMIn::stdize(dat.fit, prefix = F, omit.cols = c("Obs_class", "Obs_class2"))

m <- glmer(Obs_class2 ~ diff_weight + diff_repro + diff_habitat_width + diff_diet_width + (1|Lost) + (1|Gained),
           family = binomial, 
           data = dat.fit,
           nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                         optCtrl=list(maxfun=1e10),
                                         calc.derivs = FALSE), verbose = F)
car::Anova(m)


m_tot001 <- lmer(values ~ p.cat + (1|Gained) + (1|Lost), 
                 data = perm.test.traits %>% filter(values != 0) %>% 
                   mutate(p.cat = cut(p.value, breaks = c(0,0.001,1))))
m_tot01 <- lmer(diff_weight ~ p.cat + (1|Gained) + (1|Lost), 
                data = perm.test.traits %>% 
                  mutate(p.cat = cut(p.value, breaks = c(0,0.01,1)), diff_weight = scale(diff_weight)))
m_tot05 <- lmer(values ~ p.cat + (1|Gained) + (1|Lost), 
                data = perm.test.traits %>% filter(values != 0) %>% 
                  mutate(p.cat = cut(p.value, breaks = c(0,0.05,1))))

dat.fit <- perm.test.traits %>% 
  mutate(p.cat = cut(p.value, breaks = c(0,0.01,1)), diff_weight = scale(diff_weight), 
         pair = paste(as.numeric(Gained), as.numeric(Lost))) %>%
  select(9,15,16)

summary(aov(diff_weight ~ p.cat + Error(pair), data = dat.fit))

car::Anova(m_tot01, "III")
pairs(emmeans::emmeans(m_tot05, ~ p.cat))
plot(emmeans::emmeans(m_tot05, ~ p.cat))
effect_all <- as.data.frame(pairs(emmeans::emmeans(m_tot, ~ p.cat)))$estimate

library(gamlss)
dat.fit <- perm.test.traits %>% filter(values != 0) %>% 
  mutate(Lost = as.factor(Lost), Gained = as.factor(Gained),
         effect = Observed - Perm.median) %>% dplyr::select(-adj.p.value)

m2 <- gamlss(p.value ~  values + random(Gained) + random(Lost), 
             family = BEOI, 
             data = dat.fit)
summary(m2)

term.plot(m2)

visreg::visreg(m2, scale = "response")

## jackknife analysis
registerDoFuture()
plan(multiprocess, workers = 3)


res_jack_traits <- foreach(i = 5:85, .combine = rbind.data.frame) %dopar% {
  
  dist_temp <- FD::gowdis(traits %>% select(setdiff(5:85, i)) %>% as.data.frame)
  dist_temp <- as.matrix(dist_temp)
  rownames(dist_temp) <- traits$Species; colnames(dist_temp) <- traits$Species
  dist_temp <- as.tbl(data.frame(rows=rownames(dist_temp)[row(dist_temp)], vars=colnames(dist_temp)[col(dist_temp)],
                                 values=c(dist_temp)))
  
  perm.test.traits.temp <- perm.test %>% left_join(dist_temp, by = c("Gained" = "rows", "Lost" = "vars"))
  
  m <- lmer(values ~ p.cat + (1|Gained) + (1|Lost), weights = Observed, 
            data = perm.test.traits.temp %>% filter(values != 0) %>% 
              mutate(p.cat = cut(p.value, breaks = c(0,0.01,1))))
  cbind.data.frame(trait = names(traits)[i], effect = as.data.frame(pairs(emmeans::emmeans(m, ~ p.cat)))$estimate)
  
}

res_jack_traits <- res_jack_traits %>%
  mutate(diff_effect = (effect - effect_all)/effect_all) %>% as.tbl()

res_jack_traits %>% mutate(trait = forcats::fct_reorder(trait, diff_effect)) %>% 
  ggplot(aes(x = diff_effect, y = trait)) + geom_point() + geom_abline()
