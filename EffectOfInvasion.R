library(emmeans)
library(lmerTest)
library(visreg)
library(tidyverse)

data <- readRDS("../Data/Butterflies - Netherlands/AllSpecies_reg_gam_ind_20171206_algroutes.rds")

data %>% filter(YEAR > 1991, SPECIES == "landkaartje", regional_gam > 0) %>% mutate(SITE = as.factor(SITE)) %>%
  group_by(YEAR) %>% summarise(nSites = n()) %>%
  ggplot(aes(y = nSites, x = YEAR)) + geom_line() + geom_point() + theme(legend.position="none") + geom_smooth()

data %>% group_by(YEAR) %>% filter(YEAR > 1991) %>%
  summarise(nSites = n(), nSite_occ = length(unique(SITE[SPECIES == "landkaartje" & regional_gam > 0]))) %>%
  mutate(prop_occ = nSite_occ / nSites) %>%  
  ggplot(aes(y = prop_occ, x = YEAR)) + geom_line() + geom_point() + theme(legend.position="none") + geom_smooth()


countsFIN1 <- read.table("../Data/Butterflies - Finland/FINLAND_Records_1999-2015.txt", sep = ";", h=T)
countsFIN2 <- read.table("../Data/Butterflies - Finland/FINLAND_Records_2016.txt", sep = ";", h=T)[,-6]
countsFIN <- rbind(countsFIN1, countsFIN2)
data2 <- countsFIN %>% group_by(Species_Faunaeur, Year, Site) %>% summarise(Individuals = max(Individuals))

data2 %>% group_by(Year) %>% 
  summarise(nSites = n(), nSite_occ = length(unique(Site[Species_Faunaeur == "Araschnia levana" & Individuals > 0]))) %>%
  mutate(prop_occ = nSite_occ / nSites) %>%  
  ggplot(aes(y = prop_occ, x = Year)) + geom_line() + geom_point() + theme(legend.position="none") + geom_smooth()


res <- c()
for(i in unique(data2$Species_Faunaeur)){
  print(i)
  
  dat.temp <- data2 %>% filter(Species_Faunaeur == i) %>% ungroup()
  
    data2_invMap <- data2  %>% group_by(Site) %>% filter(Species_Faunaeur == "Araschnia levana") %>% 
      select(-1) %>% summarise(inv.Year = min(Year)) %>%
      right_join(dat.temp, by = c("Site")) %>% filter(Species_Faunaeur == i) %>%
      mutate(Invasion = ifelse(Year > inv.Year, "Yes", "No"),
             Year = Year - min(Year)) %>% 
      mutate(Invasion = factor(ifelse(is.na(Invasion), "No", Invasion), levels = c("No", "Yes")))
    
    
    if(length(unique(data2_invMap$Invasion)) > 1 & length(unique(data2_invMap$Individuals)) > 1 & 
       length(unique(data2_invMap$Invasion)) * length(unique(data2_invMap$Site)) < nrow(data2_invMap)){
    
        m <- glmer(Individuals ~ Year*Invasion + (Year|Site), family = "poisson", 
                   data = data2_invMap,
                   nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                                 optCtrl=list(maxfun=1e10),
                                                 calc.derivs = FALSE), verbose = F)
      
      trends <- as.data.frame(emtrends(m, ~ Invasion, "Year", transform = "response"))
      
      res <- rbind.data.frame(res, cbind.data.frame(Species = i, 
                                                    nData = length(unique(data2_invMap$Invasion)) * length(unique(data2_invMap$Site)), 
                                                    trends))
    }
}

ggplot(res, aes(y = Year.trend, x = Invasion)) + geom_point() + 
  geom_errorbar(aes(ymin = Year.trend- SE, ymax = Year.trend + SE))

mm <- lmer(Year.trend ~ Invasion + (1|Species), weight = 1/nData, data  = res%>%filter(nData > 50))
summary(mm)
visreg(mm, xvar = "Invasion", scale = "response")


dat.temp <-  data2_invMap %>% mutate(Site = as.factor(Site), 
                                     Year = Year - min(Year))

m <- gamm(Individuals ~ s(Year, by = A_levana), random = list(Site = ~ Year), family = "poisson", 
          data = dat.temp)
summary(m$gam)
m$gam$data <-dat.temp 
visreg(m$gam, xvar = "Year", by = "A_levana", breaks = c(0,5), rug = F, trans = exp)

m <- glmer(Individuals ~ Year*A_levana + (Year|Site), family = "poisson", 
           data = dat.temp,
           nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                         optCtrl=list(maxfun=1e10),
                                         calc.derivs = FALSE), verbose = F)
summary(m)
visreg(m, xvar = "Year", by = "A_levana", breaks = c(0,5), rug = F, trans = exp)


pairs(emtrends(m, ~ A_levana, "Year", transform = "response"))


