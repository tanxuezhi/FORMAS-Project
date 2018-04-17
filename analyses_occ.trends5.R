library(data.table)
library(visreg)
library(MuMIn)
library(lmerTest)
library(doFuture)
library(polypoly)
source("functions.R")

butterflies <- as.tbl(fread("C:/Local Folder (c)/butterflies_occ.csv"))

dat.occ.trend <- sp_site_occupancy_trend.predict(butterflies %>% group_by(Species, Site, Year) %>% summarise_all(first))
write_csv(dat.occ.trend,"../sp_site_trend_pred.csv")

dat.occ.trend <- read_csv("../sp_site_trend_pred.csv")

registerDoFuture()
options(future.globals.maxSize = Inf)
plan(multiprocess, workers = 2)

res <- foreach(i = c(1:4)/10, j = rev(c(6:9)/10), .combine = rbind.data.frame) %:% foreach(k = unique(butterflies$Scale), .combine = rbind.data.frame) %dopar%{
  
  dat.ColExt <- sp_site_colExt(dat.occ.trend, i, j)
  dat.ColExt <- left_join(dat.ColExt, butterflies %>% dplyr::select(2,3,5,6,8,10:20,23) %>% group_by(Species, Site, Scale) %>% summarise_all(first))
  dat.ColExt <- stdize(dat.ColExt, prefix = F, omit.cols = c("Extinction", "Colonisation", "Scale"))
  
  m.col <- glmer(Colonisation ~ STI_rel * PC3 + 
                   STI_rel * PC4 + 
                   STI_rel * poly_rescale(poly(PC1,2),2) + 
                   poly_rescale(poly(PC1,2),2) * PLAND * CLUMPY + 
                   Habitat + X*Y + 
                   (1|gridCell50/Site) + (1|Species), 
                 family = binomial,
                 data = subset(dat.ColExt, dat.ColExt$Scale == k),
                 nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                               optCtrl=list(maxfun=1e10),
                                               calc.derivs = FALSE), verbose = F, na.action = na.fail)
  m.ext <- glmer(Extinction ~  STI_rel * PC3 + 
                   STI_rel * PC4 + 
                   STI_rel * poly_rescale(poly(PC1,2),2) + 
                   poly_rescale(poly(PC1,2),2) * PLAND * CLUMPY + 
                   Habitat + X*Y + 
                   (1|gridCell50/Site) + (1|Species), 
                 family = binomial,
                 data = subset(dat.ColExt, dat.ColExt$Scale == k),
                 nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                               optCtrl=list(maxfun=1e10),
                                               calc.derivs = FALSE), verbose = F)
  
  rbind.data.frame(cbind.data.frame(Scale = k, 
                                    Process = "Colonisation",
                                    Range = paste(i,"-", j), 
                                    estimate = fixef(m.col), 
                                    confint(m.col, method = "Wald")[-c(1:3),]) %>% rownames_to_column("Variables"),
                   cbind.data.frame(Scale = k, 
                                    Process = "Extinction",
                                    Range = paste(i,"-", j), 
                                    estimate = fixef(m.ext), 
                                    confint(m.ext, method = "Wald")[-c(1:3),]) %>% rownames_to_column("Variables"))
}

names(res)[6:7] <- c("lwr", "upr")

write_csv(res, "../coefs_mods_ColExt.csv")

res <- read_csv("../coefs_mods_ColExt.csv")
# res$Scale <- as.factor(res$Scale)
# levels(res$Scale) <- paste(c(1,3,5,10,20,30,50), "km")
# 
# res$Variables <- gsub("poly_rescale\\(poly\\(PC1, 2\\), 2\\)2", "PC1^2", res$Variables)
# res$Variables <- gsub("poly_rescale\\(poly\\(PC1, 2\\), 2\\)1", "PC1", res$Variables)
# res$Variables <- gsub("STI_rel", "STI", res$Variables)
# res$Variables <- gsub("PLAND", "%SNH", res$Variables)
# res$Variables <- gsub("CLUMPY", "Clumpiness", res$Variables)
# res$Variables <- gsub(":", " x ", res$Variables)
# res <- res %>% mutate(Significancy = ifelse(lwr < 0 & upr < 0 | lwr > 0 & upr > 0, "Significant", "Not significant"))


ggplot(res %>% dplyr::filter(grepl("STI_rel|PC1|PC3|PC4", Variables)), 
       aes(x = Variables, y = estimate, color = Scale, alpha = Significancy)) + 
  geom_point(position = position_dodge(.8)) + facet_grid(Process ~ Range) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(.8)) + 
  coord_flip() +  scale_color_brewer(palette="YlOrRd")


#### turnover ####
dat.turn <- subset(dat.ColExt, dat.ColExt$Scale == 30000) %>% group_by(Site) %>% 
  summarise(Turnover = sum(Colonisation, Extinction),
            Extinction = sum(Extinction),
            Colonisation = sum(Colonisation),
            nSpecies = n(),
            PLAND = mean(PLAND),
            CLUMPY = mean(CLUMPY),
            X = unique(X), Y = unique(Y), gridCell50 = unique(gridCell50))

summary(dat.turn[, "Extinction"])
summary(dat.turn[, "Colonisation"])


m.turn <- glm(cbind(Turnover, nSpecies - Turnover) ~ PLAND * CLUMPY + X*Y, family = binomial, data = dat.turn)

summary(m.turn)

visreg2d(m.turn, xvar = "PLAND", yvar = "CLUMPY", scale = "response", plot.type = "image")


m.ext.num <- glm(cbind(Extinction, nSpecies - Extinction) ~ PLAND * CLUMPY + X*Y, family = binomial, data = dat.turn)

summary(m.ext.num)

visreg2d(m.ext.num, xvar = "PLAND", yvar = "CLUMPY", scale = "response", plot.type = "image")

m.col.num <- glm(cbind(Colonisation, nSpecies - Colonisation) ~ PLAND * CLUMPY + X*Y, family = binomial, data = dat.turn)

summary(m.col.num)

visreg2d(m.col.num, xvar = "PLAND", yvar = "CLUMPY", scale = "response", plot.type = "image")


### set scale at 20km and range at 0.2 - 0.8 ###

##################
## colonisation ##
##################

dat.ColExt <- sp_site_colExt(dat.occ.trend, .2, .8)
dat.ColExt <- left_join(dat.ColExt, butterflies %>% dplyr::select(2,3,5,6,8,10:20,23) %>% group_by(Species, Site, Scale) %>% summarise_all(first))
dat.ColExt <- stdize(dat.ColExt, prefix = F, omit.cols = c("Extinction", "Colonisation", "Scale"))

table(subset(dat.ColExt, dat.ColExt$Scale == 20000)[,c("Extinction")])
table(subset(dat.ColExt, dat.ColExt$Scale == 20000)[,c("Colonisation")])

dat <- subset(dat.ColExt, dat.ColExt$Scale == 30000)

m.col <- glmer(Colonisation ~ STI_rel * PC3 + 
                 STI_rel * PC4 + 
                 STI_rel * PC1 + 
                 PC1 * PLAND * CLUMPY + 
                 Habitat + X*Y + 
                 (1|gridCell50/Site) + (1|Species), 
               family = binomial,
               data = dat,
               nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                             optCtrl=list(maxfun=1e10),
                                             calc.derivs = FALSE), verbose = F, na.action = na.fail)
summary(m.col)

untrace(lme4:::predict.merMod)
visreg2d(m.col, xvar  ="CLUMPY", yvar = "PLAND", scale = "response", cond = list(PC1 = -1))
visreg2d(m.col, xvar  ="CLUMPY", yvar = "PLAND", scale = "response", cond = list(PC1 = 1))


visreg(m.col, xvar  ="PC1", by = "CLUMPY", scale = "response", breaks = 2, rug = F, overlay = T, cond = list(PLAND = 1))



#plot
par(mfrow=c(1,3))
#PC1 x %SNH
newdata1 <- dat %>% summarise_all(median) %>% dplyr::select(-PC1) %>% cbind(., data.frame(PC1 = seq(min(dat$PC1), max(dat$PC1), length.out = 100))) %>% mutate(PLAND = -1)
newdata2 <- dat %>% summarise_all(median) %>% dplyr::select(-PC1) %>% cbind(., data.frame(PC1 = seq(min(dat$PC1), max(dat$PC1), length.out = 100))) %>% mutate(PLAND = 1)

plot(predict(m.col, newdata1, re.form = NA, type = "response") ~ newdata1$PC1, type = "l", col = "blue", ylim = c(0,.05), xlab = "PC1", ylab = "Colonisation probability")
points(predict(m.col, newdata2, re.form = NA, type = "response") ~ newdata1$PC1, type = "l", col = "red")
legend("topright", legend = c("Low %SNH", "High %SNH"), lty = 1, col = c("blue", "red"), bty = "n")


#PC3 x STI
newdata1 <- dat %>% summarise_all(median) %>% dplyr::select(-PC3) %>% cbind(., data.frame(PC3 = seq(min(dat$PC3), max(dat$PC3), length.out = 300))) %>% mutate(STI_rel = -1, PC1 = rep(c(0,-1,1),each = 100))
newdata2 <- dat %>% summarise_all(median) %>% dplyr::select(-PC3) %>% cbind(., data.frame(PC3 = seq(min(dat$PC3), max(dat$PC3), length.out = 300))) %>% mutate(STI_rel = 1, PC1 = rep(c(0,-1,1),each = 100))

plot(predict(m.col, newdata1, re.form = NA, type = "response")[1:100] ~ newdata1$PC3[1:100] , type = "l", col = "blue", ylim = c(0,.05), xlab = "PC3", ylab = "Colonisation probability")
points(predict(m.col, newdata2, re.form = NA, type = "response")[1:100]  ~ newdata1$PC3[1:100] , type = "l", col = "red")
legend("topright", legend = c("Low STI", "High STI"), lty = 1, col = c("blue", "red"), bty = "n")


#PC4 x STI
newdata1 <- dat %>% summarise_all(median) %>% dplyr::select(-PC4) %>% cbind(., data.frame(PC4 = seq(min(dat$PC4), max(dat$PC4), length.out = 300))) %>% mutate(STI_rel = -1, PC1 = rep(c(0,-1,1),each = 100))
newdata2 <- dat %>% summarise_all(median) %>% dplyr::select(-PC4) %>% cbind(., data.frame(PC4 = seq(min(dat$PC3), max(dat$PC4), length.out = 300))) %>% mutate(STI_rel = 1, PC1 = rep(c(0,-1,1),each = 100))

plot(predict(m.col, newdata1, re.form = NA, type = "response")[1:100] ~ newdata1$PC4[1:100] , type = "l", col = "blue", ylim = c(0,.05), xlab = "PC4", ylab = "Colonisation probability")
points(predict(m.col, newdata2, re.form = NA, type = "response")[1:100]  ~ newdata1$PC4[1:100] , type = "l", col = "red")
legend("topright", legend = c("Low STI", "High STI"), lty = 1, col = c("blue", "red"), bty = "n")

##################
### Extinction ###
##################

m.ext <- glmer(Extinction ~ STI_rel * PC3 + 
                 STI_rel * PC4 + 
                 STI_rel * PC1 + 
                 PC1 * PLAND * CLUMPY + 
                 Habitat + X*Y + 
                 (1|gridCell50/Site) + (1|Species), 
               family = binomial,
               data = dat,
               nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                             optCtrl=list(maxfun=1e10),
                                             calc.derivs = FALSE), verbose = F, na.action = na.fail)
summary(m.ext)

par(mfrow=c(2,2))
visreg(m.ext, xvar  ="STI_rel", scale = "response", rug = F)

visreg(m.ext, xvar  ="PC1", by = "CLUMPY", scale = "response", breaks = 2, rug = F, overlay = T)
visreg(m.ext, xvar  ="PC3", by = "STI_rel", scale = "response", breaks = 2, rug = F, overlay = T)
visreg(m.ext, xvar  ="PC4", by = "STI_rel", scale = "response", breaks = 2, rug = F, overlay = T)


###########################
### show classification ###
###########################
dat.temp <- butterflies %>% dplyr::filter(Scale == 20000, Site == unique(butterflies$Site)[[600]])

par(mfrow=c(1,2))

dat.temp2 <- dat.temp %>% dplyr::filter(Species == unique(dat.temp$Species)[[15]])

m <- glm(n ~ Year, family = binomial, dat = dat.temp2)
visreg(m, xvar = "Year", scale = "response", ylim = c(0,1), band = T, rug = F)
points(n ~ Year, dat = dat.temp2)

pred <- predict(m, type = "response", newdata = data.frame(Year = c(min(dat.temp2$Year), max(dat.temp2$Year))))
abline(h = pred[[1]], col = "red")
abline(h = pred[[2]], col = "red")

abline(h = .8, col = "blue", lty = 2)
abline(h = .2, col = "blue", lty = 2)


dat.temp2 <- dat.temp %>% dplyr::filter(Species == unique(dat.temp$Species)[[8]])

m <- glm(n ~ Year, family = binomial, dat = dat.temp2)
visreg(m, xvar = "Year", scale = "response", ylim = c(0,1), band = T, rug = F)
points(n ~ Year, dat = dat.temp2)

pred <- predict(m, type = "response", newdata = data.frame(Year = c(min(dat.temp2$Year), max(dat.temp2$Year))))
abline(h = pred[[1]], col = "red")
abline(h = pred[[2]], col = "red")

abline(h = .8, col = "blue", lty = 2)
abline(h = .2, col = "blue", lty = 2)


###########################
## remove unused objects ##
###########################


gdata:::keep(butterflies, dat.occ.trend,
             sure = T)

