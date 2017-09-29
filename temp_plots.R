plot.trends.but.Ab <- as.data.frame(cld(lsmeans_but_Ab, Letters =  letters))
plot.trends.but.P <- as.data.frame(cld(lsmeans_but_P, Letters =  letters))

plot.trends.but <- rbind.data.frame(cbind.data.frame(type = "Abundance", plot.trends.but.Ab),
                                    cbind.data.frame(type = "Presence-only", plot.trends.but.P))


p1 <- ggplot(data = plot.trends.but, aes(y = Year.trend, x = as.factor(PLAND), fill = as.factor(CLUMPY))) + 
  facet_grid( type ~ ., scales = "free") + 
  geom_bar(stat = "identity", position=position_dodge(.9)) +
  geom_errorbar(aes(ymin = Year.trend - SE, ymax = Year.trend + SE), width = 0.2, position=position_dodge(.9)) +
  geom_text(aes(y = (Year.trend + SE) + .01, 
                x = as.factor(PLAND),
                label = .group), position=position_dodge(.9)) +
  scale_fill_discrete	("Fragmentation", 
                       labels = c("High fragmentation", 
                                  "Low fragmentation")) +
  scale_x_discrete("% SNH area",
                   labels = c("Low area of SNH", 
                              "High area of SNH")) +
  scale_y_continuous("CTI trend") + 
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 


###
plot.trends.bird.Ab <- as.data.frame(cld(lsmeans_bird_Ab, Letters =  letters))
plot.trends.bird.P <- as.data.frame(cld(lsmeans_bird_P, Letters =  letters))

plot.trends.bird <- rbind.data.frame(cbind.data.frame(type = "Abundance", plot.trends.bird.Ab),
                                     cbind.data.frame(type = "Presence-only", plot.trends.bird.P))


p2 <- ggplot(data = plot.trends.bird, aes(y = Year.trend, x = as.factor(PLAND), fill = as.factor(CLUMPY))) + 
  facet_grid( type ~ ., scales = "free") + 
  geom_bar(stat = "identity", position=position_dodge(.9)) +
  geom_errorbar(aes(ymin = Year.trend - SE, ymax = Year.trend + SE), width = 0.2, position=position_dodge(.9)) +
  geom_text(aes(y = (Year.trend + SE) + .01, 
                x = as.factor(PLAND),
                label = .group), position=position_dodge(.9)) +
  scale_fill_discrete	("Fragmentation", 
                       labels = c("High fragmentation", 
                                  "Low fragmentation")) +
  scale_x_discrete("% SNH area",
                   labels = c("Low area of SNH", 
                              "High area of SNH")) +
  scale_y_continuous("CTI trend") + 
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 


plot_grid(p1, p2, nrow = 2, labels = c("Butterflies", "Birds"), align = "hv", hjust = 0, vjust = 1.5)



#####
m_frag_butterflies_Ab1 <- lmer(cti ~ CLUMPY * PLAND * as.factor(Year) + X*Y +
                                 (1|country/gridCell50/Site),
                               data = data.frame(std.butterflies.data.scale_Ab %>% filter(!is.na(CLUMPY))), REML = F)

a <- lsmeans(m_frag_butterflies_Ab1,
             ~ Year | CLUMPY * PLAND,
             at= list(CLUMPY = c(-1,1), PLAND = c(-1,1)))
a1 <- summary(a)

ggplot(data = a1, aes(y = lsmean, x = Year, color = as.factor(CLUMPY))) + 
  facet_grid( ~ PLAND, scales = "free") + 
  geom_point() + geom_line()
