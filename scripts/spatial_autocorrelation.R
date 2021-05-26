# Spatial autocorrelation

d.SLA$xy <- as.numeric(d.SLA$idPot)

d.x <- NULL
for(i in 1:42)  {
  d.xi = rep(i, 12)
  d.x <- c(d.x, d.xi)
}
d.xy <- data.frame(xy=1:504, y=rep(12:1, 42), x=d.x)

pheno.design <- merge(d.SLA, d.xy, by="xy")

pheno.design %>%
  group_by(nameGen, watering) %>%
  dplyr::summarise(n = length(idPot))

ggplot(pheno.design, aes(x=x, y=y, colour=fresh_weight_mg)) + 
  geom_point(size=4)

# Model testing with LME ----
# Procedure from http://eco-stats.blogspot.com/2014/10/r-lab-inference-with-spatially.html
#We will need a dummy variable because lme requires a grouping variable, but it will be the same for all observations, hence dummy

pheno.design$dummy <- 1
#First we fit a model with no spatial correlation

nonsp.null = lme(log10(fresh_weight_mg) ~ 1, data = pheno.design, random = ~ 1 | dummy, method='ML', na.action=na.omit)
nonsp.wat = lme(log10(fresh_weight_mg) ~ watering, data = pheno.design, random = ~ 1 | dummy, method='ML', na.action=na.omit)

anova(nonsp.null,nonsp.wat) # Highly significant effect of watering


vg.null <- Variogram(nonsp.null,resType="normalized", form=~x+y)
gg.vg <- ggplot(data=vg.null, aes(y=variog, x=dist)) + geom_point(size=4) + theme_bw() + geom_smooth(se=F) +
  xlab("Distance")+ylab("Variogram")
gg.vg
vg.wat <- Variogram(nonsp.wat, resType="normalized", form=~x+y) # corGaus may provide a good fit
gg.vg %+% vg.wat


# Next we will fit some correlations functions to see what fits best
# We want the variogram to be flat when we include a spatial random effect with the correct structure.
gauss.null <- lme(log10(fresh_weight_mg) ~ 1, data = pheno.design, random = ~ 1 | dummy,
                  correlation = corGaus(form = ~ x+y| dummy),method='ML', na.action=na.omit)
vg.gauss <- Variogram(gauss.null, resType="normalized")
gg.vg %+% vg.gauss

ratio.null <- lme(log10(fresh_weight_mg) ~ 1, data = pheno.design, random = ~ 1 | dummy,
                  correlation = corRatio(form = ~ x+y| dummy),method='ML', na.action=na.omit)
vg.ratio <- Variogram(ratio.null,resType="normalized")
gg.vg %+% vg.ratio

anova(nonsp.null, gauss.null, ratio.null)

spher.null <- lme(log10(fresh_weight_mg) ~ 1, data = pheno.design, random = ~ 1 | dummy,
                  correlation = corSpher(form = ~ x+y| dummy),method='ML', na.action=na.omit)
plot(Variogram(spher.null,resType="normalized"))

anova(nonsp.null, gauss.null, spher.null) # corGauss provide a similar fit than sperical correlation structure.

# The Gaussian correlation structure seems to perform the best, so we'll use that.

gauss.wat <- lme(log10(fresh_weight_mg) ~ watering, data = pheno.design, random = ~ 1 | dummy,
                 correlation = corGaus(form = ~ x+y| dummy), method='ML', na.action=na.omit)
vg.gauss.wat <- Variogram(gauss.wat,resType="normalized")
gg.vg %+% vg.gauss.wat
anova(nonsp.wat, gauss.null, gauss.wat)
anova(nonsp.wat, gauss.wat)

ratio.wat <- lme(log10(fresh_weight_mg) ~ watering, data = pheno.design, random = ~ 1 | dummy,
                 correlation = corRatio(form = ~ x+y| dummy), method='ML', na.action=na.omit)
vg.ratio.wat <- Variogram(ratio.wat,resType="normalized")
gg.vg %+% vg.ratio.wat
anova(nonsp.wat, gauss.null, gauss.wat, ratio.wat)
anova(nonsp.wat, gauss.wat)
anova(nonsp.wat, ratio.wat)

gp.watering <- ggplot(data=pheno.design, aes(x=watering, y=fresh_weight_mg, fill=watering)) + 
  geom_boxplot() +
  # scale_y_log10() + 
  ylab("Rosette fresh weight (mg)") +
  scale_x_discrete("Watering") +
  theme_classic() + theme(legend.position="none",
                          axis.text = element_text(size=14, colour="black"),
                          axis.title = element_text(size=14))
pdf("./figures/Watering_FM.pdf", 5, 7)
gp.watering
dev.off()
system("open ./figures/Watering_FM.pdf")
# Save PDF ----
pdf("./figures/Variogram.pdf", 10, 5)
gg.vg + geom_point(data=vg.gauss, colour="red", size=4) + geom_smooth(data=vg.gauss, colour="red", se=F) 
gg.vg %+% vg.wat + geom_point(data=vg.gauss.wat, colour="red", size=4) + geom_smooth(data=vg.gauss.wat, colour="red", se=F) +
  geom_point(data=vg.ratio.wat, colour="green", size=4) + geom_smooth(data=vg.ratio.wat, colour="green", se=F) 
dev.off()
system("open ./figures/Variogram.pdf")


# Adding a Grouping variable (idGen may be tested here) ----
#Sometimes we may only want to model spatial correlation for observations within a value of some grouping variable. We will used species as a grouping variable, but this does not take into account correlation between species

nonsp.Gen.null = lme(log10(fresh_weight_mg) ~ 1, data = pheno.design, random = ~ 1 | nameGen, method='ML', na.action=na.omit)
nonsp.Gen.wat = lme(log10(fresh_weight_mg) ~ watering, data = pheno.design, random = ~ 1 | nameGen, method='ML', na.action=na.omit)

group.null <- lme(fresh_weight_mg ~ 1, data = pheno.design, random = ~ 1 | nameGen,
                  correlation = corGaus(form = ~ x+y| nameGen),method='ML', na.action=na.omit)
group.watering <- lme(fresh_weight_mg ~ watering, data = pheno.design, 
                 random = ~ 1 | nameGen,
                 correlation = corGaus(form = ~x+y| nameGen),method='ML', na.action=na.omit)
anova(nonsp.Gen.null, nonsp.Gen.wat)

anova(group.null, group.watering)



pheno.m <- pheno.design %>%
  group_by(nameGen, watering) %>%
  summarise(cv = cv((fresh_weight_mg), na.rm=T),
            n = length(na.omit(fresh_weight_mg))) %>%
  arrange(-cv)

ggplot(data=pheno.m, aes(y = cv, x = watering, fill=watering)) + geom_boxplot()
pheno.m %>%
  group_by(watering) %>%
  summarise(mean = mean(cv, na.rm=T),
            sd = sd(cv, na.rm=T))






