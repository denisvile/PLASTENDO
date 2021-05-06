# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Individual leaf growth data
# Image analyses by M1 Plant Science students (Marie, Mouad, Harkingto)
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€

library(dplyr)
library(ggplot2)
library(gdata)


load("rdata.RData")
# Data import ----

ilg <- read.xls("./data/C3M42_leafGrowth.xlsx")
str(ilg)
# day
# day 1 = 2021-02-01

ilg$idGenotype <- as.factor(ilg$idGenotype)
ilg$leafType <- factor(ilg$leafType, levels=c("F8", "F-3", "F+2"))

gp.ilg <- ggplot(data = ilg, aes(y = Area.total.mm2, x = day)) +
  geom_point()
gp.ilg

pdf("./figures/L8_growth_8genotypes.pdf", 6, 8)
gp.ilg %+% subset(ilg, leafType %in% c("F8")) +
  facet_wrap(.~ idGenotype, ncol=2) + 
  aes(colour=factor(idPot)) +
  geom_smooth(se=F) +
  theme(legend.position = "none")
dev.off()
system("open ./figures/L8_growth_8genotypes.pdf")


gp.ilg %+% subset(ilg, idGenotype %in% c("Cdm 0", "Com-1")) +
  facet_wrap(idGenotype ~ leafType) + 
  aes(colour=factor(idPot)) +
  geom_smooth(se=F)

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Sigmoid fit for leaf 8 (F8) under water deficit ----
# (decreasing soil water content following watering stop)
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Function for sigmoidal fitting ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€

source("./scripts/fc.sigmoidalFitting.r")

# Empty dataframe ----
fitResults <- data.frame()

# Try with one pot ----
fc.sigmoidalFitting(data=subset(ilg, idPot==270 & leafType=="F8"))

# Loop on idPot ----
for(i in unique(ilg$idPot)) {
  fitResults <- rbind(fitResults,
                      fc.sigmoidalFitting(data=subset(ilg, idPot==i & leafType=="F8"))$resultTemp
  )
}

# Plot fitted curves ----
pdf("./figures/growthFitCurves.pdf", 7, 6)
for(i in unique(ilg$idPot)) {
  print(fc.sigmoidalFitting(data=subset(ilg, idPot==i & leafType=="F8"))$gp.fit)
}
dev.off()
system("open ./figures/growthFitCurves.pdf")


# Plot expansion rate (ER) ----

fitResults.Pots <- ilg %>%
  group_by(idPot) %>%
  select(idGenotype, idPot, watering) %>%
  dplyr::filter(row_number() == 1) %>%
  left_join(fitResults)
fitResults.Pots$idGenotype <- factor(fitResults.Pots$idGenotype, levels = c(
  "Com-1", "Ei 02","IP-Coa","TBO 1", "Cvi-0", "Sanna-2","Vinslov","Cdm 0" ))

pdf("./figures/ER_F8_8genotypes.pdf", 8, 6)
ggplot(subset(fitResults.Pots, !(idPot %in% c(21, 316, 213))), aes(y=ER, x=idGenotype, fill=idGenotype)) +
  geom_boxplot() +
  ylab(expression(paste("Leaf 8 expansion rate (", mm^2," ", d^-1, ")"))) +
  xlab("") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"),
        legend.position = "none")
dev.off()
system("open ./figures/ER_F8_8genotypes.pdf")

m1 <- aov(ER ~ idGenotype, data=fitResults.Pots)
summary(m1)
TukeyHSD(m1)

