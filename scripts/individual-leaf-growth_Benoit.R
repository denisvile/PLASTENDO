# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Individual leaf growth data
# Image analyses by Benoît Berthet (and Maïlys Combes)
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€

# Data import ----

ilg1 <- read.xls("./data/Croissance_feuille_C3M42.xlsx")
str(ilg1)

# day = number of days from Feb 1st 2021
unique(ilg1$idGenotype) %in% unique(idPots$nameGen) 
ilg1[ilg1$idGenotype == "kulturen-1", ]$idGenotype <- "Kulturen-1"
ilg1$watering <- factor(ilg1$watering, levels = c("WW", "WD"))
ilg1$leafType <- factor(ilg1$leafType, levels = c("F8", "F-5", "F+2"))



# Arrêt de l'irrigation lorsque F6 était visible
# F8 = une feuille qui a subit le dessèchement du sol (qui a poussé pendant le dessèchement) ou qui n'a pas subit le stress en condition bien irriguée (WW)
# F-5 = la feuilles de rang -5 avant la réirrigation (qui a eu lieu 13 jours après l'atteinte de la teneur en eau du sol cible (70% RWCsoil))
# F-5 en condition WW, c'est une feuille visible au même moment (en moyenne) que la feuille F-5 de la condition WD
# F+2 = la feuilles de rang -5 après la réirrigation (qui a eu lieu 13 jours après l'atteinte de la teneur en eau du sol cible (70% RWCsoil))
# F+2 en condition WW, c'est une feuille visible au même moment (en moyenne) que la feuille F+2 de la condition WD


ilg1 %>%
  group_by(idGenotype, idPot, leafType, watering) %>%
  summarise(mean = mean(Area.total.mm2),
            n = length(Area.total.mm2))

ggplot(subset(ilg1, idGenotype=="An-1"), aes(x=day, y=Area.total.mm2, colour=watering)) + 
  geom_point() + geom_smooth(formula='y~x', aes(group=idPot), se=F) +
  facet_wrap(.~leafType)


pdf("./figures/individual_leaf_growth_Benoit.pdf", 11, 10)
for(i in unique(ilg1$idGenotype)) {
gp.ilg <- ggplot(subset(ilg1, idGenotype==i), aes(x=day, y=Area.total.mm2, colour=as.factor(idPot))) +
  geom_point() + geom_line(aes(linetype=leafType)) +
  facet_wrap(.~idGenotype+leafType+watering, ncol = 2)
print(gp.ilg)
}
dev.off()
system("open ./figures/individual_leaf_growth_Benoit.pdf")

#library(animation)
#ani.record(reset = TRUE)  # clear history before recording


# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Sigmoid fit ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€

# Function for sigmoidal fitting ----
source("./scripts/fc.sigmoidalFitting.r")
# Try with one pot ----
fc.sigmoidalFitting(data=subset(ilg1, idPot==312 & leafType=="F8"), plot=TRUE)

# Plot fitted curves ----
pdf("./figures/growthFitCurves_Benoit_F8.pdf", 7, 6)
for(i in unique(ilg1$idPot)) {
  print(fc.sigmoidalFitting(data=subset(ilg1, idPot==i & leafType=="F8"))$gp.fit)
}
dev.off()
system("open ./figures/growthFitCurves_Benoit_F8.pdf")

# Plot fitted curves ----
pdf("./figures/growthFitCurves_Benoit_F_5.pdf", 7, 6)
for(i in unique(ilg1$idPot)) {
  print(fc.sigmoidalFitting(data=subset(ilg1, idPot==i & leafType=="F-5"))$gp.fit)
}
dev.off()
system("open ./figures/growthFitCurves_Benoit_F_5.pdf")

# Plot fitted curves ----
pdf("./figures/growthFitCurves_Benoit_F+2.pdf", 7, 6)
for(i in unique(ilg1$idPot)) {
  print(fc.sigmoidalFitting(data=subset(ilg1, idPot==i & leafType=="F+2"))$gp.fit)
}
dev.off()
system("open ./figures/growthFitCurves_Benoit_F+2.pdf")


# Empty dataframes ----
fitResults_F8 <- data.frame()
fitResults_F_5 <- data.frame()
fitResults_F_2 <- data.frame()
# Loop on idPot ----
for(i in unique(ilg1$idPot)) {
  fitResults_F8 <- rbind(fitResults_F8,
                      fc.sigmoidalFitting(data=subset(ilg1, idPot==i & leafType=="F8"), plot=FALSE)$resultTemp
  )
}
for(i in unique(ilg1$idPot)) {
  fitResults_F_5 <- rbind(fitResults_F_5,
                         fc.sigmoidalFitting(data=subset(ilg1, idPot==i & leafType=="F-5"), plot=FALSE)$resultTemp
  )
}
for(i in unique(ilg1$idPot)[!(unique(ilg1$idPot) %in% c(230, 366))]) {
  fitResults_F_2 <- rbind(fitResults_F_2,
                          fc.sigmoidalFitting(data=subset(ilg1, idPot==i & leafType=="F+2"), plot=FALSE)$resultTemp
  )
}

fitResults_F8$leafType <- "F8"
fitResults_F_5$leafType <- "F-5"
fitResults_F_2$leafType <- "F+2"

fitResults.pots <- rbind(fitResults_F8, fitResults_F_5, fitResults_F_2)



# Plot ERmax ----

fitResults.pots <- fitResults.pots %>%
  left_join(idPots)
fitResults.pots$watering <- factor(fitResults.pots$watering, levels = c("WW", "WD"))

# idPot = 110 removed because of pot number mistake
pdf("./figures/ER_F8_genotypes.pdf", 10, 6)
ggplot(subset(fitResults.pots, leafType=="F8" & !(idPot %in% c(110))), aes(y=ER, x=interaction(watering, nameGen), fill=watering)) +
  geom_boxplot() +
  ylab(expression(paste("Leaf 8 max. expansion rate (", mm^2," ", d^-1, ")"))) +
  xlab("") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12, colour="black", angle=90, hjust=1),
        legend.position = "none")
dev.off()
system("open ./figures/ER_F8_genotypes.pdf")

fitResults.pots.m <- fitResults.pots %>%
  group_by(nameGen, watering, leafType) %>%
  summarise(ER = mean(ER, na.rm=T))

fitResults.pots.m <- fitResults.pots.m %>%
  left_join(CV.mean.all.3datasets)

ggplot(data=subset(fitResults.pots.m, leafType=="F8" & !(nameGen %in% "Rennes-1") & watering=="WW"), aes(x=Leaf_8_WW, y=ER, colour=watering)) +
  geom_point() + geom_smooth(method="lm", se=F) +
  facet_wrap(.~leafType) +
  theme_bw()

ggplot(data=subset(fitResults.pots.m, leafType=="F-5" & !(nameGen %in% "Rennes-1") & watering=="WW"), aes(x=Leaf_30_WW, y=ER, colour=clust)) +
  geom_point(size=4) + 
  facet_wrap(.~leafType) +
  theme_bw()

fitResults.pots.m$nameGen



