# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Individual leaf growth data
# Image analyses by Benoît Berthet (and Maïlys Combes)
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€

# Data import ----

ilg1 <- read.xls("./data/Croissance_feuille_v4.xlsx")
str(ilg1)

# day = number of days from Feb 1st 2021
unique(ilg1$idGenotype) %in% unique(idPots$nameGen) 
# ilg1[ilg1$idGenotype == "kulturen-1", ]$idGenotype <- "Kulturen-1"
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
  geom_point() + 
  geom_line(aes(linetype=leafType)) +
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

pdf("./figures/ER_F_5_genotypes.pdf", 10, 6)
ggplot(subset(fitResults.pots, leafType=="F-5" & !(idPot %in% c(110))), aes(y=ER, x=interaction(watering, nameGen), fill=watering)) +
  geom_boxplot() +
  ylab(expression(paste("Leaf -5 max. expansion rate (", mm^2," ", d^-1, ")"))) +
  xlab("") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12, colour="black", angle=90, hjust=1),
        legend.position = "none")
dev.off()
system("open ./figures/ER_F_5_genotypes.pdf")

d.SLA.m.0$leafType <- factor(d.SLA.m.0$stage, labels = c("F+2", "F8"))

fitResults.pots.m <- fitResults.pots %>%
  group_by(nameGen, watering, leafType) %>%
  summarise(ER = mean(ER, na.rm=T),
            A = mean(A, na.rm=T),
            B = mean(B, na.rm=T),
            X0 = mean(X0, na.rm=T),
            Duration = mean(Duration, na.rm=T)) %>%
  left_join(d.SLA.m.0) %>%
  left_join(CV.mean.all.3datasets)

# Plot sigmoidal growth curves ----

fun.sigm <- function(x, A, B, X0) {
  y = A / ( 1 + exp (-((x - X0) / B)))
  return(y) 
  }

fitResults.pots.m

ggplot(data.frame(x = c(0, 30)), aes(x)) + ylab("Leaf area (mm2)") +
  mapply(function(A, B, X0, nameGen, watering) {
    stat_function(fun = fun.sigm, args = list(A = A, B = B, X0=X0), aes(colour=nameGen, linetype=watering)) 
  }, 
  # enter A, B, and colors here
  A = subset(fitResults.pots.m, leafType=="F8")$A, 
  B = subset(fitResults.pots.m, leafType=="F8")$B,
  X0 = subset(fitResults.pots.m, leafType=="F8")$X0,
  nameGen = subset(fitResults.pots.m, leafType=="F8")$nameGen,
  watering = subset(fitResults.pots.m, leafType=="F8")$watering)


fitResults.pots.m.wide <- dcast(fitResults.pots.m, formula=nameGen + leafType ~ watering, mean, value.var="ER")

fitResults.pots.m.wide <- fitResults.pots.m.wide %>%
  dplyr::filter(nameGen != "Rennes-1") %>%
  mutate(ERrr = WD/WW) %>%
  left_join(CV.mean.all.3datasets)

fitResults.pots.m.wide$nameGen.ord <- factor(fitResults.pots.m.wide$nameGen, levels = subset(fitResults.pots.m.wide, leafType=="F8")$nameGen[order(subset(fitResults.pots.m.wide, leafType=="F8")$ERrr, decreasing=T)])

ggplot(data=subset(fitResults.pots.m, leafType=="F8"), aes(y=leaf.area.mm2, x=LAmax, colour=watering)) +
  geom_point() + geom_smooth(method="lm", se=F) + geom_abline(slope=1, intercept=0) +
  ylab("Leaf area of L9 (scanned individual leaf)") + xlab("Leaf area of L8 (measured on rosettes)") +
  facet_wrap(.~leafType) +
  theme_bw()

ggplot(data=subset(fitResults.pots.m, leafType %in% c("F8", "F+2")), aes(y=ER, x=SLA, colour=watering)) +
  geom_point() + geom_smooth(method="lm", se=F) + 
  facet_wrap(.~leafType)

ggplot(data=subset(fitResults.pots.m, leafType %in% c("F8")), aes(y=Leaf_8_WW, x=ER, colour=watering)) +
  geom_point() + geom_smooth(method="lm", se=F) + 
  scale_x_log10() +
  facet_wrap(.~leafType)


ggplot(data=subset(fitResults.pots.m, leafType=="F8" & !(nameGen %in% "Rennes-1") & watering=="WW"), aes(x=Leaf_8_WW, y=ER, colour=watering)) +
  geom_point() + geom_smooth(method="lm", se=F) +
  facet_wrap(.~leafType) +
  theme_bw()

ggplot(data=subset(fitResults.pots.m, leafType=="F+2" & !(nameGen %in% "Rennes-1") & watering=="WW"), aes(x=Seedling_Leaf5_WW, y=ER, colour=clust)) +
  geom_point(size=4) + 
  facet_wrap(.~leafType) +
  theme_bw()

ggplot(fitResults.pots.m.wide, aes(x=Seedling_Leaf5_WW, y=ERrr, colour=leafType)) +
  geom_point() + facet_wrap(.~leafType)
ggplot(fitResults.pots.m.wide, aes(x=Leaf_8_WW, y=ERrr, colour=leafType)) +
  geom_point() + facet_wrap(.~leafType)

ggplot(subset(fitResults.pots.m.wide, leafType=="F8" & nameGen != "Rennes-1"), aes(x=nameGen.ord, y=ERrr)) + geom_bar(stat = "identity")

gp.RR <- ggplot(subset(fitResults.pots.m.wide, leafType=="F-5"), aes(x=Leaf_8_WD/Leaf_8_WW, y=ERrr)) +
  geom_point() +
  geom_smooth(method=lm, formula = y~x, se=F) +
  geom_text_repel(aes(label=nameGen)) +
 # facet_wrap(.~leafType) +
  myTheme

pdf("./figures/change_LER_EF.pdf", 6, 6)
gp.RR + xlab("Change in leaf 8 EF") + ylab("Change in L30 ER")
gp.RR %+% subset(fitResults.pots.m.wide, leafType=="F8") + xlab("Change in L8 EF") + ylab("Change in L8  expansion rate")
gp.RR + aes(x=Leaf_30_WD/Leaf_30_WW, y=ERrr) + xlab("Change in leaf 30 EF") +
  ylab("Change in L30 expansion rate")
dev.off()
system("open ./figures/change_LER_EF.pdf")

gp.RR %+% subset(fitResults.pots.m.wide, leafType=="F8") + aes(y=ERrr, x=Seedling_Leaf5_WW)

ggplot(subset(fitResults.pots.m.wide, leafType=="F-5"), aes(x=Leaf_8_WD/Leaf_8_WW, y=ERrr)) +
  geom_text_repel(aes(label=nameGen)) +
  geom_point() + facet_wrap(.~leafType) +
  geom_smooth(method=lm, formula = y~x, se=F)

ggplot(subset(fitResults.pots.m.wide, leafType=="F+2"), aes(x=Leaf_8_WD/Leaf_8_WW, y=WD/WW, colour=clust)) + geom_text_repel(aes(label=nameGen)) +
  geom_point() + facet_wrap(.~leafType)

