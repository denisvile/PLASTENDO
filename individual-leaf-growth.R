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

gp.ilg %+% subset(ilg, leafType %in% c("F-3")) +
  facet_wrap(.~ leafType+idGenotype) + 
  aes(colour=factor(idPot)) +
  geom_smooth(se=F)


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
fc.sigmoidalFitting <- function(data) { 
  idPot <- unique(data$idPot)
  # Test of the sigmoidal fitting
  fit <- NULL
  for(Btest in 1:10) for (X0test in 2:50) {
    try(fit <- nls(Area.total.mm2  ~ A / ( 1 + exp (-((day - X0) / B))), 
                   data = data, 
                   start = list(A = max(data$Area.total.mm2, na.rm = TRUE), B = Btest , X0 = X0test)),
        silent = T)
    if(!is.null(fit)) break
  }
  # Model parameters
  try(A <- summary(fit)$coef[1], silent = T)
  try(B <- summary(fit)$coef[2], silent = T)
  try(X0 <- summary(fit)$coef[3], silent = T)
  
  # Computation of the maximum expansion rate : ER = A / (4*B)
  try(ER <- 0.25 * A / B, silent = T)
  
  # Computation of the expansion duration : Duration = X0 - B * ln ((1/0.95)-1)
  try(Duration <- X0 - B * log ((1/0.95)-1), silent = T)
  
  if(is.null(fit)) {
    resultTemp <- data.frame(idPot = idPot, A = NA, B = NA, X0 = NA, ER = NA, Duration = NA)
  }
  
  if(!is.null(fit)) {
    resultTemp <- data.frame(idPot = idPot, A = A, B = B, X0 = X0, ER = ER, Duration = Duration)
  }
  
  fun.sigm <- function(x)  A / ( 1 + exp (-((x - X0) / B)))
  
  gp.fit <- ggplot(data = data,
                   aes(x = day,
                       y = Area.total.mm2)) +
    geom_point() +
    ggtitle(idPot) +
    stat_function(fun = fun.sigm) +
    myTheme
  gp.fit
  
  return(list(resultTemp = resultTemp, gp.fit = gp.fit))
}
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

