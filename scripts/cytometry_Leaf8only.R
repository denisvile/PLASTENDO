# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Cytometry of Arabidopsis leaf samples
# Project: [ARABREED]
# Subproject: plasticity of endopolyploidy [PLASTENDO] 2021
# Benoît Berthet Master 2 
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€

# source("cytometry.libraries.r")

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Pots ids ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
idPots <- read.xls(xls = "./data/pot_C3M42.xlsx")


# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Sample ids ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
ids <- read.xls(xls = "/Users/denis/Documents/Encadrements/Stages/2021\ -\ M2\ -\ Benoit\ Berthet\ -\ Endopolyploidy/Experiment/endopolyploidy/cytometry_sample_IDs.xlsx")


# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Flow cytometry data ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
d <- read.flowSet(path=paste("/Users/denis/Documents/Encadrements/Stages/2021\ -\ M2\ -\ Benoit\ Berthet\ -\ Endopolyploidy/Experiment/endopolyploidy/data.Endopolyploidy.Nom_fichier_corrige_Leaf8only"))
d.flowset <- d

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Flow cytometry data transformation ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
dat_L8 <- transform(d, "lgDAPI"=log10(`DAPI`))

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Automated filtering of flow cytometry data by curve1Filter ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Global ----
resL8 <- filter(dat_L8[,], curv1Filter("lgDAPI", bwFac=0.48)) 
resSumL8 <- summary(resL8)

dfResL8 <- toTable(resSumL8)
names(dfResL8); head(dfResL8); str(dfResL8)
dfResL8$population <- as.factor(dfResL8$population)
dfResL8$population <- factor(dfResL8$population, labels=c("peak1", "peak2", "peak3", "peak4", "peak5","peak6", "peak7", "peak8", "rest"))
dfResL8.wide <- reshape(dfResL8[, c("p", "sample", "population")], v.names = "p",  idvar = "sample", timevar = "population", direction = "wide")
names(dfResL8.wide)
dfResL8.wide <- do.call(data.frame, lapply(dfResL8.wide, function(x) replace(x, is.na(x), 0)))

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Calculate cycle value (endoreplication factor) ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
dfCV_L8 <- within(data = dfResL8.wide, { cycleValue <- (0*p.peak1 + 1*p.peak2 + 2*p.peak3 + 3*p.peak4 + 4*p.peak5 + 5*p.peak6 + 6*p.peak6 + 7*p.peak8)/(p.peak1 + p.peak2 + p.peak3 + p.peak4 + p.peak5 + p.peak6 + p.peak7 + p.peak8) })

ggplot(subset(dfResL8), aes(y=percent, x=population)) + geom_boxplot()

ggplot(subset(dfCV_L8), aes(x=cycleValue)) + geom_histogram()

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Split frames per filtering ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# split(dat, res)

split.resL8 <- split(dat_L8, resL8, population=list(keep=c("peak 1")))
Lth <- as.numeric(summary(split.resL8L8)[1])

df.resL8 <- data.frame(min1=rep(NA, Lth), max1=rep(NA, Lth), 
                     min2=rep(NA, Lth), max2=rep(NA, Lth),
                     min3=rep(NA, Lth), max3=rep(NA, Lth),
                     min4=rep(NA, Lth), max4=rep(NA, Lth),
                     min5=rep(NA, Lth), max5=rep(NA, Lth),
                     min6=rep(NA, Lth), max6=rep(NA, Lth),
                     min7=rep(NA, Lth), max7=rep(NA, Lth)
                     )

for(i in 1:Lth){
  df.resL8[i, c("min1", "max1")] <- 10^(range(exprs(split.resL8$keep[[i]]$"lgDAPI")))
}
split.resL8 <- split(dat_L8, resL8, population=list(keep=c("peak 2")))
for(i in 1:Lth){
  df.resL8[i, c("min2", "max2")] <- 10^(range(exprs(split.resL8$keep[[i]]$"lgDAPI")))
}
split.resL8 <- split(dat_L8, resL8, population=list(keep=c("peak 3")))
for(i in 1:Lth){
  df.resL8[i, c("min3", "max3")] <- 10^(range(exprs(split.resL8$keep[[i]]$"lgDAPI")))
}
split.resL8 <- split(dat_L8, resL8, population=list(keep=c("peak 4")))
for(i in 1:Lth){
  df.resL8[i, c("min4", "max4")] <- 10^(range(exprs(split.resL8$keep[[i]]$"lgDAPI")))
}
split.resL8 <- split(dat_L8, resL8, population=list(keep=c("peak 5")))
for(i in 1:Lth){
  df.resL8[i, c("min5", "max5")] <- 10^(range(exprs(split.resL8$keep[[i]]$"lgDAPI")))
}
split.resL8 <- split(dat_L8, resL8, population=list(keep=c("peak 6")))
for(i in 1:Lth){
  df.resL8[i, c("min6", "max6")] <- 10^(range(exprs(split.resL8$keep[[i]]$"lgDAPI")))
}
split.resL8 <- split(dat_L8, resL8, population=list(keep=c("peak 7")))
for(i in 1:Lth){
  df.resL8[i, c("min7", "max7")] <- 10^(range(exprs(split.resL8$keep[[i]]$"lgDAPI")))
}
split.resL8 <- split(dat_L8, resL8, population=list(keep=c("peak 8")))
for(i in 1:Lth){
  df.resL8[i, c("min8", "max8")] <- 10^(range(exprs(split.resL8$keep[[i]]$"lgDAPI")))
}


# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
#  ave df.res into d1.L8 for safe handling ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€

d1.L8 <- df.resL8

d1.L8 <- do.call(data.frame,lapply(df.res, function(x) replace(x, is.infinite(x),NA)))
d1.L8 <- do.call(data.frame,lapply(d1.L8, function(x) replace(x, x==0,NA)))

d1.L8$names <- pData(d[1:622,])
d1.L8$num <- 1:Lth

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Plot fluorescence peaks automatically detected by curve1filter() ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
pdf("./Figures/Fluorescence_peaks_L8.pdf", 7, 7)
gg1 <- ggplot(d1.L8, aes(x=min1, y=num)) + scale_x_log10() + 
  geom_segment(aes(xend=max1, ystart=num, yend=num)) + 
  geom_segment(aes(x=min2, xend=max2, ystart=num, yend=num), col="red") +
  geom_segment(aes(x=min3, xend=max3, ystart=num, yend=num), col="blue") +
  geom_segment(aes(x=min4, xend=max4, ystart=num, yend=num), col="green") +
  geom_segment(aes(x=min5, xend=max5, ystart=num, yend=num), col="violet") +
  geom_segment(aes(x=min6, xend=max6, ystart=num, yend=num), col="grey") +
  geom_segment(aes(x=min7, xend=max7, ystart=num, yend=num), col="yellow") +
  geom_segment(aes(x=min8, xend=max8, ystart=num, yend=num), col="orange") +
  xlab("Fluorescence") + ylab("Sample") + theme_bw()
gg1
dev.off()
system(paste("open", "./Figures/Fluorescence_peaks_L8.pdf"))

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 1/ calculate mean limits of each peak / day of measurement ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
df.m.L8 <- melt(d1.L8[, c("min1", "max1", "min2", "max2", "min3", "max3", "min4", "max4", "min5", "max5", "min6", "max6", "min7", "max7", "min8", "max8")], 
             measure.vars = c("min1", "max1", "min2", "max2", "min3", "max3", "min4", "max4", "min5", "max5", "min6", "max6", "min7", "max7", "min8", "max8"))
df.m.L8$mM <- factor(df.m.L8$variable, labels=rep(c("min","max"), 8))

gp2.L8 <- ggplot(df.m.L8, aes(y=value, x=variable, fill=mM)) + geom_boxplot() + scale_y_log10(limits=c(1, 1000)) + ylab("") + xlab("Peak limits") + coord_flip()

peak.limits.L8 <- aggregate(value ~ variable, data=df.m.L8, FUN = summary)

# Plot limits 
gp2.L8 + geom_hline(data=data.frame(yint=peak.limits.L8[, 2][, "1st Qu."]), aes(yintercept=yint), col="gray") + 
  geom_hline(data=data.frame(yint=peak.limits.L8[, 2][, "3rd Qu."]), aes(yintercept=yint), col="gray") + 
  theme_bw()

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 1.a/ Define the limits using 1st and 3rd quartiles ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
pl.L8 <- peak.limits.L8[, 2]
peak.lim.L8 <- list(
  peak1=c(pl.L8[1, "Min."], mean(c(pl.L8[2, "3rd Qu."], pl.L8[3,"1st Qu."]))),
  peak2=c(mean(c(pl.L8[2, "3rd Qu."], pl.L8[3,"1st Qu."])), 
          mean(c(pl.L8[4, "3rd Qu."], pl.L8[5,"1st Qu."]))),
  peak3=c(mean(c(pl.L8[4, "3rd Qu."], pl.L8[5,"1st Qu."])), 
          mean(c(pl.L8[6, "3rd Qu."], pl.L8[7,"1st Qu."]))),
  peak4=c(mean(c(pl.L8[6, "3rd Qu."], pl.L8[7,"1st Qu."])),
          mean(c(pl.L8[8, "3rd Qu."], pl.L8[9,"1st Qu."]))),
  peak5=c(mean(c(pl.L8[8, "3rd Qu."], pl.L8[9,"1st Qu."])),
          mean(c(pl.L8[10, "3rd Qu."], pl.L8[11,"1st Qu."]))),
  peak6=c(mean(c(pl.L8[10, "3rd Qu."], pl.L8[11,"1st Qu."])),
          mean(c(pl.L8[12, "3rd Qu."], pl.L8[13,"1st Qu."]))),
  peak7=c(mean(c(pl.L8[12, "3rd Qu."], pl.L8[13,"1st Qu."])), pl.L8[15, "Max."]),
  peak8=c(pl.L8[15, "Max."], pl.L8[16, "Max."])
)
pl1.L8 <- melt(as.data.frame(peak.lim.L8))
gp2.L8 + geom_hline(data=pl1.L8, aes(yintercept=value), col="blue") + theme_bw()

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 1.b/ Plot global peak limits ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
pdf("./figures/estimatedPeakLimits_boxplots.L8.pdf", 7, 7)
gp2.L8 + geom_hline(data=pl1.L8, aes(yintercept=value), col="blue") + theme_bw() 
dev.off()
system("open ./figures/estimatedPeakLimits_boxplots.L8.pdf")

system("open ../figures/detected_and_estimatedPeakLimits.pdf")
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 2/ use these limits as filters in filter() ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
rectGate1.L8 <- rectangleGate(filterId = "peak1", "DAPI"=peak.lim.L8$peak1)
#rectGate1.1 <- rectangleGate(filterId = "peak1", "FL1-A"=peak.lim1.manual$peak1)
#rectGate1.2 <- rectangleGate(filterId = "peak1", "FL1-A"=peak.lim2$peak1)

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 3/ use these limits to extract peak % ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 
fres00.L8 <- filter(dat_L8, filter=rectGate1.L8)
cycleProportion.L8 <- toTable(summary(fres00.L8))
for(i in c("peak2", "peak3", "peak4", "peak5", "peak6", "peak7", "peak8")) {
  rectGate.i <- rectangleGate(filterId = i, "DAPI"=peak.lim.L8[[i]])
  fres00.L8 <- filter(dat_L8, filter=rectGate.i)
  cycleProportion.L8 <- rbind(cycleProportion.L8, toTable(summary(fres00.L8)))
}

cycleProportion.L8$ploidy <- factor(cycleProportion.L8$population, labels=c("debris", "x2C", "x4C", "x8C", "x16C", "x32C", "x64C", "x128C"))
cycleProportion.L8.wide <- dcast(cycleProportion.L8[, c("sample", "ploidy", "p")], formula=sample~ploidy)
cycleProportion.L8.wide[is.na(cycleProportion.L8.wide$x128C), "x128C"] <- 0
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Calculate new cycle value (endoreduplication factor) ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# The two periods of measurement treated conjointly ----
dfCV.L8 <- within(data = cycleProportion.L8.wide, {
  cycleValue <- (0*x2C + 1*x4C + 2*x8C + 3*x16C + 4*x32C + 5*x64C + 6*x128C)/(x2C + x4C + x8C + x16C + x32C + x64C + x128C)
})
dfCV.L8$cycleValue
dfCV.L8 <- dfCV.L8 %>%
  left_join(ids, by=c("sample"="fileName")) %>%
  left_join(idPots)

dfCV.L8$tissueType.ord <- factor(dfCV.L8$tissueType, 
                             levels = c("f6", "f8", "F8", "f30"), 
                             labels = c("Seedling_Leaf5","Leaf_8", "Leaf_8", "Leaf_30"))

dfCV.L8[dfCV.L8$tissueType.ord=="Seedling@Leaf5", "watering"] <- "WW"
dfCV.L8$watering <- factor(dfCV.L8$watering, levels = c("WW", "WD"))

CV.mean.L8 <- subset(dfCV.L8) %>%
  select(nameGen, tissueType.ord, watering, cycleValue) %>%
  group_by(nameGen, tissueType.ord, watering) %>%
  summarize(CVmean = mean(cycleValue))

CV.mean.L8.wide <- dcast(CV.mean.L8, formula=nameGen~tissueType.ord + watering, mean)

dfCV.L8$nameGen.OrderedL8_WW <- as.factor(dfCV.L8$nameGen)
dfCV.L8$nameGen.OrderedL8_WW <- factor(dfCV.L8$nameGen.OrderedL8_WW,
                                          levels =levels(dfCV.L8$nameGen.OrderedL8_WW)[order(CV.mean.L8.wide$Leaf_8_WW)])

pdf("./figures/cycleValue_30_genotypes.L8.pdf", 12, 8)
ggplot(data=dfCV.L8[1:622,], aes(y=cycleValue, x=tissueType.ord, colour=watering)) + geom_boxplot() + 
  facet_wrap(.~nameGen)

ggplot(data=subset(dfCV.L8, tissueType.ord%in%c("Leaf_8")), 
       aes(y=cycleValue, x=nameGen.OrderedL8_WW, fill=watering)) +
  geom_boxplot(outlier.alpha = 0) +
  xlab("Accessions (ordered by L8_WW CV)") +
  ylab("Cycle value") +
  scale_fill_manual(values=c("#386FA4", "#84D2F6")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))

dev.off()
system("open ./figures/cycleValue_30_genotypes.L8.pdf")

gp.corr <- ggplot(data=CV.mean.wide, aes(x = Leaf_8_WW, y = Leaf_8_WD)) +
  geom_point() + geom_smooth(method = lm, se=F)

pdf("./figures/cycleValue_30_genotypes_correlations.L8_vs_all.pdf", 8, 7)
gp.corr + geom_abline(slope = 1, intercept = 0)
dev.off()
system("open ./figures/cycleValue_30_genotypes_correlations.L8_vs_all.pdf")

dfCV.L8_all <- dfCV.L8[, c("nameGen", "watering", "tissueType.ord", "idPot", "cycleValue")] %>%
  left_join(subset(dfCVall[, c("nameGen", "watering", "tissueType.ord", "idPot", "cycleValue")], tissueType.ord=="Leaf_8"), by="idPot")

ggplot(dfCV.L8_all, aes(x=cycleValue.x, y=cycleValue.y)) + geom_point()

