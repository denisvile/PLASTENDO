# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Cytometry of Arabidopsis leaf samples
# Project: [ARABREED]
# Subproject: plasticity of endopolyploidy [PLASTENDO] 2021
# Benoît Berthet Master 2 
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€

source("./scripts/cytometry.libraries.r")

source("./scripts/cytometry_IDsamples.r")

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Flow cytometry data ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
d.Lsdlg <- read.flowSet(path=paste("/Users/denis/Documents/Encadrements/Stages/2021\ -\ M2\ -\ Benoit\ Berthet\ -\ Endopolyploidy/Experiment/endopolyploidy/data.Endopolyploidy.Nom_fichier_corrige_LeavesSdlgOnly"))
d.flowset.Lsdlg <- d.Lsdlg

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Flow cytometry data transformation ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
dat_Lsdlg <- transform(d.Lsdlg, "lgDAPI"=log10(`DAPI`))

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Automated filtering of flow cytometry data by curve1Filter ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
resLsdlg <- filter(dat_Lsdlg[,], curv1Filter("lgDAPI", bwFac=0.5)) 
resSumLsdlg <- summary(resLsdlg)

dfResLsdlg <- toTable(resSumLsdlg)
names(dfResLsdlg); head(dfResLsdlg); str(dfResLsdlg)
dfResLsdlg$population <- as.factor(dfResLsdlg$population)
dfResLsdlg$population <- factor(dfResLsdlg$population, labels=c("peak1", "peak2", "peak3", "peak4", "peak5","peak6", "peak7", "rest"))
dfResLsdlg.wide <- reshape(dfResLsdlg[, c("p", "sample", "population")], v.names = "p",  idvar = "sample", timevar = "population", direction = "wide")
names(dfResLsdlg.wide)
dfResLsdlg.wide <- do.call(data.frame, lapply(dfResLsdlg.wide, function(x) replace(x, is.na(x), 0)))

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Calculate cycle value (endoreplication factor) ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
dfCV_Lsdlg <- within(data = dfResLsdlg.wide, { cycleValue <- (0*p.peak2 + 1*p.peak3 + 2*p.peak4 + 3*p.peak5 + 4*p.peak6 + 5*p.peak7)/(p.peak2 + p.peak3 + p.peak4 + p.peak5 + p.peak6 + p.peak7) })

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Split frames per filtering ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
split.resLsdlg <- split(dat_Lsdlg, resLsdlg, population=list(keep=c("peak 1")))
Lth_Lsdlg <- as.numeric(summary(split.resLsdlg)[1])

df.resLsdlg <- data.frame(min1=rep(NA, Lth_Lsdlg), max1=rep(NA, Lth_Lsdlg), 
                     min2=rep(NA, Lth_Lsdlg), max2=rep(NA, Lth_Lsdlg),
                     min3=rep(NA, Lth_Lsdlg), max3=rep(NA, Lth_Lsdlg),
                     min4=rep(NA, Lth_Lsdlg), max4=rep(NA, Lth_Lsdlg),
                     min5=rep(NA, Lth_Lsdlg), max5=rep(NA, Lth_Lsdlg),
                     min6=rep(NA, Lth_Lsdlg), max6=rep(NA, Lth_Lsdlg),
                     min7=rep(NA, Lth_Lsdlg), max7=rep(NA, Lth_Lsdlg)
                     )
for(i in 1:Lth_Lsdlg){
  df.resLsdlg[i, c("min1", "max1")] <- 10^(range(exprs(split.resLsdlg$keep[[i]]$"lgDAPI")))}
split.resLsdlg <- split(dat_Lsdlg, resLsdlg, population=list(keep=c("peak 2")))
for(i in 1:Lth_Lsdlg){
  df.resLsdlg[i, c("min2", "max2")] <- 10^(range(exprs(split.resLsdlg$keep[[i]]$"lgDAPI")))}
split.resLsdlg <- split(dat_Lsdlg, resLsdlg, population=list(keep=c("peak 3")))
for(i in 1:Lth_Lsdlg){
  df.resLsdlg[i, c("min3", "max3")] <- 10^(range(exprs(split.resLsdlg$keep[[i]]$"lgDAPI")))}
split.resLsdlg <- split(dat_Lsdlg, resLsdlg, population=list(keep=c("peak 4")))
for(i in 1:Lth_Lsdlg){
  df.resLsdlg[i, c("min4", "max4")] <- 10^(range(exprs(split.resLsdlg$keep[[i]]$"lgDAPI")))}
split.resLsdlg <- split(dat_Lsdlg, resLsdlg, population=list(keep=c("peak 5")))
for(i in 1:Lth_Lsdlg){
  df.resLsdlg[i, c("min5", "max5")] <- 10^(range(exprs(split.resLsdlg$keep[[i]]$"lgDAPI")))}
split.resLsdlg <- split(dat_Lsdlg, resLsdlg, population=list(keep=c("peak 6")))
for(i in 1:Lth_Lsdlg){
  df.resLsdlg[i, c("min6", "max6")] <- 10^(range(exprs(split.resLsdlg$keep[[i]]$"lgDAPI")))}
split.resLsdlg <- split(dat_Lsdlg, resLsdlg, population=list(keep=c("peak 7")))
for(i in 1:Lth_Lsdlg){
  df.resLsdlg[i, c("min7", "max7")] <- 10^(range(exprs(split.resLsdlg$keep[[i]]$"lgDAPI")))}

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
#  Save df.res into d1.Lsdlg for safe handling ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
d1.Lsdlg <- df.resLsdlg
d1.Lsdlg <- do.call(data.frame,lapply(df.resLsdlg, function(x) replace(x, is.infinite(x),NA)))
d1.Lsdlg <- do.call(data.frame,lapply(d1.Lsdlg, function(x) replace(x, x==0,NA)))

d1.Lsdlg$names <- pData(d.Lsdlg[,])
d1.Lsdlg$num <- 1:Lth_Lsdlg

d1.Lsdlg <- d1.Lsdlg %>% 
  mutate(name=as.character(names$name)) %>%
  left_join(IDS, by=c("name"="fileName"))

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Plot fluorescence peaks automatically detected by curve1filter() ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€

gg1 <- ggplot(d1.Lsdlg, aes(x=min1, y=num)) + scale_x_log10() + 
  geom_segment(aes(xend=max1, ystart=num, yend=num)) + 
  geom_segment(aes(x=min2, xend=max2, ystart=num, yend=num), col="red") +
  geom_segment(aes(x=min3, xend=max3, ystart=num, yend=num), col="blue") +
  geom_segment(aes(x=min4, xend=max4, ystart=num, yend=num), col="green") +
  geom_segment(aes(x=min5, xend=max5, ystart=num, yend=num), col="violet") +
  geom_segment(aes(x=min6, xend=max6, ystart=num, yend=num), col="grey") +
  geom_segment(aes(x=min7, xend=max7, ystart=num, yend=num), col="yellow") +
  xlab("Fluorescence") + ylab("Sample") + theme_bw() +
  facet_wrap(.~watering, ncol=1)

pdf("./Figures/Fluorescence_peaks_Lsdlg.pdf", 7, 7)
gg1
dev.off()
system(paste("open", "./Figures/Fluorescence_peaks_Lsdlg.pdf"))

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 1/ calculate mean limits of each peak / day of measurement ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
df.m.Lsdlg <- melt(d1.Lsdlg[, c("min1", "max1", "min2", "max2", "min3", "max3",
                            "min4", "max4", "min5", "max5", "min6", "max6", 
                            "min7", "max7")], 
             measure.vars = c("min1", "max1", "min2", "max2", "min3", "max3",
                              "min4", "max4", "min5", "max5", "min6", "max6",
                              "min7", "max7"))
df.m.Lsdlg$mM <- factor(df.m.Lsdlg$variable, labels=rep(c("min","max"), 7))

gp2.Lsdlg <- ggplot(df.m.Lsdlg, aes(y=value, x=variable, fill=mM)) + geom_boxplot() + scale_y_log10(limits=c(1, 1000)) + ylab("") + xlab("Peak limits") + coord_flip()

peak.limits.Lsdlg <- aggregate(value ~ variable, data=df.m.Lsdlg, FUN = summary)

# Plot limits 
gp2.Lsdlg + geom_hline(data=data.frame(yint=peak.limits.Lsdlg[, 2][, "1st Qu."]), aes(yintercept=yint), col="gray") + 
  geom_hline(data=data.frame(yint=peak.limits.Lsdlg[, 2][, "3rd Qu."]), aes(yintercept=yint), col="gray") + 
  theme_bw()

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 1.a/ Define the limits using 1st and 3rd quartiles ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
pl.Lsdlg <- peak.limits.Lsdlg[, 2]
peak.lim.Lsdlg <- list(
  peak1=c(pl.Lsdlg[1, "Min."], mean(c(pl.Lsdlg[2, "3rd Qu."], pl.Lsdlg[3,"1st Qu."]))),
  peak2=c(mean(c(pl.Lsdlg[2, "3rd Qu."], pl.Lsdlg[3,"1st Qu."])), 
          mean(c(pl.Lsdlg[4, "3rd Qu."], pl.Lsdlg[5,"1st Qu."]))),
  peak3=c(mean(c(pl.Lsdlg[4, "3rd Qu."], pl.Lsdlg[5,"1st Qu."])), 
          mean(c(pl.Lsdlg[6, "3rd Qu."], pl.Lsdlg[7,"1st Qu."]))),
  peak4=c(mean(c(pl.Lsdlg[6, "3rd Qu."], pl.Lsdlg[7,"1st Qu."])),
          mean(c(pl.Lsdlg[8, "3rd Qu."], pl.Lsdlg[9,"1st Qu."]))),
  peak5=c(mean(c(pl.Lsdlg[8, "3rd Qu."], pl.Lsdlg[9,"1st Qu."])),
          mean(c(pl.Lsdlg[10, "3rd Qu."], pl.Lsdlg[11,"1st Qu."]))),
  peak6=c(mean(c(pl.Lsdlg[10, "3rd Qu."], pl.Lsdlg[11,"1st Qu."])),
          mean(c(pl.Lsdlg[12, "3rd Qu."], pl.Lsdlg[13,"1st Qu."]))),
  peak7=c(mean(c(pl.Lsdlg[12, "3rd Qu."], pl.Lsdlg[13,"1st Qu."])),
          pl.Lsdlg[14, "Max."])
  )
pl1.Lsdlg <- melt(as.data.frame(peak.lim.Lsdlg))
gp2.Lsdlg + geom_hline(data=pl1.Lsdlg, aes(yintercept=value), col="blue") + theme_bw()

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 1.b/ Plot global peak limits ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
pdf("./figures/estimatedPeakLimits_boxplots.Lsdlg.pdf", 7, 7)
gp2.Lsdlg + geom_hline(data=pl1.Lsdlg, aes(yintercept=value), col="blue") + theme_bw() 
dev.off()
system("open ./figures/estimatedPeakLimits_boxplots.Lsdlg.pdf")

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 2/ use these limits as filters in filter() ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
rectGate1.Lsdlg <- rectangleGate(filterId = "peak1", "DAPI"=peak.lim.Lsdlg$peak1)

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 3/ use these limits to extract peak % ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 
fres00.Lsdlg <- filter(dat_Lsdlg, filter=rectGate1.Lsdlg)
cycleProportion.Lsdlg <- toTable(summary(fres00.Lsdlg))
for(i in c("peak2", "peak3", "peak4", "peak5", "peak6", "peak7")) {
  rectGate.i <- rectangleGate(filterId = i, "DAPI"=peak.lim.Lsdlg[[i]])
  fres00.Lsdlg <- filter(dat_Lsdlg, filter=rectGate.i)
  cycleProportion.Lsdlg <- rbind(cycleProportion.Lsdlg, toTable(summary(fres00.Lsdlg)))
  }

cycleProportion.Lsdlg$ploidy <- factor(cycleProportion.Lsdlg$population, labels=c("debris", "x2C", "x4C", "x8C", "x16C", "x32C", "x64C"))
cycleProportion.Lsdlg.wide <- dcast(cycleProportion.Lsdlg[, c("sample", "ploidy", "p")], formula=sample~ploidy, value.var="p")
cycleProportion.Lsdlg.wide[is.na(cycleProportion.Lsdlg.wide$x128C), "x128C"] <- 0
cycleProportion.Lsdlg.wide[is.na(cycleProportion.Lsdlg.wide$x128C), "x256C"] <- 0
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Calculate new cycle value (endoreduplication factor) ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
dfCV.Lsdlg <- within(data = cycleProportion.Lsdlg.wide, {
  cycleValue <- (0*x2C + 1*x4C + 2*x8C + 3*x16C + 4*x32C + 5*x64C)/(x2C + x4C + x8C + x16C + x32C + x64C)
})
dfCV.Lsdlg$cycleValue
dfCV.Lsdlg <- dfCV.Lsdlg %>%
  left_join(ids, by=c("sample"="fileName")) %>%
  left_join(idPots, by=c("idPot"))

dfCV.Lsdlg$tissueType.ord <- factor(dfCV.Lsdlg$tissueType, 
                             levels = c("f6", "f8", "F8", "f30"), 
                             labels = c("Seedling_Leaf5","Leaf_8", "Leaf_8", "Leaf_30"))

dfCV.Lsdlg[dfCV.Lsdlg$tissueType.ord=="Seedling_Leaf5", "watering"] <- "WW"
dfCV.Lsdlg$watering <- factor(dfCV.Lsdlg$watering, levels = c("WW", "WD"))

CV.mean.Lsdlg <- dfCV.Lsdlg %>%
  select(nameGen, tissueType.ord, watering, cycleValue) %>%
  group_by(nameGen, tissueType.ord, watering) %>%
  summarize(CVmean = mean(cycleValue))

CV.mean.Lsdlg.wide <- dcast(CV.mean.Lsdlg, formula=nameGen~tissueType.ord + watering, mean, value.var = "CVmean")

dfCV.Lsdlg$nameGen.OrderedLsdlg_WW <- as.factor(dfCV.Lsdlg$nameGen)
dfCV.Lsdlg$nameGen.OrderedLsdlg_WW <- factor(dfCV.Lsdlg$nameGen.OrderedLsdlg_WW,
                                          levels =levels(dfCV.Lsdlg$nameGen.OrderedLsdlg_WW)[order(CV.mean.Lsdlg.wide$Seedling_Leaf5_WW)])


CV.mean.Lsdlg.wide$clust <- as.factor(kmeans(CV.mean.Lsdlg.wide[, "Seedling_Leaf5_WW"], centers=5)$cluster)
dfCV.Lsdlg <- dfCV.Lsdlg %>%
  left_join(CV.mean.Lsdlg.wide)





