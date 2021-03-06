# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Cytometry of Arabidopsis leaf samples
# Project: [ARABREED]
# Subproject: plasticity of endopolyploidy [PLASTENDO] 2021
# Benoît Berthet Master 2 
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€

# source("cytometry.libraries.r")

source("./scripts/cytometry_IDsamples.r")

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Flow cytometry data ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
d.L30 <- read.flowSet(path=paste("/Users/denis/Documents/Encadrements/Stages/2021\ -\ M2\ -\ Benoit\ Berthet\ -\ Endopolyploidy/Experiment/endopolyploidy/data.Endopolyploidy.Nom_fichier_corrige_Leaf30only"))
d.flowset.L30 <- d.L30

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Flow cytometry data transformation ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
dat_L30 <- transform(d.L30, "lgDAPI"=log10(`DAPI`))

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Automated filtering of flow cytometry data by curve1Filter ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
resL30 <- filter(dat_L30[,], curv1Filter("lgDAPI", bwFac=0.5)) 
resSumL30 <- summary(resL30)

dfResL30 <- toTable(resSumL30)
names(dfResL30); head(dfResL30); str(dfResL30)
dfResL30$population <- as.factor(dfResL30$population)
dfResL30$population <- factor(dfResL30$population, labels=c("peak1", "peak2", "peak3", "peak4", "peak5","peak6", "peak7", "peak8", "peak9", "rest"))
dfResL30.wide <- reshape(dfResL30[, c("p", "sample", "population")], v.names = "p",  idvar = "sample", timevar = "population", direction = "wide")
names(dfResL30.wide)
dfResL30.wide <- do.call(data.frame, lapply(dfResL30.wide, function(x) replace(x, is.na(x), 0)))

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Calculate cycle value (endoreplication factor) ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
dfCV_L30 <- within(data = dfResL30.wide, { cycleValue <- (0*p.peak2 + 1*p.peak3 + 2*p.peak4 + 3*p.peak5 + 4*p.peak6 + 5*p.peak7 + 6*p.peak8 + 7*p.peak9)/(p.peak2 + p.peak3 + p.peak4 + p.peak5 + p.peak6 + p.peak7 + p.peak8 + p.peak9) })

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Split frames per filtering ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
split.resL30 <- split(dat_L30, resL30, population=list(keep=c("peak 1")))
Lth <- as.numeric(summary(split.resL30)[1])

df.resL30 <- data.frame(min1=rep(NA, Lth), max1=rep(NA, Lth), 
                     min2=rep(NA, Lth), max2=rep(NA, Lth),
                     min3=rep(NA, Lth), max3=rep(NA, Lth),
                     min4=rep(NA, Lth), max4=rep(NA, Lth),
                     min5=rep(NA, Lth), max5=rep(NA, Lth),
                     min6=rep(NA, Lth), max6=rep(NA, Lth),
                     min7=rep(NA, Lth), max7=rep(NA, Lth),
                     min8=rep(NA, Lth), max8=rep(NA, Lth)
                     )
for(i in 1:Lth){
  df.resL30[i, c("min1", "max1")] <- 10^(range(exprs(split.resL30$keep[[i]]$"lgDAPI")))
  }
split.resL30 <- split(dat_L30, resL30, population=list(keep=c("peak 2")))
for(i in 1:Lth){
  df.resL30[i, c("min2", "max2")] <- 10^(range(exprs(split.resL30$keep[[i]]$"lgDAPI")))}
split.resL30 <- split(dat_L30, resL30, population=list(keep=c("peak 3")))
for(i in 1:Lth){
  df.resL30[i, c("min3", "max3")] <- 10^(range(exprs(split.resL30$keep[[i]]$"lgDAPI")))}
split.resL30 <- split(dat_L30, resL30, population=list(keep=c("peak 4")))
for(i in 1:Lth){
  df.resL30[i, c("min4", "max4")] <- 10^(range(exprs(split.resL30$keep[[i]]$"lgDAPI")))}
split.resL30 <- split(dat_L30, resL30, population=list(keep=c("peak 5")))
for(i in 1:Lth){
  df.resL30[i, c("min5", "max5")] <- 10^(range(exprs(split.resL30$keep[[i]]$"lgDAPI")))}
split.resL30 <- split(dat_L30, resL30, population=list(keep=c("peak 6")))
for(i in 1:Lth){
  df.resL30[i, c("min6", "max6")] <- 10^(range(exprs(split.resL30$keep[[i]]$"lgDAPI")))}
split.resL30 <- split(dat_L30, resL30, population=list(keep=c("peak 7")))
for(i in 1:Lth){
  df.resL30[i, c("min7", "max7")] <- 10^(range(exprs(split.resL30$keep[[i]]$"lgDAPI")))}
split.resL30 <- split(dat_L30, resL30, population=list(keep=c("peak 8")))
for(i in 1:Lth){
  df.resL30[i, c("min8", "max8")] <- 10^(range(exprs(split.resL30$keep[[i]]$"lgDAPI")))}
split.resL30 <- split(dat_L30, resL30, population=list(keep=c("peak 9")))
for(i in 1:Lth){
  df.resL30[i, c("min9", "max9")] <- 10^(range(exprs(split.resL30$keep[[i]]$"lgDAPI")))}

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
#  Save df.res into d1.L30 for safe handling ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
d1.L30 <- df.resL30
d1.L30 <- do.call(data.frame,lapply(df.resL30, function(x) replace(x, is.infinite(x),NA)))
d1.L30 <- do.call(data.frame,lapply(d1.L30, function(x) replace(x, x==0,NA)))

d1.L30$names <- pData(d.L30[,])
d1.L30$num <- 1:Lth

d1.L30 <- d1.L30 %>% 
  mutate(name=as.character(names$name)) %>%
  left_join(IDS, by=c("name"="fileName"))

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Plot fluorescence peaks automatically detected by curve1filter() ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€

gg1 <- ggplot(d1.L30, aes(x=min1, y=num)) + scale_x_log10() + 
  geom_segment(aes(xend=max1, ystart=num, yend=num)) + 
  geom_segment(aes(x=min2, xend=max2, ystart=num, yend=num), col="red") +
  geom_segment(aes(x=min3, xend=max3, ystart=num, yend=num), col="blue") +
  geom_segment(aes(x=min4, xend=max4, ystart=num, yend=num), col="green") +
  geom_segment(aes(x=min5, xend=max5, ystart=num, yend=num), col="violet") +
  geom_segment(aes(x=min6, xend=max6, ystart=num, yend=num), col="grey") +
  geom_segment(aes(x=min7, xend=max7, ystart=num, yend=num), col="yellow") +
  geom_segment(aes(x=min8, xend=max8, ystart=num, yend=num), col="orange") +
  geom_segment(aes(x=min9, xend=max9, ystart=num, yend=num), col="black") +
  xlab("Fluorescence") + ylab("Sample") + theme_bw() +
  facet_wrap(.~watering, ncol=1)

pdf("./Figures/Fluorescence_peaks_L30.pdf", 7, 7)
gg1
dev.off()
system(paste("open", "./Figures/Fluorescence_peaks_L30.pdf"))

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 1/ calculate mean limits of each peak / day of measurement ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
df.m.L30 <- melt(d1.L30[, c("min1", "max1", "min2", "max2", "min3", "max3",
                            "min4", "max4", "min5", "max5", "min6", "max6", 
                            "min7", "max7", "min8", "max8", "min9", "max9")], 
             measure.vars = c("min1", "max1", "min2", "max2", "min3", "max3",
                              "min4", "max4", "min5", "max5", "min6", "max6",
                              "min7", "max7", "min8", "max8", "min9", "max9"))
df.m.L30$mM <- factor(df.m.L30$variable, labels=rep(c("min","max"), 9))

gp2.L30 <- ggplot(df.m.L30, aes(y=value, x=variable, fill=mM)) + geom_boxplot() + scale_y_log10(limits=c(1, 1000)) + ylab("") + xlab("Peak limits") + coord_flip()

peak.limits.L30 <- aggregate(value ~ variable, data=df.m.L30, FUN = summary)

# Plot limits 
gp2.L30 + geom_hline(data=data.frame(yint=peak.limits.L30[, 2][, "1st Qu."]), aes(yintercept=yint), col="gray") + 
  geom_hline(data=data.frame(yint=peak.limits.L30[, 2][, "3rd Qu."]), aes(yintercept=yint), col="gray") + 
  theme_bw()

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 1.a/ Define the limits using 1st and 3rd quartiles ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
pl.L30 <- peak.limits.L30[, 2]
peak.lim.L30 <- list(
  peak1=c(pl.L30[1, "Min."], mean(c(pl.L30[2, "3rd Qu."], pl.L30[3,"1st Qu."]))),
  peak2=c(mean(c(pl.L30[2, "3rd Qu."], pl.L30[3,"1st Qu."])), 
          mean(c(pl.L30[4, "3rd Qu."], pl.L30[5,"1st Qu."]))),
  peak3=c(mean(c(pl.L30[4, "3rd Qu."], pl.L30[5,"1st Qu."])), 
          mean(c(pl.L30[6, "3rd Qu."], pl.L30[7,"1st Qu."]))),
  peak4=c(mean(c(pl.L30[6, "3rd Qu."], pl.L30[7,"1st Qu."])),
          mean(c(pl.L30[8, "3rd Qu."], pl.L30[9,"1st Qu."]))),
  peak5=c(mean(c(pl.L30[8, "3rd Qu."], pl.L30[9,"1st Qu."])),
          mean(c(pl.L30[10, "3rd Qu."], pl.L30[11,"1st Qu."]))),
  peak6=c(mean(c(pl.L30[10, "3rd Qu."], pl.L30[11,"1st Qu."])),
          mean(c(pl.L30[12, "3rd Qu."], pl.L30[13,"1st Qu."]))),
  peak7=c(mean(c(pl.L30[12, "3rd Qu."], pl.L30[13,"1st Qu."])), 
          mean(c(pl.L30[14, "3rd Qu."], pl.L30[15,"1st Qu."]))),
  peak8=c(mean(c(pl.L30[14, "3rd Qu."], pl.L30[15,"1st Qu."])), 
          mean(c(pl.L30[16, "3rd Qu."], pl.L30[17,"1st Qu."]))),
  peak9=c(mean(c(pl.L30[16, "3rd Qu."], pl.L30[17,"1st Qu."])),
          pl.L30[18, "Max."])
  )
pl1.L30 <- melt(as.data.frame(peak.lim.L30))
gp2.L30 + geom_hline(data=pl1.L30, aes(yintercept=value), col="blue") + theme_bw()

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 1.b/ Plot global peak limits ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
pdf("./figures/estimatedPeakLimits_boxplots.L30.pdf", 7, 7)
gp2.L30 + geom_hline(data=pl1.L30, aes(yintercept=value), col="blue") + theme_bw() 
dev.off()
system("open ./figures/estimatedPeakLimits_boxplots.L30.pdf")

system("open ../figures/detected_and_estimatedPeakLimits.pdf")
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 2/ use these limits as filters in filter() ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
rectGate1.L30 <- rectangleGate(filterId = "peak1", "DAPI"=peak.lim.L30$peak1)

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 3/ use these limits to extract peak % ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 
fres00.L30 <- filter(dat_L30, filter=rectGate1.L30)
cycleProportion.L30 <- toTable(summary(fres00.L30))
for(i in c("peak2", "peak3", "peak4", "peak5", "peak6", "peak7", "peak8", "peak9")) {
  rectGate.i <- rectangleGate(filterId = i, "DAPI"=peak.lim.L30[[i]])
  fres00.L30 <- filter(dat_L30, filter=rectGate.i)
  cycleProportion.L30 <- rbind(cycleProportion.L30, toTable(summary(fres00.L30)))
  }

cycleProportion.L30$ploidy <- factor(cycleProportion.L30$population, labels=c("debris", "x2C", "x4C", "x8C", "x16C", "x32C", "x64C", "x128C", "x256C"))
cycleProportion.L30.wide <- dcast(cycleProportion.L30[, c("sample", "ploidy", "p")], formula=sample~ploidy, value.var="p")
cycleProportion.L30.wide[is.na(cycleProportion.L30.wide$x128C), "x128C"] <- 0
cycleProportion.L30.wide[is.na(cycleProportion.L30.wide$x128C), "x256C"] <- 0
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Calculate new cycle value (endoreduplication factor) ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
dfCV.L30 <- within(data = cycleProportion.L30.wide, {
  cycleValue <- (0*x2C + 1*x4C + 2*x8C + 3*x16C + 4*x32C + 5*x64C + 6*x128C + 7*x256C)/(x2C + x4C + x8C + x16C + x32C + x64C + x128C + x256C)
})
dfCV.L30$cycleValue
dfCV.L30 <- dfCV.L30 %>%
  left_join(ids, by=c("sample"="fileName")) %>%
  left_join(idPots, by=c("idPot"))

dfCV.L30$tissueType.ord <- factor(dfCV.L30$tissueType, 
                             levels = c("f6", "f8", "F8", "f30"), 
                             labels = c("Seedling_Leaf5","Leaf_8", "Leaf_8", "Leaf_30"))

dfCV.L30[dfCV.L30$tissueType.ord=="Seedling_Leaf5", "watering"] <- "WW"
dfCV.L30$watering <- factor(dfCV.L30$watering, levels = c("WW", "WD"))

CV.mean.L30 <- subset(dfCV.L30) %>%
  mutate(p2C=x2C/(x2C+x4C+x8C+x16C+x32C+x64C+x128C+x256C),
         p4C=x4C/(x2C+x4C+x8C+x16C+x32C+x64C+x128C+x256C),
         p8C=x8C/(x2C+x4C+x8C+x16C+x32C+x64C+x128C+x256C),
         p16C=x16C/(x2C+x4C+x8C+x16C+x32C+x64C+x128C+x256C),
         p32C=x32C/(x2C+x4C+x8C+x16C+x32C+x64C+x128C+x256C),
         p64C=x64C/(x2C+x4C+x8C+x16C+x32C+x64C+x128C+x256C),
         p128C=x128C/(x2C+x4C+x8C+x16C+x32C+x64C+x128C+x256C),
         p256C=x256C/(x2C+x4C+x8C+x16C+x32C+x64C+x128C+x256C)) %>%
  select(idPot, nameGen, tissueType.ord, watering, cycleValue, p2C,p4C,p8C,p16C,p32C,p64C,p128C,p256C) %>%
  group_by(nameGen, tissueType.ord, watering) %>%
  summarize(CVmean = mean(cycleValue),
            p2C = mean(p2C/(p2C+p4C+p8C+p16C+p32C+p64C+p128C+p256C)),
            p4C = mean(p4C/(p2C+p4C+p8C+p16C+p32C+p64C+p128C+p256C)),
            p8C = mean(p8C/(p2C+p4C+p8C+p16C+p32C+p64C+p128C+p256C)),
            p16C = mean(p16C/(p2C+p4C+p8C+p16C+p32C+p64C+p128C+p256C)),
            p32C = mean(p32C/(p2C+p4C+p8C+p16C+p32C+p64C+p128C+p256C)),
            p64C = mean(p64C/(p2C+p4C+p8C+p16C+p32C+p64C+p128C+p256C)),
            p128C = mean(p128C/(p2C+p4C+p8C+p16C+p32C+p64C+p128C+p256C)),
            p256C = mean(p256C/(p2C+p4C+p8C+p16C+p32C+p64C+p128C+p256C))) %>%
  mutate(p2C=p2C/(p2C+p4C+p8C+p16C+p32C+p64C+p128C+p256C),
         p4C=p4C/(p2C+p4C+p8C+p16C+p32C+p64C+p128C+p256C),
         p8C=p8C/(p2C+p4C+p8C+p16C+p32C+p64C+p128C+p256C),
         p16C=p16C/(p2C+p4C+p8C+p16C+p32C+p64C+p128C+p256C),
         p32C=p32C/(p2C+p4C+p8C+p16C+p32C+p64C+p128C+p256C),
         p64C=p64C/(p2C+p4C+p8C+p16C+p32C+p64C+p128C+p256C),
         p128C=p128C/(p2C+p4C+p8C+p16C+p32C+p64C+p128C+p256C),
         p256C=p256C/(p2C+p4C+p8C+p16C+p32C+p64C+p128C+p256C)) 

CV.mean.L30.wide <- dcast(CV.mean.L30, formula=nameGen~tissueType.ord + watering, mean, value.var = "CVmean")

CV.mean.L30.long <- CV.mean.L30 %>%
  select(nameGen, tissueType.ord, watering, p2C,p4C,p8C,p16C,p32C,p64C,p128C,p256C) %>%
  tidyr::pivot_longer(cols=c(p2C,p4C,p8C,p16C,p32C,p64C,p128C,p256C)) %>%
  mutate(name=factor(name, levels=c("p2C","p4C","p8C","p16C","p32C","p64C","p128C","p256C"),
                     labels=c("2C","4C","8C","16C","32C","64C","128C","256C")))

CV.mean.L30.long$nameGen.OrderedL30_WW <- as.factor(CV.mean.L30.long$nameGen)
CV.mean.L30.long$nameGen.OrderedL30_WW <- factor(CV.mean.L30.long$nameGen.OrderedL30_WW,
                              levels =levels(CV.mean.L30.long$nameGen.OrderedL30_WW)[order(CV.mean.L30.wide$Leaf_30_WW)])
         
gp.percentNuclei <- ggplot(subset(CV.mean.L30.long, watering=="WW"), aes(y=value, x=nameGen.OrderedL30_WW, fill=forcats::fct_rev(name))) + 
  geom_bar(position="fill", stat="identity", width=1) +
  ylab("% of nuclei in leaf #30") + xlab("") +
  coord_cartesian(ylim=c(0,1), expand = F) +
  facet_wrap(vars(
    # Change factor level name
    fct_recode(watering, "Well-watered condition" = "WW"))
    ) +
  scale_fill_grey("Size class") +
  myTheme + 
  theme(legend.position = "right",
        axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5),
        strip.background = element_rect(fill="transparent"),
        strip.text = element_text(size=14),
        panel.border=element_rect(fill="transparent", size=0.75))
pdf("./figures/percent-nuclei_L30_WW.pdf", 8, 7)
gp.percentNuclei
dev.off()
system("open ./figures/percent-nuclei_L30_WW.pdf")

dfCV.L30$nameGen.OrderedL30_WW <- as.factor(dfCV.L30$nameGen)
dfCV.L30$nameGen.OrderedL30_WW <- factor(dfCV.L30$nameGen.OrderedL30_WW,
                                          levels =levels(dfCV.L30$nameGen.OrderedL30_WW)[order(CV.mean.L30.wide$Leaf_30_WW)])

pdf("./figures/cycleValue_30_genotypes.L30.pdf", 12, 8)

ggplot(data=subset(dfCV.L30, tissueType.ord%in%c("Leaf_30")), 
       aes(y=cycleValue, x=nameGen.OrderedL30_WW, fill=watering)) +
  geom_boxplot(outlier.alpha = 0) +
  xlab("Accessions (ordered by L30_WW CV)") +
  ylab("Cycle value") +
  scale_fill_manual(values=c("#386FA4", "#84D2F6")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))

dev.off()
system("open ./figures/cycleValue_30_genotypes.L30.pdf")

gp.corr <- ggplot(data=CV.mean.wide, aes(x = Leaf_30_WW, y = Leaf_30_WD)) +
  geom_point() + geom_smooth(method = lm, se=F)

pdf("./figures/cycleValue_30_genotypes_correlations.L30_vs_all.pdf", 8, 7)
gp.corr + geom_abline(slope = 1, intercept = 0)
dev.off()
system("open ./figures/cycleValue_30_genotypes_correlations.L30_vs_all.pdf")

dfCV.L30_all <- dfCV.L30[, c("nameGen", "watering", "tissueType.ord", "idPot", "cycleValue")] %>%
  left_join(subset(dfCVall[1:622, c("nameGen", "watering", "tissueType.ord", "idPot", "cycleValue")], tissueType.ord=="Leaf_30"), by="idPot")

pdf("./figures/cycleValue_correlations.L30_vs_all_bw_fac_0.8.pdf", 8, 7)
ggplot(dfCV.L30_all, aes(x=cycleValue.x, y=cycleValue.y)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  xlab("Cycle value with gating on L30 samples") +
  ylab("Cycle value with gating on all samples") +
  theme(text = element_text(size=16))
dev.off()
system("open ./figures/cycleValue_correlations.L30_vs_all_bw_fac_0.8.pdf")

dfCV.L30 %>%
  select(cycleValue, idPot, nameGen, tissueType.ord) %>%
  left_join(subset(d.SLA, stage == "L30" ), by=c("idPot"="idPot")) %>%
  ggplot(., aes(y=leaf.area.mm2, x=cycleValue)) + geom_point() + geom_smooth(method=lm)
