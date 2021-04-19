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
d <- read.flowSet(path=paste("/Users/denis/Documents/Encadrements/Stages/2021\ -\ M2\ -\ Benoit\ Berthet\ -\ Endopolyploidy/Experiment/endopolyploidy/data.Endopolyploidy.Nom_fichier_corrige"))

d.flowset <- d
str(d.flowset)
str(pData(d.flowset))
#write.table(dir("../data/FCSexports/rosette_Leaf/"), file="idSamples.csv")
phenoData(d.flowset)

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Test, data handling... ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
dff <-  d[1, "DAPI"][[1]]
# Modifying the minRange of the data
pData(parameters(dff))[,4] <- min(exprs(dff))

head(dff)
dff$"DAPI"  # Equivalent to dff[, "DAPI"]

# Extract date of measurement
keyword(dff)$`$DATE`

#summary(dff)
#plot(dff)
ddf <- as.data.frame(exprs(dff))
str(ddf)
names(ddf) <- "DAPI"
log10(ddf$DAPI)
ggplot(ddf, aes(x=log(DAPI))) + geom_histogram(bins = 100) + 
  scale_x_log10(limits=c(1, 10^3))
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Flow cytometry data transformation ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
dat <- transform(d, "lgDAPI"=log10(`DAPI`))

pdf("./figures/Distribution_DAPI.pdf")
for(i in 1:length(dir("/Users/denis/Documents/Encadrements/Stages/2021\ -\ M2\ -\ Benoit\ Berthet\ -\ Endopolyploidy/Experiment/endopolyploidy/data.Endopolyploidy/Ag-0/F6/"))) {
  dd <- dat[i, c("DAPI", "lgDAPI")][[1]]
  minR <- min(exprs(dd))
  #pData(parameters(dd))[4] <- log10(minR)
  df <- as.data.frame(exprs(dd))
  names(df) <- c("DAPI", "lgDAPI")
  gg <- ggplot(df, aes(x=DAPI)) + geom_histogram(aes(y=..density..), bins = 100) + 
    geom_density() + scale_x_log10(limits=c(1, 10^3)) + ggtitle(pData(d)[i, ])
  print(gg)
  print(densityplot(~ `lgDAPI`, dd, filter=curv1Filter("lgDAPI")))
}
dev.off()
system(paste("open","./figures/Distribution_DAPI.pdf"))

km <- kmeans(df[, "lgDAPI"],centers=5)
df$clust <- as.factor(km$cluster)
ggplot(df, aes(x=lgDAPI))  + geom_density(aes(fill=clust)) + geom_density()#, binwidth=0.05, color="grey50")# + stat_density(geom="line", color="red")

ggplot(df, aes(x=lgDAPI)) + 
  geom_histogram(aes(fill=clust,y=..count../sum(..count..)), 
                 binwidth=0.06, color="grey50") + 
  stat_density(geom="line", color="red")

pdf("./figures/Distribution_DAPI_filter.pdf")
densityplot(~ `lgDAPI`, dat[4,][[1]], filter=curv1Filter("lgDAPI", bwFac=0.8, gridsize = rep(401,2)), xlim=c(0, 4))
dev.off()
system(paste("open","./figures/Distribution_DAPI_filter.pdf"))

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Automated filtering of flow cytometry data by curve1Filter ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Global ----
res <- filter(dat[1:622,], curv1Filter("lgDAPI", bwFac=0.80)) # ID 623 excluded because not obtained with same scale (and it is a duplicate sample)
resSum <- summary(res)

dfRes <- toTable(resSum)
names(dfRes); head(dfRes); str(dfRes)
dfRes$population <- as.factor(dfRes$population)
dfRes$population <- factor(dfRes$population, labels=c("peak1", "peak2", "peak3", "peak4", "peak5","peak6", "peak7", "peak8", "rest"))
dfRes.wide <- reshape(dfRes[, c("p", "sample", "population")], v.names = "p",  idvar = "sample", timevar = "population", direction = "wide")
names(dfRes.wide)
dfRes.wide <- do.call(data.frame, lapply(dfRes.wide, function(x) replace(x, is.na(x), 0)))

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Calculate cycle value (endoreplication factor) ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
dfCV <- within(data = dfRes.wide, { cycleValue <- (0*p.peak1 + 1*p.peak2 + 2*p.peak3 + 3*p.peak4 + 4*p.peak5 + 5*p.peak6)/(p.peak1 + p.peak2 + p.peak3 + p.peak4 + p.peak5 + p.peak6) })

ggplot(subset(dfRes), aes(y=percent, x=population)) + geom_boxplot()

ggplot(subset(dfCV), aes(x=cycleValue)) + geom_histogram()

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Split frames per filtering ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# split(dat, res)

split.res <- split(dat[1:622,], res, population=list(keep=c("peak 1")))
Lth <- as.numeric(summary(split.res)[1])

df.res <- data.frame(min1=rep(NA, Lth), max1=rep(NA, Lth), 
                     min2=rep(NA, Lth), max2=rep(NA, Lth),
                     min3=rep(NA, Lth), max3=rep(NA, Lth),
                     min4=rep(NA, Lth), max4=rep(NA, Lth),
                     min5=rep(NA, Lth), max5=rep(NA, Lth),
                     min6=rep(NA, Lth), max6=rep(NA, Lth),
                     min7=rep(NA, Lth), max7=rep(NA, Lth)
)

for(i in 1:Lth){
  df.res[i, c("min1", "max1")] <- 10^(range(exprs(split.res$keep[[i]]$"lgDAPI")))
}
split.res <- split(dat[1:622,], res, population=list(keep=c("peak 2")))
for(i in 1:Lth){
  df.res[i, c("min2", "max2")] <- 10^(range(exprs(split.res$keep[[i]]$"lgDAPI")))
}
split.res <- split(dat[1:622,], res, population=list(keep=c("peak 3")))
for(i in 1:Lth){
  df.res[i, c("min3", "max3")] <- 10^(range(exprs(split.res$keep[[i]]$"lgDAPI")))
}
split.res <- split(dat[1:622,], res, population=list(keep=c("peak 4")))
for(i in 1:Lth){
  df.res[i, c("min4", "max4")] <- 10^(range(exprs(split.res$keep[[i]]$"lgDAPI")))
}
split.res <- split(dat[1:622,], res, population=list(keep=c("peak 5")))
for(i in 1:Lth){
  df.res[i, c("min5", "max5")] <- 10^(range(exprs(split.res$keep[[i]]$"lgDAPI")))
}
split.res <- split(dat[1:622,], res, population=list(keep=c("peak 6")))
for(i in 1:Lth){
  df.res[i, c("min6", "max6")] <- 10^(range(exprs(split.res$keep[[i]]$"lgDAPI")))
}
split.res <- split(dat[1:622,], res, population=list(keep=c("peak 7")))
for(i in 1:Lth){
  df.res[i, c("min7", "max7")] <- 10^(range(exprs(split.res$keep[[i]]$"lgDAPI")))
}
split.res <- split(dat[1:622,], res, population=list(keep=c("peak 8")))
for(i in 1:Lth){
  df.res[i, c("min8", "max8")] <- 10^(range(exprs(split.res$keep[[i]]$"lgDAPI")))
}
#split.res <- split(dat[1:622,], res, population=list(keep=c("peak 9")))
#for(i in 1:Lth){
#  df.res[i, c("min9", "max9")] <- 10^(range(exprs(split.res$keep[[i]]$"lgDAPI")))
#}

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
#  ave df.res into d1 for safe handling ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€

d1 <- df.res

d1 <- do.call(data.frame,lapply(df.res, function(x) replace(x, is.infinite(x),NA)))
d1 <- do.call(data.frame,lapply(d1, function(x) replace(x, x==0,NA)))

d1$names <- pData(d[1:622,])
d1$num <- 1:Lth
###### Les identifiants ne sont pas encore faits !!!!!!!!!!!!!! A FAIRE
#d1 <- cbind(d1[1:Lth, ], ids)
#d1$idAccession <- as.factor(d1$idAccession)
#d1 <- merge(d1, measureDate, by.x="num", by.y="id")

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Plot fluorescence peaks automatically detected by curve1filter() ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
pdf("./Figures/Fluorescence_peaks.pdf", 7, 7)
gg1 <- ggplot(d1, aes(x=min1, y=num)) + scale_x_log10() + 
  geom_segment(aes(xend=max1, ystart=num, yend=num)) + 
  geom_segment(aes(x=min2, xend=max2, ystart=num, yend=num), col="red") +
  geom_segment(aes(x=min3, xend=max3, ystart=num, yend=num), col="blue") +
  geom_segment(aes(x=min4, xend=max4, ystart=num, yend=num), col="green") +
  geom_segment(aes(x=min5, xend=max5, ystart=num, yend=num), col="violet") +
  geom_segment(aes(x=min6, xend=max6, ystart=num, yend=num), col="grey") +
  geom_segment(aes(x=min7, xend=max7, ystart=num, yend=num), col="yellow") +
  geom_segment(aes(x=min8, xend=max8, ystart=num, yend=num), col="orange") +
 # geom_segment(aes(x=min9, xend=max9, ystart=num, yend=num), col="black") +
  xlab("Fluorescence") + ylab("Sample") + theme_bw()
print(gg1)
#print(gg1 + facet_wrap(~idAccession))
#print(gg1 + facet_wrap(~bloc))
#print(gg1 + facet_grid(date~.) + theme_bw() + ylab("Sample (panels ordered by date of measurement)"))
dev.off()
system(paste("open", "./Figures/Fluorescence_peaks.pdf"))

gg2 <- ggplot(d1, aes(x=min1, y=num, colour=date)) + scale_x_log10() + 
  geom_segment(aes(xend=max1, ystart=num, yend=num)) + 
  geom_segment(aes(x=min2, xend=max2, ystart=num, yend=num)) +
  geom_segment(aes(x=min3, xend=max3, ystart=num, yend=num)) +
  geom_segment(aes(x=min4, xend=max4, ystart=num, yend=num)) +
  geom_segment(aes(x=min5, xend=max5, ystart=num, yend=num)) +
  geom_segment(aes(x=min6, xend=max6, ystart=num, yend=num)) +
  geom_segment(aes(x=min7, xend=max7, ystart=num, yend=num)) +
  xlab("Fluorescence") + ylab("Sample") + theme_bw()
print(gg2)


# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 1/ calculate mean limits of each peak / day of measurement ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
df.m <- melt(d1[, c("min1", "max1", "min2", "max2", "min3", "max3", "min4", "max4", "min5", "max5", "min6", "max6", "min7", "max7", "min8", "max8")],#, "min9", "max9")], 
             measure.vars = c("min1", "max1", "min2", "max2", "min3", "max3", "min4", "max4", "min5", "max5", "min6", "max6", "min7", "max7", "min8", "max8"#, "min9", "max9"
                              ))
df.m$mM <- factor(df.m$variable, labels=rep(c("min","max"), 8))
droplevels(df.m)

gp2.all <- ggplot(df.m, aes(y=value, x=variable, fill=mM)) + geom_boxplot() + scale_y_log10(limits=c(1, 1000)) + ylab("") + xlab("Peak limits") + coord_flip()

peak.limits <- aggregate(value ~ variable, data=df.m, FUN = summary)

# Plot limits 
gp2.all + geom_hline(data=data.frame(yint=peak.limits[, 2][, "1st Qu."]), aes(yintercept=yint), col="gray") + geom_hline(data=data.frame(yint=peak.limits[, 2][, "3rd Qu."]), aes(yintercept=yint), col="gray") + theme_bw()


# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 1.a/ Define the limits using 1st and 3rd quartiles ----
# 
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
pl.all <- peak.limits[, 2]
peak.lim.all <- list(
  peak1=c(pl.all[1, "Min."], mean(c(pl.all[2, "3rd Qu."], pl.all[3,"1st Qu."]))),
  peak2=c(mean(c(pl.all[2, "3rd Qu."], pl.all[3,"1st Qu."])), 
          mean(c(pl.all[4, "3rd Qu."], pl.all[5,"1st Qu."]))),
  peak3=c(mean(c(pl.all[4, "3rd Qu."], pl.all[5,"1st Qu."])), 
          mean(c(pl.all[6, "3rd Qu."], pl.all[7,"1st Qu."]))),
  peak4=c(mean(c(pl.all[6, "3rd Qu."], pl.all[7,"1st Qu."])),
          mean(c(pl.all[8, "3rd Qu."], pl.all[9,"1st Qu."]))),
  peak5=c(mean(c(pl.all[8, "3rd Qu."], pl.all[9,"1st Qu."])),
          mean(c(pl.all[10, "3rd Qu."], pl.all[11,"1st Qu."]))),
  peak6=c(mean(c(pl.all[10, "3rd Qu."], pl.all[11,"1st Qu."])),
          mean(c(pl.all[12, "3rd Qu."], pl.all[13,"1st Qu."]))),
  peak7=c(mean(c(pl.L30[12, "3rd Qu."], pl.L30[13,"1st Qu."])), 
          mean(c(pl.L30[14, "3rd Qu."], pl.L30[15,"1st Qu."]))),
  peak8=c(mean(c(pl.L30[14, "3rd Qu."], pl.L30[15,"1st Qu."])),
          pl.L30[16, "Max."])
#  peak7=c(mean(c(pl.all[12, "3rd Qu."], pl.all[13,"1st Qu."])),
#          mean(c(pl.all[14, "3rd Qu."], pl.all[15,"1st Qu."]))),
#  peak8=c(mean(c(pl.all[14, "3rd Qu."], pl.all[15,"1st Qu."])), pl.all[17, "Max."]),
#  peak9=c(pl.all[17, "Max."], pl.all[18, "Max."])
)
pl1.all <- melt(as.data.frame(peak.lim.all))
gp2.all + geom_hline(data=pl1.all, aes(yintercept=value), col="blue") + theme_bw()


# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 1.b/ Plot global peak limits ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
pdf("./figures/estimatedPeakLimits_boxplots.pdf", 7, 7)
gp2.all + geom_hline(data=pl1.all, aes(yintercept=value), col="blue") + theme_bw() 
dev.off()
system("open ./figures/estimatedPeakLimits_boxplots.pdf")

system("open ../figures/detected_and_estimatedPeakLimits.pdf")
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 2/ use these limits as filters in filter() ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
rectGate1.all <- rectangleGate(filterId = "peak1", "DAPI"=peak.lim.all$peak1)
#rectGate1.1 <- rectangleGate(filterId = "peak1", "FL1-A"=peak.lim1.manual$peak1)
#rectGate1.2 <- rectangleGate(filterId = "peak1", "FL1-A"=peak.lim2$peak1)

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 3/ use these limits to extract peak % ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# 
fres00 <- filter(dat[1:622, ], filter=rectGate1.all)
cycleProportion.all <- toTable(summary(fres00))
for(i in c("peak2", "peak3", "peak4", "peak5", "peak6", "peak7", "peak8")) {
  rectGate.i <- rectangleGate(filterId = i, "DAPI"=peak.lim.all[[i]])
  fres00 <- filter(dat[1:622, ], filter=rectGate.i)
  cycleProportion.all <- rbind(cycleProportion.all, toTable(summary(fres00)))
}

cycleProportion.all$ploidy <- factor(cycleProportion.all$population, labels=c("debris", "x2C", "x4C", "x8C", "x16C", "x32C", "x64C", "x128C"))
cycleProportion.all.wide <- dcast(cycleProportion.all[, c("sample", "ploidy", "p")], formula=sample~ploidy)
cycleProportion.all.wide[is.na(cycleProportion.all.wide$x128C), "x128C"] <- 0
#cycleProportion.all.wide[is.na(cycleProportion.all.wide$x128C), "x256C"] <- 0
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Calculate new cycle value (endoreduplication factor) ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# The two periods of measurement treated conjointly ----
dfCVall <- within(data = cycleProportion.all.wide, {
  cycleValue <- (0*x2C + 1*x4C + 2*x8C + 3*x16C + 4*x32C + 5*x64C + 6*x128C)/(x2C + x4C + x8C + x16C + x32C + x64C + x128C)
})
dfCVall$cycleValue
dfCVall <- dfCVall %>%
  left_join(ids, by=c("sample"="fileName")) %>%
  left_join(idPots)

dfCVall$tissueType.ord <- factor(dfCVall$tissueType, 
                             levels = c("f6", "f8", "F8", "f30"), 
                             labels = c("Seedling_Leaf5","Leaf_8", "Leaf_8", "Leaf_30"))

dfCVall[dfCVall$tissueType.ord=="Seedling_Leaf5", "watering"] <- "WW"
dfCVall$watering <- factor(dfCVall$watering, levels = c("WW", "WD"))

CV.mean <- subset(dfCVall[1:622,]) %>%
  select(nameGen, tissueType.ord, watering, cycleValue) %>%
  group_by(nameGen, tissueType.ord, watering) %>%
  summarize(CVmean = mean(cycleValue))

CV.mean.wide <- dcast(CV.mean, formula=nameGen~tissueType.ord + watering, mean) 

dfCVall$nameGen.OrderedSeedling <- as.factor(dfCVall$nameGen)

dfCVall$nameGen.OrderedSeedling <- factor(dfCVall$nameGen.OrderedSeedling,
                                          levels =levels(dfCVall$nameGen.OrderedSeedling)[order(CV.mean.wide$Seedling_Leaf5_WW)])

dfCVall$nameGen.OrderedL8_WW <- as.factor(dfCVall$nameGen)
dfCVall$nameGen.OrderedL8_WW <- factor(dfCVall$nameGen.OrderedL8_WW,
                                          levels =levels(dfCVall$nameGen.OrderedL8_WW)[order(CV.mean.wide$Leaf_8_WW)])

CV.mean.wide$clust <- as.factor(kmeans(CV.mean.wide[, "Seedling_Leaf5_WW"], centers=5)$cluster)
dfCVall <- dfCVall %>%
  left_join(CV.mean.wide)


pdf("./figures/cycleValue_30_genotypes.pdf", 12, 8)
ggplot(data=dfCVall[1:622,], aes(y=cycleValue, x=tissueType.ord, colour=watering)) + 
  geom_boxplot() + 
  ylab("Endoreduplication factor (EF)") + xlab("Leaf sample") +
  facet_wrap(.~nameGen) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))

ggplot(data=subset(dfCVall[1:622,], tissueType.ord%in%c("Seedling_Leaf5")), 
       aes(y=cycleValue, x=nameGen.OrderedSeedling, fill = clust)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(size = 1) +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))

ggplot(data=subset(dfCVall[1:622,], tissueType.ord%in%c("Leaf_8")), 
       aes(y=cycleValue, x=nameGen.OrderedL8_WW, fill=watering)) +
  geom_boxplot(outlier.alpha = 0) +
  xlab("Accessions (ordered by L8_WW CV)") +
  ylab("Cycle value") +
  scale_fill_manual(values=c("#386FA4", "#84D2F6")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))

ggplot(data=subset(dfCVall[1:622,], tissueType.ord%in%c("Leaf_30")), 
       aes(y=cycleValue, x=nameGen.OrderedL8_WW, fill=watering)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("#386FA4", "#84D2F6")) +theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))

dev.off()
system("open ./figures/cycleValue_30_genotypes.pdf")

subset(dfCVall[1:622,], x256C > 0)

library(ggrepel)
gp.corr <- ggplot(data=CV.mean.wide, aes(x = Seedling_Leaf5_WW, y = Leaf_8_WW)) +
  geom_point() + geom_smooth(method = lm, se=F) + geom_text_repel(aes(label=nameGen)) +
  theme_bw() #+ geom_abline(slope = 1, intercept = 0)

pdf("./figures/cycleValue_30_genotypes_correlations.pdf", 8, 7)
gp.corr + geom_abline(slope = 1, intercept = 0) 
gp.corr + aes(x = Seedling_Leaf5_WW, y = Leaf_8_WD)
gp.corr + aes(x = Leaf_8_WW, y = Leaf_8_WD) + geom_abline(slope = 1, intercept = 0)
gp.corr + aes(x = Leaf_30_WW, y = Leaf_30_WD) + geom_abline(slope = 1, intercept = 0)
gp.corr + aes(x = Seedling_Leaf5_WW, y = Leaf_30_WW) 
gp.corr + aes(x = Seedling_Leaf5_WW, y = Leaf_30_WD) 
dev.off()
system("open ./figures/cycleValue_30_genotypes_correlations.pdf")


mean(subset(dfCVall[1:622,], tissueType.ord%in%c("Seedling_Leaf5"))$cycleValue)
