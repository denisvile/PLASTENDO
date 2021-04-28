# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Cytometry of Arabidopsis leaf samples
# Project: [ARABREED]
# Subproject: plasticity of endopolyploidy [PLASTENDO] 2021
# Benoît Berthet Master 2 
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€

### PLOTS AND ANALYSES

gp.corr <- ggplot(data=CV.mean.L30.wide, aes(x = Leaf_30_WW, y = Leaf_30_WD)) +
  geom_point() + geom_smooth(method = lm, se=F) 
gp.corr + geom_abline(slope = 1, intercept = 0)

dfCV.Lsdlg_all <- dfCV.Lsdlg[, c("nameGen", "watering", "tissueType.ord", "idPot", "cycleValue")] %>%
  left_join(subset(dfCVall[1:622, c("nameGen", "watering", "tissueType.ord", "idPot", "cycleValue")], tissueType.ord=="Leaf_30"), by="idPot")

pdf("./figures/cycleValue_correlations.Lsdlg_vs_all_bw_fac_0.5.pdf", 8, 7)
ggplot(dfCV.Lsdlg_all, aes(x=cycleValue.x, y=cycleValue.y)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  xlab("Cycle value with gating on Lsdlg samples") +
  ylab("Cycle value with gating on all samples") +
  theme(text = element_text(size=16))
dev.off()
system("open ./figures/cycleValue_correlations.Lsdlg_vs_all_bw_fac_0.5.pdf")

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Boxplot of cycle value for seelings leaves ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
pdf("./figures/cycleValue_30_genotypes.Lsdlg.pdf", 12, 8)
gp.CV.sdlg <- ggplot(data=subset(dfCV.Lsdlg, tissueType.ord%in%c("Seedling_Leaf5")), 
       aes(y=cycleValue, x=nameGen.OrderedLsdlg_WW, fill=clust)) +
  geom_boxplot(outlier.alpha = 0) +
  xlab("Accessions (ordered by Lsdlg_WW CV)") +
  ylab("Cycle value") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))
dev.off()
system("open ./figures/cycleValue_30_genotypes.Lsdlg.pdf")

# Cycle values of the three harvests merged together -----
CV.mean.all.3datasets <- CV.mean.L30.wide %>%
  left_join(CV.mean.L8.wide) %>%
  left_join(CV.mean.Lsdlg.wide)

gp.corr.3datasets <- ggplot(data=CV.mean.all.3datasets, aes(x = Seedling_Leaf5_WW, y = Leaf_8_WW)) +
  geom_point() + geom_smooth(method = lm, se=F)

pdf("./figures/cycleValue_30_genotypes_correlations.3datasets.pdf", 8, 7)
gp.corr.3datasets + geom_abline(slope = 1, intercept = 0) 
gp.corr.3datasets + aes(x = Seedling_Leaf5_WW, y = Leaf_8_WD)
gp.corr.3datasets + aes(x = Leaf_8_WW, y = Leaf_8_WD) + geom_abline(slope = 1, intercept = 0)
gp.corr.3datasets + aes(x = Leaf_30_WW, y = Leaf_30_WD) + geom_abline(slope = 1, intercept = 0)
gp.corr.3datasets + aes(x = Seedling_Leaf5_WW, y = Leaf_30_WW) 
gp.corr.3datasets + aes(x = Seedling_Leaf5_WW, y = Leaf_30_WD) 
dev.off()
system("open ./figures/cycleValue_30_genotypes_correlations.3datasets.pdf")

# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Flow cytometry data for contrasted accessions ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
unique(ids$idGen)
An1_f6 <- ids %>% dplyr::filter(idGen=="an1" & tissueType=="f6")

IP_Hum_2_f6 <- ids %>% dplyr::filter(idGen=="iphum2" & tissueType=="f6")

tal03_f6 <- ids %>% dplyr::filter(idGen=="tal03" & tissueType=="f6")

IP_Coa_0_f6 <- ids %>% dplyr::filter(idGen=="ipcoa0" & tissueType=="f6")

rowN_An1 <- which(rownames(phenoData(d.Lsdlg)) %in% An1_f6$fileName)
rowN_IP_Hum_2 <- which(rownames(phenoData(d.Lsdlg)) %in% IP_Hum_2_f6$fileName)
rowN_tal03 <- which(rownames(phenoData(d.Lsdlg)) %in% tal03_f6$fileName)
rowN_IP_Coa_0 <- which(rownames(phenoData(d.Lsdlg)) %in% IP_Coa_0_f6$fileName)

dat_An1_f6 <- transform(d.Lsdlg[rowN_An1, ], "lgDAPI"=log10(`DAPI`))
dat_IP_Hum_2_f6 <- transform(d.Lsdlg[rowN_IP_Hum_2, ], "lgDAPI"=log10(`DAPI`))
dat_tal03_f6 <- transform(d.Lsdlg[rowN_tal03, ], "lgDAPI"=log10(`DAPI`))
dat_IP_Coa_0_f6 <- transform(d.Lsdlg[rowN_IP_Coa_0, ], "lgDAPI"=log10(`DAPI`))

ddf_An1_f6 <- data.frame(
  lgDAPI= c(as.data.frame(exprs(dat_An1_f6[1, "lgDAPI"][[1]]))$lgDAPI,
             as.data.frame(exprs(dat_An1_f6[2, "lgDAPI"][[1]]))$lgDAPI,
             as.data.frame(exprs(dat_An1_f6[3, "lgDAPI"][[1]]))$lgDAPI,
             as.data.frame(exprs(dat_An1_f6[4, "lgDAPI"][[1]]))$lgDAPI),
  idPot = c(rep(rowN_An1[1], dim(as.data.frame(exprs(dat_An1_f6[1, "lgDAPI"][[1]])))[1]),
            rep(rowN_An1[2], dim(as.data.frame(exprs(dat_An1_f6[2, "lgDAPI"][[1]])))[1]),
            rep(rowN_An1[3], dim(as.data.frame(exprs(dat_An1_f6[3, "lgDAPI"][[1]])))[1]),
            rep(rowN_An1[4], dim(as.data.frame(exprs(dat_An1_f6[4, "lgDAPI"][[1]])))[1])
            ))

ddf_IP_Hum_2_f6 <- data.frame(lgDAPI= c(as.data.frame(exprs(dat_IP_Hum_2_f6[1, "lgDAPI"][[1]]))$lgDAPI, as.data.frame(exprs(dat_IP_Hum_2_f6[2, "lgDAPI"][[1]]))$lgDAPI, as.data.frame(exprs(dat_IP_Hum_2_f6[3, "lgDAPI"][[1]]))$lgDAPI, as.data.frame(exprs(dat_IP_Hum_2_f6[4, "lgDAPI"][[1]]))$lgDAPI),
                              idPot = c(rep(rowN_IP_Hum_2[1], dim(as.data.frame(exprs(dat_IP_Hum_2_f6[1, "lgDAPI"][[1]])))[1]),
                                        rep(rowN_IP_Hum_2[2], dim(as.data.frame(exprs(dat_IP_Hum_2_f6[2, "lgDAPI"][[1]])))[1]),
                                        rep(rowN_IP_Hum_2[3], dim(as.data.frame(exprs(dat_IP_Hum_2_f6[3, "lgDAPI"][[1]])))[1]),
                                        rep(rowN_IP_Hum_2[4], dim(as.data.frame(exprs(dat_IP_Hum_2_f6[4, "lgDAPI"][[1]])))[1])
)
)

ddf_Tal_03_f6 <- data.frame(
  lgDAPI= c(
  as.data.frame(exprs(dat_tal03_f6[1, "lgDAPI"][[1]]))$lgDAPI, 
  as.data.frame(exprs(dat_tal03_f6[2, "lgDAPI"][[1]]))$lgDAPI, 
  as.data.frame(exprs(dat_tal03_f6[3, "lgDAPI"][[1]]))$lgDAPI, 
  as.data.frame(exprs(dat_tal03_f6[4, "lgDAPI"][[1]]))$lgDAPI, 
  as.data.frame(exprs(dat_tal03_f6[5, "lgDAPI"][[1]]))$lgDAPI),
  idPot = c(rep(rowN_tal03[1], dim(as.data.frame(exprs(dat_tal03_f6[1, "lgDAPI"][[1]])))[1]),
                                        rep(rowN_tal03[2], dim(as.data.frame(exprs(dat_tal03_f6[2, "lgDAPI"][[1]])))[1]),
                                        rep(rowN_tal03[3], dim(as.data.frame(exprs(dat_tal03_f6[3, "lgDAPI"][[1]])))[1]),
                                        rep(rowN_tal03[4], dim(as.data.frame(exprs(dat_tal03_f6[4, "lgDAPI"][[1]])))[1]),
            rep(rowN_tal03[5], dim(as.data.frame(exprs(dat_tal03_f6[5, "lgDAPI"][[1]])))[1]))
)

ddf_IP_Coa_0_f6 <- data.frame(
  lgDAPI= c(
    as.data.frame(exprs(dat_IP_Coa_0_f6[1, "lgDAPI"][[1]]))$lgDAPI, 
    as.data.frame(exprs(dat_IP_Coa_0_f6[2, "lgDAPI"][[1]]))$lgDAPI, 
    as.data.frame(exprs(dat_IP_Coa_0_f6[3, "lgDAPI"][[1]]))$lgDAPI, 
    as.data.frame(exprs(dat_IP_Coa_0_f6[4, "lgDAPI"][[1]]))$lgDAPI, 
    as.data.frame(exprs(dat_IP_Coa_0_f6[5, "lgDAPI"][[1]]))$lgDAPI),
  idPot = c(rep(rowN_IP_Coa_0[1], dim(as.data.frame(exprs(dat_IP_Coa_0_f6[1, "lgDAPI"][[1]])))[1]),
            rep(rowN_IP_Coa_0[2], dim(as.data.frame(exprs(dat_IP_Coa_0_f6[2, "lgDAPI"][[1]])))[1]),
            rep(rowN_IP_Coa_0[3], dim(as.data.frame(exprs(dat_IP_Coa_0_f6[3, "lgDAPI"][[1]])))[1]),
            rep(rowN_IP_Coa_0[4], dim(as.data.frame(exprs(dat_IP_Coa_0_f6[4, "lgDAPI"][[1]])))[1]),
            rep(rowN_IP_Coa_0[5], dim(as.data.frame(exprs(dat_IP_Coa_0_f6[5, "lgDAPI"][[1]])))[1]))
)

ddf_An1_f6$idGen <- "An-1"
ddf_IP_Hum_2_f6$idGen <- "IP-Hum-2"
ddf_Tal_03_f6$idGen <- "Tal-03"
ddf_IP_Coa_0_f6$idGen <- "IP-Coa-0"

ddf_genoSelect <- rbind(ddf_An1_f6, 
                        ddf_IP_Hum_2_f6, 
                        ddf_Tal_03_f6,
                        ddf_IP_Coa_0_f6)

pdf("./figures/Distribution_lgDAPI_selected_geno.pdf", 10, 8)
ggplot(subset(ddf_genoSelect), 
       aes(x=lgDAPI, color=as.factor(idPot))) + 
  geom_histogram(aes(y=..density..), bins = 150) + 
  geom_density(aes(color=as.factor(idPot))) +
  facet_wrap(.~idGen, nrow=2, ncol=2) +
  theme(legend.position = "none")
dev.off()
system(paste("open","./figures/Distribution_lgDAPI_selected_geno.pdf"))

pdf("./figures/CV_selected_geno.pdf", 5, 5)
gp.CV.sdlg %+% subset(dfCV.Lsdlg, tissueType.ord%in%c("Seedling_Leaf5") &
                        nameGen.OrderedLsdlg_WW %in% c("An-1", "IP-Coa-0", "IP-Hum-2", "TAL 03"))
dev.off()
system(paste("open","./figures/CV_selected_geno.pdf"))

