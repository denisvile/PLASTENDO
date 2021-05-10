#
# C3M42 
# Leaf area and mass of leaves and rosettes
#

d.mass <- read.xls("./data/SLA_C3M42_v3.xlsx")
d.area <- read.xls("./data/Leaf_area_v3.xlsx")

d.SLA <- d.mass %>%
  select(-name, -ecotypeID, -watering) %>%
  left_join(d.area) %>%
  left_join(idPots, by=c("idPot"="idPot")) %>%
  mutate(watering=factor(watering, levels=c("WW", "WD")))
head(d.SLA)


# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Diagnose plots for identifcation of outliers (measurement errors...) ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Area vs. dry weight ----
gg.diagnose <- ggplot(d.SLA, aes(x=leaf.area.mm2, y=dry_weight_mg, color=stage)) + geom_point() +
  scale_x_log10() + scale_y_log10()
gg.diagnose
linmod_y_x <- lm(data=d.SLA, log10(dry_weight_mg) ~ log10(leaf.area.mm2))
diagnose_plots <-lindia::gg_diagnose(linmod_y_x, plot.all=F)
diagnose_plots$cooksd$
lm_matrix <- fortify(linmod_y_x)
lm_matrix[, "rowname"] <- 1:nrow(lm_matrix)
cooksd = lm_matrix[, ".cooksd"]
n = nrow(lm_matrix)
idPot_above_threshold <- subset(diagnose_plots$cooksd$data, .cooksd > 4/n ) # 4/n is the convention threshold
# na.omit(d.SLA[, c("leaf_area_mm2", "dry_weight_mg", "idPot")])[, ]
idPot.exclude.SLA <- idPot_above_threshold %>%
  mutate(idPot=rownames(.)) %>%
  arrange(desc(.cooksd),)
idPot.exclude.SLA

gp.diag1 <- gg.diagnose + geom_text_repel(data=d.SLA[as.numeric(rownames((idPot_above_threshold))),], aes(label=idPot))

# Area vs. fresh weight ----
gg.diagnose + aes(x=leaf.area.mm2, y=fresh_weight_mg)
linmod_y_x <- lm(data=d.SLA, log10(fresh_weight_mg) ~ log10(leaf.area.mm2))
diagnose_plots <-lindia::gg_diagnose(linmod_y_x, plot.all=F)
diagnose_plots$cooksd
lm_matrix <- fortify(linmod_y_x)
lm_matrix[, "rowname"] <- 1:nrow(lm_matrix)
cooksd = lm_matrix[, ".cooksd"]
n = nrow(lm_matrix)
idPot_above_threshold <- subset(diagnose_plots$cooksd$data, .cooksd > 4/n ) 
idPot_above_threshold %>%
  mutate(idPot=rownames(.)) %>%
  arrange(desc(.cooksd),)

gp.diag2 <-gg.diagnose + aes(x=leaf.area.mm2, y=fresh_weight_mg) + 
  geom_text_repel(data=d.SLA[as.numeric(rownames((idPot_above_threshold))),], aes(label=idPot))


# Fresh vs. dry weight ----
gg.diagnose + aes(x=fresh_weight_mg, y=dry_weight_mg)
linmod_y_x <- lm(data=d.SLA, log10(dry_weight_mg) ~ log10(fresh_weight_mg))
diagnose_plots <-lindia::gg_diagnose(linmod_y_x, plot.all=F)
diagnose_plots$cooksd
lm_matrix <- fortify(linmod_y_x)
lm_matrix[, "rowname"] <- 1:nrow(lm_matrix)
cooksd = lm_matrix[, ".cooksd"]
n = nrow(lm_matrix)
idPot_above_threshold <- subset(diagnose_plots$cooksd$data, .cooksd > 4/n ) 
idPot_above_threshold %>%
  mutate(idPot=rownames(.)) %>%
  arrange(desc(.cooksd),)

gp.diag3 <-gg.diagnose + aes(x=fresh_weight_mg, y=dry_weight_mg) + 
  geom_text_repel(data=d.SLA[as.numeric(rownames((idPot_above_threshold))),], aes(label=idPot))

# Rehydrated wt vs. fresh weight ----
gg.diagnose + aes(x=fresh_weight_mg, y=rehydrated_weight_mg)
linmod_y_x <- lm(data=d.SLA, rehydrated_weight_mg ~ fresh_weight_mg)
diagnose_plots <-lindia::gg_diagnose(linmod_y_x, plot.all=F)
diagnose_plots$cooksd
lm_matrix <- fortify(linmod_y_x)
lm_matrix[, "rowname"] <- 1:nrow(lm_matrix)
cooksd = lm_matrix[, ".cooksd"]
n = nrow(lm_matrix)
idPot_above_threshold <- subset(diagnose_plots$cooksd$data, .cooksd > 4/n ) 
idPot_above_threshold %>%
  mutate(idPot=rownames(.)) %>%
  arrange(desc(.cooksd),)

gp.diag4 <- ggplot(d.SLA, aes(x=fresh_weight_mg, y=rehydrated_weight_mg, color=stage)) + geom_point() + 
  geom_text_repel(data=d.SLA[as.numeric(rownames((idPot_above_threshold))),], aes(label=idPot))

pdf("./figures/diagnose_plots_area_mass.pdf", 7, 5)
gp.diag1
gp.diag2
gp.diag3
gp.diag4
dev.off()
system("open ./figures/diagnose_plots_area_mass.pdf")


# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# SLA calculation with exclusion of above identified outliers
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
d.SLA <- d.SLA %>% 
  mutate(SLA = (leaf.area.mm2)/(dry_weight_mg))

gp.SLA <- ggplot(subset(d.SLA, idPot != 312 & stage %in%c("L30", "L9")), 
       aes(y=SLA, x=interaction(watering, nameGen), fill=watering)) +
  geom_boxplot() +
#  ylab(expression(paste("SLA"))) +
  xlab("") +
#  scale_y_log10() +
  facet_wrap(.~stage, nrow=2) +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12, colour="black", angle=90, hjust=1),
        legend.position = "none")

pdf("./figures/SLA.pdf", 10, 8)
gp.SLA
gp.SLA + aes(y=leaf.area.mm2) + scale_y_log10()
gp.SLA + aes(y=dry_weight_mg) + scale_y_log10()
dev.off()
system("open ./figures/SLA.pdf")


d.SLA.m.0 <- subset(d.SLA, idPot != 312 & stage %in% c("L9", "L30")) %>%
  group_by(nameGen, watering, stage) %>%
  summarise(SLA = mean(SLA, na.rm=T),
            leaf.area.mm2 = mean(leaf.area.mm2, na.rm=T))


CV.mean.watering.long <- data.frame(nameGen=rep(CV.mean.all.3datasets$nameGen, 4),
           watering=c(rep("WW", length(CV.mean.all.3datasets$Leaf_8_WW)), rep("WD", length(CV.mean.all.3datasets$Leaf_8_WD)),
                      rep("WW", length(CV.mean.all.3datasets$Leaf_30_WW)), rep("WD", length(CV.mean.all.3datasets$Leaf_30_WD))),
           stage=c(rep("L9", length(CV.mean.all.3datasets$Leaf_8_WW)), rep("L9", length(CV.mean.all.3datasets$Leaf_8_WD)),
                      rep("L30", length(CV.mean.all.3datasets$Leaf_30_WW)), rep("L30", length(CV.mean.all.3datasets$Leaf_30_WD))),
           Leaf_CV=c(CV.mean.all.3datasets$Leaf_8_WW, CV.mean.all.3datasets$Leaf_8_WD, 
                     CV.mean.all.3datasets$Leaf_30_WW, CV.mean.all.3datasets$Leaf_30_WD))

d.SLA.m <- subset(d.SLA, idPot != 312 & stage %in% c("L9", "L30")) %>% 
  mutate(SLA = leaf.area.mm2/dry_weight_mg) %>%
  # dplyr::filter(!(idPot %in% as.numeric(idPot.exclude.SLA$idPot[1:10]))) %>%
  group_by(nameGen, stage, watering) %>%
  dplyr::summarise(SLA = mean(SLA, na.omit=T),
                   leaf.area.mm2 = mean(leaf.area.mm2, na.rm=T),
                   leaf.dry.mass = mean(dry_weight_mg, na.rm=T)) %>%
  left_join(CV.mean.watering.long)

ggplot(subset(d.SLA.m, stage=="L9"),
       aes(x=Leaf_CV, y=leaf.area.mm2)) + 
  geom_point(aes(colour=watering)) +
  geom_smooth(method=lm, se=F) + geom_smooth(aes(colour=watering), method=lm, se=F) +
  geom_text_repel(aes(label=nameGen)) +
  facet_wrap(.~stage)

ggplot(subset(d.SLA.m, stage=="L9" & !(nameGen %in% c("IP-Vis-0", "IP-Hum-2", "Kulturen-1", "Com-1"))),
       aes(x=Leaf_CV, y=leaf.area.mm2)) + 
  geom_point(aes(colour=watering)) +
  geom_point(data=subset(d.SLA.m, stage=="L9" & (nameGen %in% c("IP-Vis-0", "IP-Hum-2", "Kulturen-1", "Com-1"))), col="grey") +
  geom_smooth(method=lm, se=F, col="black") + geom_smooth(aes(colour=watering), method=lm, se=F) +
 # geom_text_repel(aes(label=nameGen)) +
  facet_wrap(.~stage) +
  theme_bw()

d.SLA.m %>%
  group_by(nameGen, stage) %>%
  dplyr::filter(stage=="L30" & !(nameGen %in% c("IP-Vis-0", "IP-Hum-2"))) %>%
  summarise(leaf.area.mm2 = mean(leaf.area.mm2, na.rm=T),
            leaf.dry.mass = mean(leaf.dry.mass, na.rm=T),
            Leaf_CV = mean(Leaf_CV, na.rm=T)) %>%
ggplot(., aes(x=Leaf_CV, y=leaf.area.mm2)) + 
  geom_point() +
  geom_smooth(method=lm, se=F) +
  geom_text_repel(aes(label=nameGen)) +
  facet_wrap(.~stage)


ggplot(subset(d.SLA.m), aes(x=Seedling_Leaf5_WW, y=SLA)) + geom_point() 
       