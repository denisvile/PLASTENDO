#
# C3M42 
# Leaf area and mass of leaves and rosettes
#

d.mass <- read.xls("./data/mass_SLA_C3M42.xlsx")
d.area <- read.xls("./data/Surface_feuille_C3M42.xlsx")

d.SLA <- d.mass %>% 
  left_join(d.area) %>%
  left_join(idPots, by=c("idPot"="idPot"))


# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Diagnose plots for identifcation of outliers (measurement errors...) ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Area vs. dry weight ----
gg.diagnose <- ggplot(d.SLA, aes(x=leaf_area_mm2, y=dry_weight_mg, color=stage)) + geom_point() +
  scale_x_log10() + scale_y_log10()
gg.diagnose
linmod_y_x <- lm(data=d.SLA, log10(dry_weight_mg) ~ log10(leaf_area_mm2))
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
gg.diagnose + aes(x=leaf_area_mm2, y=fresh_weight_mg)
linmod_y_x <- lm(data=d.SLA, log10(fresh_weight_mg) ~ log10(leaf_area_mm2))
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

gp.diag2 <-gg.diagnose + aes(x=leaf_area_mm2, y=fresh_weight_mg) + 
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

d.SLA.m <- d.SLA %>% 
  mutate(SLA = leaf_area_mm2/dry_weight_mg) %>%
  dplyr::filter(!(idPot %in% as.numeric(idPot.exclude.SLA$idPot[1:10]))) %>%
  group_by(nameGen, stage, watering.x) %>%
  dplyr::summarise(SLA = mean(SLA, na.omit=T)) %>%
  left_join(CV.mean.all.3datasets, by=c("nameGen" = "nameGen"))

ggplot(subset(d.SLA.m, stage=="L9" & watering.x=="WW"), aes(x=Leaf_8_WW, y=SLA)) + geom_point() 
ggplot(subset(d.SLA.m), aes(x=Seedling_Leaf5_WW, y=SLA)) + geom_point() 
       