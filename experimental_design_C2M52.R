# *************************************************
# Experimental design ----
# PHENOPSIS C2M52
# *************************************************

geno <- d.geno[order(d.geno$latitude), 1:6]
geno$numGen <- 1:30
names(geno)[1] <- "ecotypeID"
geno

set.seed(3) # Set (to 3, here) for replication purpose.
idPot <- sample(1:504)

idPot <- 1:504

# Sample 12 genotypes for 9 replicates per watering condition ----
set.seed(3) # Set (to 3, here) for replication purpose.
geno9rep <- sample(1:30, 12, replace = F)
# Remaining genotypes with 8 replicates per watering condition ----
geno8rep <- geno$numGen[!(geno$numGen %in% geno9rep)]

# Create data frame with 504 pots (randomly sampled), replicated genotypes, and watering conditions (WW = well-watered, WD=water deficit) ----
pots <- data.frame(idPot = idPot,
                   numGen = rep(c(rep(geno8rep, 8), rep(geno9rep, 9)), 2),
                   watering = factor(c(rep("WW", 504/2), rep("WD", 504/2)), levels = c("WW", "WD"))
                   )

pots <- pots %>%
  left_join(geno[, c(7,1,2)]) %>%
  group_by(idPot) %>%
  arrange(idPot) %>%
  dplyr::select(idPot, numGen, ecotypeID, name, watering)
pots

# Verify replicate number per genotype per condition ----
pots %>%
  group_by(numGen, ecotypeID, watering) %>%
  summarise(pot_n = n())

# Save data frame for Excel handling and sowing procedure ----

write.csv2(pots, file="./data/pots_C2M52.csv", row.names = F)

