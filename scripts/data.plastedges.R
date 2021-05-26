

d.plastedges <- read.csv("./data/mean_data_plastedges.csv")[,-1]
names(d.plastedges)

d.plastedges <- read.xls("./data/mean_data_plastedges.xlsx")
names(d.plastedges)

d.geno <- d.plastedges %>%
  select(Genotype, idGen, latitude, ER, BTDAY) %>%
  unique()

dim(d.geno)
d.geno

write.csv(d.geno, "./data/data.genotypes.csv")

ggplot(d.geno, aes(x = latitude, y = ER)) + 
  geom_point()

length(d.SLA.m.wide$nameGen[d.SLA.m.wide$nameGen %in% d.geno$idGen])

d.geno1 <- d.SLA.m.wide %>%
  left_join(d.geno, by=c("nameGen"="idGen"))

gp.rel.latitude.EndoSeedlings <- ggplot(subset(d.geno1, stage=="L9" & !(nameGen %in% c("IP-Vis-0", "IP-Hum-2", "Cvi-0"))), 
                                        aes(x=latitude, y =  Seedling_Leaf5_WW)) +
  geom_point() + geom_text_repel(aes(label=nameGen)) + geom_smooth(method="lm", formula=y~x, se=F, col="black") +
  geom_point(data=subset(d.geno1, stage=="L9" & (nameGen %in% c("IP-Vis-0", "IP-Hum-2", "Cvi-0"))), col="grey") +
  geom_text_repel(data=subset(d.geno1, stage=="L9" & (nameGen %in% c("IP-Vis-0", "IP-Hum-2", "Cvi-0"))), col="grey", aes(label=nameGen)) +
  xlab("Latitude") + ylab("Endopolyploidy factor (seedlings)") + xlim(10,80) +
  myTheme

pdf("./figures/latitude.Endo.pdf", 18, 5)
ggarrange(
  gp.rel.latitude.EndoSeedlings + annotate(geom="text", x=15, y=1, label=c("r = 0.50**"), size=5),
  gp.rel.latitude.EndoSeedlings + aes(y = Leaf_CV_WW) + ylab("Endopolyploidy factor (leaf #8)"),
  gp.rel.latitude.EndoSeedlings %+% subset(d.geno1, stage=="L30" & !(nameGen %in% c("IP-Vis-0", "IP-Hum-2", "Cvi-0"))) + aes(y = Leaf_CV_WW) + ylab("Endopolyploidy factor (leaf #30)"),
  ncol=3
)
dev.off()
system("open ./figures/latitude.Endo.pdf")

summary(lm(data=subset(d.geno1, stage=="L9" & !(nameGen %in% c("IP-Vis-0", "IP-Hum-2",  "Cvi-0"))),
   Seedling_Leaf5_WW ~ latitude))
