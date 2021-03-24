# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Individual leaf growth data
# Image analyses by M1 Plant Science students (Marie, Mouad, Harkingto)
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€


# Data import ----

ilg <- read.xls("./data/C3M42_leafGrowth.xlsx")
str(ilg)
# day
# day 1 = 2021-02-01

gp.ilg <- ggplot(data = ilg, aes(y = Area.total.mm2, x = day)) +
  geom_point()
gp.ilg

gp.ilg %+% subset(ilg, leafType %in% c("F8")) + facet_wrap(.~ leafType+idGenotype) + aes(colour=factor(idPot))


