

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
