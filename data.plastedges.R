

d.plastedges <- read.csv("./data/mean_data_plastedges.csv")[,-1]
names(d.plastedges)

d.geno <- d.plastedges %>%
  select(Genotype, name, latitude, ER, BTDAY) %>%
  unique()

dim(d.geno)
d.geno

ggplot(d.geno, aes(x = latitude, y = ER)) + 
  geom_point()
