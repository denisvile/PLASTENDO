# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Map of accessions location ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€


tempcol <- colorRampPalette(c("purple", "blue", "skyblue", "green", "lightgreen", "yellow", "orange", "red", "darkred"))

# Countries
some.eu.countries <- c("Portugal", "Spain", "France", "Morocco", "Algeria", "Sweden", "Germany", "Austria", "Nederland", "Cape Verde", "Italy", "")

some.eu.countries <-c("Austria","Belgium","Bulgaria","Croatia","Cyprus",
                      "Czech Rep.","Denmark","Estonia","Finland","France",
                      "Germany","Greece","Hungary","Ireland","Italy","Latvia",
                      "Lithuania","Luxembourg","Malta","Netherlands","Poland",
                      "Portugal","Romania","Slovakia","Slovenia","Spain",
                      "Sweden","United Kingdom")
# Retrieve the map data
some.eu.maps <- map_data("world") #, region = some.eu.countries)


region.lab.data <- some.eu.maps %>%
  group_by(region) %>%
  summarise(long = mean(long), lat = mean(lat))

gp.map <- ggplot(some.eu.maps, aes(x = long, y = lat)) +
  geom_polygon(aes(group = group), fill = "grey", col="grey50") +
  geom_point(data = d.geno,
             aes(x = longitude, y = latitude), shape=21, size=5, fill = "red") +
  geom_text_repel(data = d.geno,
                  aes(x = longitude, y = latitude, label=name), segment.size = 0, segment.alpha=0, size=4, colour="black") +
#  geom_text(aes(label = region), data = region.lab.data,  size = 5, hjust = 0.5)+
  #scale_fill_gradient("Altitude", low = "#EBB1A6", high = "red") +
  coord_cartesian(xlim = c(-30, 30), 
                  ylim = c(10,70)) + 
  xlab("Longitude") + ylab("Latitude") +
  theme(legend.position = c(0.95,0.05),
        legend.justification = c(1,0),
        axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank(),
        panel.background = element_blank())
gp.map

pdf("./figures/map.pdf", 8,7)
gp.map
dev.off()
system("open ./figures/map.pdf")
