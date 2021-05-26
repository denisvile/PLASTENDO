
# Individual leaf growth

# Relative growth rate

str(ilg1)

ilg1.try <- subset(ilg1, idPot==5 & leafType=="F8")

ilg1.try.RER <- ilg1.try %>%
  select(idPot, leafType, day, Area.total.mm2) %>%
  mutate(rn = rownames(.)) %>%
  tidyr::pivot_wider(names_from=rn, values_from=c(Area.total.mm2, day)) %>%
  mutate(RER.1_2 = (log(Area.total.mm2_2) - log(Area.total.mm2_1)/(day_2 - day_1)),
         )

ilg1.try.RER$RER.1_2


library(imageData)
library(growthPheno)

?intervalGRaverage

ilg1.splines <- splitSplines(ilg1, response="Area.total.mm2", x="day", 
             INDICES = "idPot", 
             df = 4, deriv=1, suffices.deriv="AGRdv", RGR="RGRdv")

data(exampleData)
longi.dat <- splitSplines(longi.dat, response="Area", x="xDays", 
                          INDICES = "Snapshot.ID.Tag", 
                          df = 4, deriv=1, suffices.deriv="AGRdv", RGR="RGRdv")

Area.smooth.GR <- intervalGRaverage("Area.smooth", which.rates = c("AGR","RGR"), 
                                    suffices.rates = c("AGRdv","RGRdv"), 
                                    start.time = 31, end.time = 35, 
                                    suffix.interval = "31to35",
                                    data = longi.dat)

plt <- plotLongitudinal(data = longi.dat, response = "Area.smooth", x.title = "DAP",  
                        y.title = "Area.smooth", x="xDays+35.42857143", printPlot=FALSE)
plt <- plt + ggplot2::geom_vline(xintercept=29, linetype="longdash", size=1) +
  ggplot2::scale_x_continuous(breaks=seq(28, 42, by=2)) + 
  ggplot2::scale_y_continuous(limits=c(0,750))
print(plt)

plotLongitudinal(data = longi.dat, response = "Area.smooth", x.title = "DAP",  
                 y.title = "Area.smooth", x="xDays+35.42857143", 
                 ggplotFuncs = list(ggplot2::geom_vline(xintercept=29, linetype="longdash", 
                                                        size=1), 
                                    ggplot2::scale_x_continuous(breaks=seq(28, 42, by=2)), 
                                    ggplot2::scale_y_continuous(limits=c(0,750))))


