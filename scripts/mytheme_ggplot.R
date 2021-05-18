# ggplot theme ----

textSize <- 14

myTheme <- theme(
  panel.background=element_rect(fill="transparent", color= "black", size = 1),  
  plot.background = element_rect(fill = "transparent",colour = NA),
  axis.line=element_blank(), # line(color="black", size = 1), 
  axis.title=element_text(size=textSize),
  axis.text.y=element_text(size=textSize, colour="black"), 
  axis.text.x=element_text(size=textSize, colour="black", angle=0, hjust=0.5),
  panel.grid.minor = element_blank(), 
  panel.grid.major = element_blank(), # element_line(colour="grey90", size=0.2), 
  legend.position="right", 
  legend.text=element_text(size=textSize),
  legend.title=element_text(size=textSize),
  strip.background=element_rect(fill="transparent", color=NA), 
  strip.text = element_text(color="black", size=textSize)
)
