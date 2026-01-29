# Daniel Castaneda Mogollon, PhD
# September 4th, 2025
# Purpose: This script was made to geneerate a bar and a gradient legend associated with that bar. Parameters can change to 
# adjust colour, size, font, etc.

library(ggplot2)
library(scales)
df = data.frame(x=seq(0,50,length.out=100),
                y=1, #change axis values if you want the bar to be vertical
                fill=(seq(0,10,length.out=100)))
ggplot(df,aes(x=x,y=y,fill=fill)) + geom_tile() +
  scale_fill_gradient(low="white",high="red", breaks = seq(0,10,by=2),name="Contamination", #Change this for the numbering under the legend and title
                      labels = scales::label_number(accuracy=1)) +
  theme_minimal() +
  guides(fill = guide_colorbar(
    direction = "horizontal",
    barwidth = 24,
    barheight = 6,
    ticks = TRUE,
    ticks.colour = "black"
  )) +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        legend.position="top",
        legend.title=element_text(size=48, face="bold"),
        legend.title.position = "top",
        legend.text = element_text(size=48),
        axis.text=element_blank(),
        legend.key.height=unit(5,"cm"),
        legend.ticks.length = unit(0.30,"cm"),
        panel.border=element_rect(color="black",fill=NA))
       
