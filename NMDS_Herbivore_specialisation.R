# https://environmentalcomputing.net/graphics/multivariate-vis/mds/
# A rule of thumb is that stress values should ideally be less than 0.2 or even 0.1.

library(vegan)
library(ggplot2)


Herbivores <- read.csv(file = "/Users/songweizhi/PycharmProjects/SMP/demo_data/NMDS/Herbivore_specialisation.csv", header = TRUE)

Herb_community <- Herbivores[5:11]

Herb_community.mds <- metaMDS(comm = Herb_community, distance = "bray", trace = FALSE, autotransform = FALSE)

MDS_xy <- data.frame(Herb_community.mds$points)
MDS_xy$Habitat <- Herbivores$Habitat
MDS_xy$DayNight <- Herbivores$DayNight

# get stress_value
stress_value = Herb_community.mds$stress
stress_value = round(stress_value, digits=3)

ggplot(MDS_xy, aes(MDS1, MDS2, color = Habitat)) +
  annotate("text", x = max(MDS_xy$MDS1), y = min(MDS_xy$MDS1) - 0.5, label = paste("Stress: ", stress_value), hjust = 1, vjust = 0, size = 5, color = "black") +  
  geom_point() +
  theme_bw() + 
  theme(
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks.x = element_line(linewidth =0.2, color = 'black'),
    axis.ticks.y = element_line(linewidth =0.2, color = 'black'),
    axis.title.y = element_text(size = 12, angle = 90,face = "plain", color = 'black'),
    axis.title.x = element_text(size = 12, face = "plain", color = 'black'),
    axis.text.x = element_text(size = 12, face = "plain", color = 'black'),
    axis.text.y = element_text(size = 12,hjust = 1, face = "plain", color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_text(size=10, face = "plain",color = 'black'),
    legend.text = element_text(size=10, face = "plain",color = 'black')
  )
