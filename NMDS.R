# https://environmentalcomputing.net/graphics/multivariate-vis/mds/
# A rule of thumb is that stress values should ideally be less than 0.2 or even 0.1.

suppressMessages(suppressWarnings(library(vegan)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(optparse)))
library(dplyr)
library(ggplot2)
library(ggalt)
library(ggforce)
library(concaveman)

################################################################################

option_list = list(
  make_option(c("-i", "--file_in"),   type="character", default=NULL, help="otu table"),
  make_option(c("-o", "--file_out"),  type="character", default=NULL, help="output plot"));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
input_table_txt    = opt$file_in
output_plot        = opt$file_out

# input_table_txt    = '/Users/songweizhi/Desktop/SMP/Analysis_3_NMDS/Sponge_Water_Sediment_otu_table_subset_T_with_group.txt'
# output_plot        = '/Users/songweizhi/Desktop/SMP/Analysis_3_NMDS/Sponge_Water_Sediment_nmds.pdf'

################################################################################

df_in <- read.csv(file = input_table_txt, sep = '\t', header = TRUE)
otu_count <- df_in[5:(ncol(df_in))]
otu_count.mds <- metaMDS(comm = otu_count, distance = "bray", trace = FALSE, autotransform = FALSE)

MDS_xy <- data.frame(otu_count.mds$points)
MDS_xy$Group <- df_in$Group
MDS_xy$Shape <- df_in$Shape
MDS_xy$Color <- df_in$Color

# get stress_value
stress_value = round(otu_count.mds$stress, digits=3)
color_uniq = unique(df_in$Color)
shape_uniq = unique(df_in$Shape)
# shape_group = rep(1:3, length.out = length(df_in$Group))

grp_uniq = unique(df_in$Group)

p = ggplot(MDS_xy, aes(MDS1, MDS2, color=Group)) + #
  annotate("text", x = max(MDS_xy$MDS1), y = min(MDS_xy$MDS2), label = paste("Stress: ", stress_value), hjust = 1, vjust = 0, size = 5, color = "black") +
  geom_point() +
  # geom_point(aes(shape = Group, color = Group), size =1.5) +
  # scale_shape_manual(values = c(15, 16, 17, 15, 16, 17, 15, 16, 17, 3, 8)) +
  # scale_color_manual(values = c("#959595", "#9f399f", "#3787c0", "#ffb939", "#959595", "#9f399f", "#3787c0", "#ffb939", "#959595", "black", "black")) +

  #geom_point(aes(shape = c(15, 16, 17, 15, 16, 17, 15, 16, 17, 3, 8), color = c("#959595", "#9f399f", "#3787c0", "#ffb939", "#959595", "#9f399f", "#3787c0", "#ffb939", "#959595", "black", "black")), size = 1.5) +
  #geom_point(aes(shape = c(15, 16, 17, 15, 16, 17, 15, 16, 17, 3, 8), color = Group), size = 1.5) +
  #scale_shape_manual(values = c(15, 16, 17, 15, 16, 17, 15, 16, 17, 3, 8)) +
  #scale_color_manual(values = c("#959595", "#9f399f", "#3787c0", "#ffb939", "#959595", "#9f399f", "#3787c0", "#ffb939", "#959595", "black", "black")) +
  #scale_color_manual(values = color_uniq) +
  #print(MDS_xy$Group)
  # geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=Group)) +   # good
  geom_mark_ellipse(expand=0, linewidth=0, aes(fill=Group))+            # very good
  #geom_mark_ellipse(expand=0, linetype=1, linewidth=0.3, aes(fill=Group))+            # very good
  # stat_ellipse(level=0.95, linetype=2, lwd=0.3) +
  # stat_ellipse(geom = "polygon", aes(fill = Group), alpha = 0.2) +
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

ggsave(filename=output_plot, plot = p, device = cairo_pdf, width =180, height =120, units = "mm")
  
