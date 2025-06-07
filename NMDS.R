# https://environmentalcomputing.net/graphics/multivariate-vis/mds/
# A rule of thumb is that stress values should ideally be less than 0.2 or even 0.1.

suppressMessages(suppressWarnings(library(plyr)))
suppressMessages(suppressWarnings(library(vegan)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggalt)))
suppressMessages(suppressWarnings(library(ggforce)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(concaveman)))

################################################################################

option_list = list(
  make_option(c("-i", "--file_in"),   type="character", default=NULL, help="otu table"),
  make_option(c("-o", "--file_out"),  type="character", default=NULL, help="output plot"));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input_table_txt = opt$file_in
output_plot     = opt$file_out

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

chulls_for_polygon <- ddply(MDS_xy, .(Group), function(MDS_xy) MDS_xy[chull(MDS_xy$MDS1, MDS_xy$MDS2), ])

# NMDS plot with customized colors, shapes, and outlines
p = ggplot(MDS_xy, aes(MDS1, MDS2, color=Group)) +
  labs(x = "NMDS1", y = "NMDS2") +
  annotate("text", x = max(MDS_xy$MDS1), y = min(MDS_xy$MDS2), label = paste("Stress: ", stress_value), hjust = 1, vjust = 0, size = 5, color = "black") +
  geom_point() +
  scale_colour_manual(values=setNames(df_in$Color, df_in$Group)) +
  
  # Or:
  # scale_color_manual(values = c("black_coral" = "#8E44AD", "octocoral" = "#28B463", "water" = "#3498DB", "sediment" = "#a77260"))+
  # geom_point(aes(shape = Group, color = Group), size =1.5) +
  # scale_shape_manual(values = c(15, 16, 17, 15, 16, 17, 15, 16, 17, 3, 8)) +
  # scale_color_manual(values = c("#959595", "#9f399f", "#3787c0", "#ffb939", "#959595", "#9f399f", "#3787c0", "#ffb939", "#959595", "black", "black")) +
  #geom_point(aes(shape = c(15, 16, 17, 15, 16, 17, 15, 16, 17, 3, 8), color = c("#959595", "#9f399f", "#3787c0", "#ffb939", "#959595", "#9f399f", "#3787c0", "#ffb939", "#959595", "black", "black")), size = 1.5) +
  #geom_point(aes(shape = c(15, 16, 17, 15, 16, 17, 15, 16, 17, 3, 8), color = Group), size = 1.5) +
  #scale_shape_manual(values = c(15, 16, 17, 15, 16, 17, 15, 16, 17, 3, 8)) +
  #scale_color_manual(values = c("#959595", "#9f399f", "#3787c0", "#ffb939", "#959595", "#9f399f", "#3787c0", "#ffb939", "#959595", "black", "black")) +
  #scale_color_manual(values = color_uniq) +

  ######################## draw polygon and/or ellipse #########################

  # geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=Group)) +           # good
  # geom_mark_polygon(expand=0, linewidth=0, aes(fill=Group))+                  # very good !!!
  # geom_polygon(aes(fill = value, group = id, subgroup = subid)) +
  # geom_polygon(aes(fill=Group, group=Group, alpha=0.1, rule="winding")) +
  # geom_encircle(aes(fill=Group)）+
  # geom_polygon(data=chulls, linewidth=0.5, linetype="dashed", fill=Group, alpha=0.3, aes(fill = Group, group=Group)) +
  
  geom_polygon(data=chulls_for_polygon, linewidth=0.5, alpha=0.35, aes(fill = Group, group=Group)) +
  scale_fill_manual(values = setNames(df_in$Color, df_in$Group)) +
  
  #geom_mark_ellipse(expand=0, linetype=1, linewidth=0.3, aes(fill=Group))+     # very good
  # stat_ellipse(level=0.95, linetype=2, lwd=0.3) +
  # stat_ellipse(geom = "polygon", aes(fill = Group), alpha = 0.2) +

  ##############################################################################

  theme_bw() +
  theme(
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks.x = element_line(linewidth = 0.2, color = 'black'),
    axis.ticks.y = element_line(linewidth = 0.2, color = 'black'),
    axis.title.x = element_text(size = 12, face = "plain", color = 'black'),
    axis.title.y = element_text(size = 12, angle = 90, face = "plain", color = 'black'),
    axis.text.x  = element_text(size = 12, face = "plain", color = 'black'),
    axis.text.y  = element_text(size = 12, hjust = 1, face = "plain", color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_text(size=10, face = "plain", color = 'black'),
    legend.text = element_text(size=10, face = "plain", color = 'black'))

ggsave(filename = output_plot, plot = p, device = cairo_pdf, width = 180, height = 120, units = "mm")
