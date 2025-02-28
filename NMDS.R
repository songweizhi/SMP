# https://environmentalcomputing.net/graphics/multivariate-vis/mds/
# A rule of thumb is that stress values should ideally be less than 0.2 or even 0.1.

suppressMessages(suppressWarnings(library(vegan)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(optparse)))

################################################################################

option_list = list(
  make_option(c("-i", "--in"),   type="character", default=NULL, help="otu table"),
  make_option(c("-o", "--out"),   type="character", default=NULL, help="output plot"));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
input_table_txt    = opt$in
output_plot        = opt$out

# input_table_txt = '/Users/songweizhi/Desktop/nmds/otu_table_subset_T_with_group.txt'
# output_plot     = '/Users/songweizhi/Desktop/nmds/otu_table_subset_T_with_group.pdf'

#input_table_txt = '/Users/songweizhi/Desktop/nmds/otu_table_subset_T_with_Dive.txt'
#output_plot     = '/Users/songweizhi/Desktop/nmds/otu_table_subset_T_with_Dive.pdf'

################################################################################

df_in <- read.csv(file = input_table_txt, sep = '\t', header = TRUE)

otu_count <- df_in[3:ncol(df_in)]

otu_count.mds <- metaMDS(comm = otu_count, distance = "bray", trace = FALSE, autotransform = FALSE)

MDS_xy <- data.frame(otu_count.mds$points)
MDS_xy$Group <- df_in$Group

# get stress_value
stress_value = round(otu_count.mds$stress, digits=3)

p = ggplot(MDS_xy, aes(MDS1, MDS2, color = Group)) +
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

ggsave(filename=output_plot, plot = p, device = cairo_pdf, width =250, height =200, units = "mm")
