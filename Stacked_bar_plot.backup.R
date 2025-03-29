suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(optparse)))

################################################################################

option_list = list(
  make_option(c("-i", "--table_in"), type="character", default=NULL,    help="input table"),
  make_option(c("-o", "--pdf_out"),  type="character", default=NULL,    help="output plot"));
  make_option(c("-x", "--width"),    type="double",    default=NULL,    help="the width of sankey plot"),
  make_option(c("-y", "--height"),   type="double",    default=NULL,    help="the height of sankey plot"),
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input_table = opt$table_in
output_plot = opt$pdf_out
plot_width = opt$width
plot_height = opt$height

################################################################################

df = read.csv(input_table, sep='\t', header = T)

pdf(output_plot, width =plot_width, height =plot_height)

ggplot(df, aes(x = Sample, y = Abundance, fill = Taxa)) + 
  geom_bar(stat = "identity", position = 'stack') +
  facet_grid(~ Group, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
  # panel.grid = element_blank()
  )

dev.off()
