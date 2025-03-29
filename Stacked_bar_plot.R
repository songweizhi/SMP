suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(optparse)))

################################################################################

option_list = list(
  make_option(c("-i", "--table_in"), type="character", default=NULL,    help="input table"),
  make_option(c("-o", "--pdf_out"),  type="character", default=NULL,    help="output plot"),
  make_option(c("-x", "--width"),    type="double",    default=NULL,    help="the width of sankey plot"));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input_table = opt$table_in
output_plot = opt$pdf_out
plot_width = opt$width

################################################################################

df = read.csv(input_table, sep='\t', header = T)

sample_num = length(unique(df$Sample))
plot_height = sample_num*0.15
plot_height <- ifelse(plot_height < 12, 8, plot_height)
pdf(output_plot, width =plot_width,  height =plot_height)
ggplot(df, aes(x = Sample, y = Abundance, fill = Taxa)) + 
  geom_bar(stat = "identity", position = 'stack') +
  # scale_y_discrete(expand = c(0.01,0.01))+
  facet_grid(Group ~ ., scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5,vjust = 0.5),
        axis.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5),
        axis.ticks.y = element_blank(), # hide tick
        strip.text.x = element_text(angle = 90),
        strip.text.y = element_text(angle = 0, hjust = 0), panel.background = element_blank()) +
  coord_flip() # plot vertically
dev.off()

