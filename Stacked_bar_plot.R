suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(optparse)))

################################################################################

option_list = list(
  make_option(c("-i", "--table_in"),  type="character", default=NULL, help="input table"),
  make_option(c("-o", "--pdf_out"),   type="character", default=NULL, help="output plot"));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input_table = opt$table_in
output_plot = opt$pdf_out

################################################################################

df = read.csv(input_table, sep='\t', header = T)

sample_num = length(unique(df$Sample))
plot_height = sample_num*0.25
plot_height <- ifelse(plot_height < 12, 8, plot_height)
pdf(output_plot, width =12,  height =plot_height)
ggplot(df, aes(x = Sample, y = Abundance, fill = Taxa)) + 
  geom_bar(stat = "identity", position = 'stack') +
  facet_grid(Group ~ ., scales = "free", space = "free") +
  # facet_grid( ~ Group, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5,vjust = 0.5),
        strip.text.x = element_text(angle = 90),
        strip.text.y = element_text(angle = 0)) +
  # panel.grid = element_blank()
  coord_flip() # plot vertically
dev.off()
