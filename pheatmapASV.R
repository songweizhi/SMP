# https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/

suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(pheatmap)))
suppressMessages(suppressWarnings(library(dendextend)))

################################################################################

option_list = list(
  make_option(c("-i", "--table_in"),  type="character", default=NULL,    help="input table"),
  make_option(c("-o", "--plot_out"),  type="character", default=NULL,    help="output plot"),
  make_option(c("-s", "--split_col"), type="character", default=NULL,    help="positions to split column"),
  make_option(c("-r", "--rowa"),      type="character", default=NULL,    help="row annotation file"),
  make_option(c("-c", "--cola"),      type="character", default=NULL,    help="column annotation file"),
  make_option(c("-x", "--width"),     type="double",    default=NULL,    help="the width of plot"),
  make_option(c("-y", "--height"),    type="double",    default=NULL,   help="the height of plot"));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

otu_table_txt       = opt$table_in
output_plot         = opt$plot_out
row_annotation_txt  = opt$rowa
col_annotation_txt  = opt$cola
plot_width          = opt$width
plot_height         = opt$height
number_str          = opt$split_col

number_vector = as.integer(unlist(strsplit(number_str, ",")))

################################################################################

otu_df = read.csv(file = otu_table_txt, sep = '\t', header = TRUE, row.names=1)
col_annotation_df = read.csv(file = col_annotation_txt, sep = '\t', header = TRUE, row.names=1)
otu_annotation_df = read.csv(file = row_annotation_txt, sep = '\t', header = TRUE, row.names=1)

# plot and save
pdf(output_plot, width =plot_width,  height =plot_height)
my_heatmap <- pheatmap(otu_df,
                       na_col = '#000000',
                       cluster_rows=FALSE, # cutree_rows = 2,
                       cluster_cols=FALSE, # cutree_cols = 10,
                       # gaps_row = '',
                       gaps_col = number_vector,
                       #annotation_row = otu_annotation_df,
                       #annotation_col = col_annotation_df
                       )
invisible(dev.off())
rm(list=ls())
