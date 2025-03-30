# https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/

suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(pheatmap)))
suppressMessages(suppressWarnings(library(dendextend)))

################################################################################

option_list = list(
  make_option(c("-i", "--table_in"), type="character", default=NULL,    help="input table"),
  make_option(c("-o", "--plot_out"), type="character", default=NULL,    help="output plot"),
  make_option(c("-r", "--rowa"),     type="character", default=NULL,    help="row annotation file"),
  make_option(c("-c", "--cola"),     type="character", default=NULL,    help="column annotation file"));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

otu_table_txt       = opt$table_in
output_plot         = opt$plot_out
row_annotation_txt  = opt$rowa
col_annotation_txt  = opt$cola

################################################################################

otu_df = read.csv(file = otu_table_txt, sep = '\t', header = TRUE, row.names=1)
col_annotation_df = read.csv(file = col_annotation_txt, sep = '\t', header = TRUE, row.names=1)
otu_annotation_df = read.csv(file = row_annotation_txt, sep = '\t', header = TRUE, row.names=1)

# plot and save
pdf(output_plot, width =12,  height =6)
my_heatmap <- pheatmap(otu_df,
                       cluster_rows=T, cutree_rows = 2,
                       cluster_cols=T,
                       annotation_row = otu_annotation_df,
                       annotation_col = col_annotation_df)
dev.off()
