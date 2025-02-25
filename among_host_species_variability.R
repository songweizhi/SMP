suppressMessages(suppressWarnings(library(vegan)))
suppressMessages(suppressWarnings(library(purrr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(optparse)))

################################################################################

option_list = list(
  make_option(c("-i", "--otu"),     type="character", default=NULL, help="otu table"),
  make_option(c("-g", "--group"),   type="character", default=NULL, help="sample group"),
  make_option(c("-c", "--color"),   type="character", default=NULL, help="group color"),
  make_option(c("-b", "--boxplot"), type="character", default=NULL, help="output plot"),
  make_option(c("-p", "--plot"),    type="character", default=NULL, help="output plot"));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

otu_table_txt    = opt$otu
grouping_txt     = opt$group
color_code_txt   = opt$color
output_plot      = opt$plot
output_boxplot   = opt$boxplot

# otu_table_txt  = '/Users/songweizhi/Desktop/SMP/among_host_species_variability/otu_table.txt'
# grouping_txt   = '/Users/songweizhi/Desktop/SMP/among_host_species_variability/grouping.txt'
# color_code_txt = '/Users/songweizhi/Desktop/SMP/01_metadata/color_code_sample_type.txt'
# output_plot    = '/Users/songweizhi/Desktop/1.pdf'
# output_boxplot = '/Users/songweizhi/Desktop/2.pdf'

# otu_table_txt  = '/Users/songweizhi/Desktop/among_host_species_variability/otu_table_subset.txt'
# grouping_txt   = '/Users/songweizhi/Desktop/among_host_species_variability/grouping.txt'
# color_code_txt = '/Users/songweizhi/Desktop/among_host_species_variability/color.txt'
# output_plot    = '/Users/songweizhi/Desktop/among_host_species_variability/plot.pdf'
# output_boxplot = '/Users/songweizhi/Desktop/among_host_species_variability/boxplot.pdf'

################################################################################

otu_table = read.csv(otu_table_txt, sep='\t', row.names = 1, header = T)
otu_table_T = t(otu_table)
grouping_df = read.delim(grouping_txt)
color_df    = read.delim(color_code_txt)
sample_list = as.character(rownames(otu_table_T))

# Use map() to get corresponding metadata from the SampleGroup column
sample_group_list <- map(sample_list, ~ grouping_df$SampleGroup[grouping_df$SampleID == .x])
sample_group_list = as.character(sample_group_list)

# get color list for group
sample_group_list_unique = unique(sample_group_list)
sample_group_list_unique_sorted = sort(sample_group_list_unique)
sample_group_list_unique_sorted_color <- map(sample_group_list_unique_sorted, ~ color_df$GroupColor[color_df$GroupID == .x])
sample_group_list_unique_sorted_color = as.character(sample_group_list_unique_sorted_color)

# Bray-Curtis distances between samples
dis = vegdist(otu_table_T)

## Calculate multivariate dispersions
mod = betadisper(dis, sample_group_list)

# Extract average distance to median per group
average_dist_to_median <- tapply(mod$distances, sample_group_list, mean)

sorted_groups <- names(sort(average_dist_to_median))

sorted_groups_color <- map(sorted_groups, ~ color_df$GroupColor[color_df$GroupID == .x])
sorted_groups_color = as.character(sorted_groups_color)

##################################### plot #####################################

# Plot the groups and distances to centroids on the first two PCoA axes
pdf(output_plot, width =9, height =9)
par(cex.axis=1, mar=c(5, 5, 5, 5))
plot(mod)
dev.off()

pdf(output_boxplot, width =12, height =9)
par(cex.axis=1, mar=c(5, 10, 5, 5))
boxplot(mod, 
        #names = sorted_groups,
        col = sample_group_list_unique_sorted_color,
        xlab = 'Distance to group centroids', 
        ylab = '', # Sponge species (environments) 
        las=1,     # control direction of group label
        horizontal=TRUE)

dev.off()

################################################################################

#rm(list=ls())
