# plot a rarefaction curve with ggplot
# https://stackoverflow.com/questions/47234809/coloring-rarefaction-curve-lines-by-metadata-vegan-package-phyloseq-package

# install phyloseq package from biocmanager
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install(c("phyloseq"))

suppressMessages(suppressWarnings(library(vegan)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(phyloseq)))
suppressMessages(suppressWarnings(library(tidyverse)))

################################################################################

option_list = list(
  make_option(c("-i", "--otu"),   type="character", default=NULL, help="otu table"),
  make_option(c("-g", "--group"), type="character", default=NULL, help="sample group"),
  make_option(c("-c", "--color"), type="character", default=NULL, help="group color"),
  make_option(c("-o", "--out"),   type="character", default=NULL, help="output plot"));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

otu_table_txt    = opt$otu
sample_group_txt = opt$group
color_code_txt   = opt$color
output_plot      = opt$out

#otu_table_txt    = '/Users/songweizhi/PycharmProjects/SpongeMicrobiomeProject/demo_data/rarefaction/otu_table.txt'
#sample_group_txt = '/Users/songweizhi/PycharmProjects/SpongeMicrobiomeProject/demo_data/rarefaction/sample_group.txt'
#color_code_txt   = '/Users/songweizhi/PycharmProjects/SpongeMicrobiomeProject/demo_data/rarefaction/color_code.txt'
#output_plot      = '/Users/songweizhi/PycharmProjects/SpongeMicrobiomeProject/demo_data/rarefaction/demo.pdf'

################################################################################

# data preparation, OTU as rows
otu_raw <- read.delim(otu_table_txt, row.names = 1)
otu_raw_t <-t(otu_raw)

# check the lowest reads among all samples
# rowSums(otu_raw_t)
raremax_raw <- min(rowSums(otu_raw_t))
# raremax_raw

message("running rarecurve()")
# fit into rarecurve matrix and have a general look to make sure everything is correct
out1 <- rarecurve(otu_raw_t, step=100, sample=raremax_raw, label=T)
message("finished rarecurve()")

# have a list each element corresponds to one sample
rare1 <- lapply(out1, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

# label the list
sample_names1 <- rownames(otu_raw_t)
names(rare1) <- sample_names1

# convert to data frame:
rare1 <- map_dfr(rare1, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

########################## add a grouping column to df #########################

color_df    = read.delim(color_code_txt)
grouping_df = read.delim(sample_group_txt)

#head(rare1)
# add a grouping column to df
rare1 <- map_dfr(grouping_df$SampleID, function(x){
  z <- rare1[rare1$sample == x,]                             
  SampleGroup <- grouping_df$SampleGroup[grouping_df$SampleID == x]
  z <- data.frame(z, SampleGroup)
  return(z)
  })
#head(rare1)

##################################### plot #####################################

message("running ggplot()")
p <- ggplot(data = rare1)+
  geom_line(aes(x = raw.read, y = OTU, group = sample, color = SampleGroup),linewidth =0.1)+ # color="black",, linetype=Compartment
  theme_bw() + theme(
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks.x = element_line(linewidth =0.2, color = 'black'),
    axis.ticks.y = element_line(linewidth =0.2, color = 'black'),
    axis.title.y = element_text(size = 7, angle = 90,face = "plain", color = 'black'),
    axis.title.x = element_text(size = 7, face = "plain", color = 'black'),
    axis.text.x = element_text(size = 6, face = "plain", color = 'black'),
    axis.text.y = element_text(size = 6,hjust = 1, face = "plain", color = 'black'),
    panel.spacing = unit(0.4, "lines"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    
    # legend style
    legend.title=element_text(size=7, face = "plain",color = 'black'),
    legend.text = element_text(size=7, face = "plain",color = 'black'),
    legend.margin=margin(1.5,3,1.5,3),
    legend.background=element_blank(),
    # legend.box.background = element_rect(color = 'black',linewidth  = 0.4),
    legend.box.background = element_blank(),
    # legend.box.margin=margin(0,0,-10,0),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.x = unit(0.1, 'cm'),
    legend.position  = 'right',
    legend.key = element_rect(colour = "transparent", fill = NA),
    legend.spacing.y = unit(0.1, 'cm'),
    aspect.ratio=1)+
  scale_colour_manual(breaks=color_df$GroupID, values = color_df$GroupColor) +
  
  labs( y = "Number of OTUs", x='Number of sequences', color = "")

ggsave(filename=output_plot, plot = p, device = cairo_pdf, width =100, height =70, units = "mm")

################################################################################

rm(list=ls())

# cd /Users/songweizhi/PycharmProjects/SpongeMicrobiomeProject/demo_data/rarefaction
# Rscript /Users/songweizhi/PycharmProjects/SpongeMicrobiomeProject/rarefaction.R -i otu_table.txt -g sample_group.txt -c color_code.txt -o rarefaction.pdf
