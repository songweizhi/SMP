import os
import pandas as pd


otu_table_txt       = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s07_AllSamples_unoise_otu_table_NonEU.txt'
otu_table_stats_txt = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s07_AllSamples_unoise_otu_table_NonEU_count.txt'


otu_df                = pd.read_csv(otu_table_txt, sep='\t', header=0, index_col=0)
otu_count_dict        = otu_df.sum().to_dict()

otu_table_stats_txt_handle = open(otu_table_stats_txt, 'w')
for sample in {k: v for k, v in sorted(otu_count_dict.items(), key=lambda item: item[1])}:
    value = otu_count_dict[sample]
    otu_table_stats_txt_handle.write('%s\t%s\n' % (sample, value))
otu_table_stats_txt_handle.close()
