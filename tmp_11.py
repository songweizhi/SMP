import os
import math
import math
import os.path
import random
import argparse
import pandas as pd
import seaborn as sns


def transpose_csv(file_in, file_out, sep_symbol, column_name_pos, row_name_pos):

    csv = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_csv = pd.DataFrame(data=csv)
    transposed_csv = df_csv.T
    transposed_csv.to_csv(file_out, sep=sep_symbol)


########################################################################################################################

# file in
otu_table_txt       = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s07_AllSamples_unoise_otu_table_noEU_n10_pct10_min20000.txt'
sample_txt          = '/Users/songweizhi/Desktop/10.txt'
classification_txt  = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s08_AllSamples_unoise_nc.blca.GTDB.2.txt'
otu_table_txt_t     = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s07_AllSamples_unoise_otu_table_noEU_n10_pct10_min20000_T.txt'

########################################################################################################################

otu_tax_dict = dict()
for otu_tax in open(classification_txt):
    otu_tax_split = otu_tax.strip().split('\t')
    otu_id = otu_tax_split[0]
    tax_str = otu_tax_split[1]
    otu_tax_dict[otu_id] = tax_str

transpose_csv(otu_table_txt, otu_table_txt_t, '\t', 0, 0)

interested_sample_set = set()
for sample in open(sample_txt):
    interested_sample_set.add(sample.strip())
print(interested_sample_set)


print('Sample\tzotu1\tzotu27\tzotu28')
col_header_list = []
line_index = 0
for line in open(otu_table_txt_t):
    line_split = line.strip().split('\t')
    if line_index == 0:
        col_header_list = line_split
    else:
        sample_id = line_split[0]
        otu_count_list = line_split[1:]
        if sample_id in interested_sample_set:
            otu_count_dict = dict()
            for (otu_id, otu_count) in zip(col_header_list, otu_count_list):
                if int(otu_count) != 0:
                    otu_count_dict[otu_id] = int(otu_count)

            total_ar_count = 0
            for each_otu in sorted(list(otu_count_dict.keys())):
                otu_tax = otu_tax_dict[each_otu]
                otu_count = otu_count_dict[each_otu]
                if 'd__Archaea' in otu_tax:
                    total_ar_count += otu_count
                    #print(each_otu, otu_count, otu_tax)

            zotu1_count = otu_count_dict.get('Zotu1', 0)
            zotu27_count = otu_count_dict.get('Zotu27', 0)
            zotu28_count = otu_count_dict.get('Zotu28', 0)

            zotu1_pct = float("{0:.2f}".format(zotu1_count*100/total_ar_count))
            zotu27_pct = float("{0:.2f}".format(zotu27_count*100/total_ar_count))
            zotu28_pct = float("{0:.2f}".format(zotu28_count*100/total_ar_count))

            print('%s\tZotu1\t%s/%s\t%s'  % (sample_id, zotu1_count,  total_ar_count, zotu1_pct))
            print('%s\tZotu27\t%s/%s\t%s' % (sample_id, zotu27_count, total_ar_count, zotu27_pct))
            print('%s\tZotu28\t%s/%s\t%s' % (sample_id, zotu28_count, total_ar_count, zotu28_pct))
            print('%s\t%s\t%s\t%s' % (sample_id, zotu1_pct, zotu27_pct, zotu28_pct))

            print()


    line_index += 1




