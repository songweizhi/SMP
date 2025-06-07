import os
import argparse
import pandas as pd

########################################################################################################################

# file in
otu_table_txt           = '/Users/songweizhi/Desktop/SMP/02_Usearch16S_20250526_356/s07_AllSamples_unoise_otu_table_noEU_min10000.txt'
interested_sample_txt   = '/Users/songweizhi/Desktop/SMP/00_metadata/Sponge_samples_20250524.txt'
sample_metadata_txt     = '/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250528.txt'
host_taxon_rank         = 'g'

# file out
otu_table_txt_subset    = '/Users/songweizhi/Desktop/SMP/Analysis_4_Co-occurrence_network/s07_AllSamples_unoise_otu_table_noEU_min10000_subset.txt'
sample_group_txt        = '/Users/songweizhi/Desktop/SMP/Analysis_4_Co-occurrence_network/sample_group.txt'

########################################################################################################################

# subset otu table and remove otus absent from all sample subset
interested_sample_set = set()
for sample in open(interested_sample_txt):
    interested_sample_set.add(sample.strip().split()[0])
otu_table_txt_subset_tmp = '%s.tmp.txt' % otu_table_txt_subset
subset_df_cmd = 'BioSAK subset_df -i %s -c %s -o %s' % (otu_table_txt, interested_sample_txt, otu_table_txt_subset_tmp)
os.system(subset_df_cmd)
df_tmp1 = pd.read_csv(otu_table_txt_subset_tmp, sep='\t', header=0, index_col=0)
df_tmp2 = df_tmp1[df_tmp1.sum(axis=1) != 0]
df_tmp2.to_csv(otu_table_txt_subset, sep='\t')
os.system('rm %s' % otu_table_txt_subset_tmp)



samples_in_otu_table = open(otu_table_txt_subset).readline().strip().split()

# read in metadata_txt
sample_group_dict = dict()
col_index = dict()
line_num_index = 0
for each_line in open(sample_metadata_txt):
    line_num_index += 1
    line_split = each_line.strip().split('\t')
    if line_num_index == 1:
        col_index = {key: i for i, key in enumerate(line_split)}
    else:
        sample_id = line_split[col_index['Sample_ID']]
        sample_source = line_split[col_index['Source']]
        sample_host_tax_str = line_split[col_index['Host_Taxonomy_NCBI']]

        if sample_id in samples_in_otu_table:
            if sample_source == 'Water':
                sample_group_dict[sample_id] = 'Water'
            elif sample_source == 'Sediment':
                sample_group_dict[sample_id] = 'Sediment'
            else:
                sample_host_tax_str_split = sample_host_tax_str.split(';')
                print(sample_host_tax_str)

                lowest_known_rank = ''
                for each_rank in sample_host_tax_str_split:
                    if each_rank.endswith('__') is False:
                        lowest_known_rank = each_rank
                sample_group_dict[sample_id] = sample_host_tax_str

sample_group_txt_handle = open(sample_group_txt, 'w')
for sample in sorted(sample_group_dict.keys()):
    sam_grp = sample_group_dict[sample]
    sample_group_txt_handle.write('%s\t%s\n' % (sample, sam_grp))
sample_group_txt_handle.close()
