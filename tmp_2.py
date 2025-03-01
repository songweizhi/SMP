import os
import pandas as pd

########################################################################################################################

otu_table_txt           = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s07_AllSamples_unoise_otu_table_nonEU.txt'
classification_txt      = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s08_AllSamples_unoise_nc.blca.gtdb.2.txt'
interested_sample_txt   = '/Users/songweizhi/Desktop/SMP/sample_Corals_with_abundant_archaea.txt'

########################################################################################################################

otu_tax_dict = dict()
for otu_tax in open(classification_txt):
    otu_tax_split = otu_tax.strip().split('\t')
    otu_id = otu_tax_split[0]
    tax_str = otu_tax_split[1]
    otu_tax_dict[otu_id] = tax_str

otu_df = pd.read_csv(otu_table_txt, sep='\t', header=0, index_col=0)

otu_count_dod = dict()
for each_sample in open(interested_sample_txt):
    sample_id = each_sample.strip()
    otu_count_dict = otu_df[sample_id].to_dict()
    otu_count_dict = {x:y for x,y in otu_count_dict.items() if y!=0}
    otu_count_dod[sample_id] = otu_count_dict

    for otu in {k: v for k, v in sorted(otu_count_dict.items(), key=lambda item: item[1])[::-1]}:
        otu_count = otu_count_dict[otu]
        otu_tax = otu_tax_dict[otu]
        if otu_count >= 5:
            print('%s\t%s\t%s\t%s' % (sample_id, otu, otu_count, otu_tax))
    print()
