import os
import pandas as pd


def remove_df_small_values(table_in, value_cutoff, decimal_round_at, table_out):

    df_in           = pd.read_csv(table_in, sep='\t', header=0, index_col=0)
    df_pct          = df_in.div(df_in.sum(axis=0), axis=1)
    df_pct          = df_pct.round(decimal_round_at)
    df_pct_filtered = df_pct.where(df_pct >= value_cutoff, other=0)
    df_pct_no_zero  = df_pct_filtered[df_pct_filtered.sum(axis=1) != 0]
    df_pct_no_zero.to_csv(table_out, sep='\t')


########################################################################################################################

otu_table_txt           = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s07_AllSamples_unoise_otu_table_nonEU.txt'
abundance_cutoff        = 0.0001
decimal_round_at        = 4
otu_table_txt_filtered  = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s07_AllSamples_unoise_otu_table_nonEU_0.0001.txt'

########################################################################################################################

remove_df_small_values(otu_table_txt, abundance_cutoff, decimal_round_at, otu_table_txt_filtered)
