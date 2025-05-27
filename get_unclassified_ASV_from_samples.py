import pandas as pd


def df_take_percentage_by_col(df_in_file, df_out_file):
    df_in = pd.read_csv(df_in_file, sep='\t', header=0, index_col=0)
    df_out = df_in.div(df_in.sum()) * 100
    df_out = df_out.round(3)
    df_out.to_csv(df_out_file, sep='\t')


def subset_df(file_in, file_out, rows_to_keep_set, cols_to_keep_set):

    column_name_pos     = 0
    row_name_pos        = 0
    sep_symbol          = '\t'
    exclude_rows_cols   = False

    ###################################### get the id of rows and cols to subset #######################################

    # put all row and col headers in list
    df = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_row_header_list = df.index.values.tolist()
    df_col_header_list = df.columns.values.tolist()

    # read in rows_to_keep_file
    rows_found_set = set()
    rows_missing_set = set()
    for each_r in rows_to_keep_set:
        if each_r in df_row_header_list:
            rows_found_set.add(each_r)
        else:
            rows_missing_set.add(each_r)

    # read in cols_to_keep_file
    cols_found_set = set()
    cols_missing_set = set()
    for each_c in cols_to_keep_set:
        if each_c in df_col_header_list:
            cols_found_set.add(each_c)
        else:
            cols_missing_set.add(each_c)

    # report
    if len(rows_missing_set) > 0:
        print('The following rows are missing from the dataframe:\n%s'    % ','.join(sorted(list(rows_missing_set))))

    if len(cols_missing_set) > 0:
        print('The following columns are missing from the dataframe:\n%s' % ','.join(sorted(list(cols_missing_set))))

    ####################################################################################################################

    rows_to_keep_set = set()
    cols_to_keep_set = set()
    if exclude_rows_cols is False:
        rows_to_keep_set = rows_found_set
        cols_to_keep_set = cols_found_set
    else:
        for each_row in df_row_header_list:
            if each_row not in rows_found_set:
                rows_to_keep_set.add(each_row)
        for each_col in df_col_header_list:
            if each_col not in cols_found_set:
                cols_to_keep_set.add(each_col)

    # turn sets into lists
    rows_to_keep_list_sorted = sorted(list(rows_to_keep_set))
    cols_to_keep_list_sorted = sorted(list(cols_to_keep_set))

    if len(rows_to_keep_list_sorted) == 0:
        if len(cols_to_keep_list_sorted) == 0:
            subset_df = df.loc[:, :]
        else:
            subset_df = df.loc[:, cols_to_keep_list_sorted]
    else:
        if len(cols_to_keep_list_sorted) == 0:
            subset_df = df.loc[rows_to_keep_list_sorted, :]
        else:
            subset_df = df.loc[rows_to_keep_list_sorted, cols_to_keep_list_sorted]

    subset_df_filtered = subset_df[subset_df.sum(axis=1) != 0]
    subset_df_filtered.to_csv(file_out, sep=sep_symbol)


########################################################################################################################

# file in
otu_tax_txt     = '/Users/songweizhi/Desktop/SMP/02_Usearch16S_20250526_356/s08_AllSamples_unoise_nc.BLCA.combined.updated.txt'
otu_table_txt   = '/Users/songweizhi/Desktop/SMP/02_Usearch16S_20250526_356/s07_AllSamples_unoise_otu_table_noEU_min10000.txt'
sample_txt      = 'JL316_B16_3,JL306_B14_4,JL306_B14_5,JL311_B03,JL313_B28_3'

# file out
op_otu_table_txt   = '/Users/songweizhi/Desktop/SMP/02_Usearch16S_20250526_356/s07_AllSamples_unoise_otu_table_noEU_min10000_demo.txt'

########################################################################################################################

sample_list = sample_txt.split(',')
sample_list_sorted = sorted(sample_list)

unclassified_otu_set = set()
for each_line in open(otu_tax_txt):
    each_line_split = each_line.strip().split('\t')
    otu_id  = each_line_split[0]
    otu_tax = each_line_split[1]
    if otu_tax == 'Unclassified':
        unclassified_otu_set.add(otu_id)
unclassified_otu_list_sorted = sorted(list(unclassified_otu_set))

subset_df(otu_table_txt, op_otu_table_txt, unclassified_otu_list_sorted, sample_list_sorted)
