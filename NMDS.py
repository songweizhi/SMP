import os
import math
import random
import os.path
import pandas as pd
import seaborn as sns


def transpose_csv(file_in, file_out, sep_symbol, column_name_pos, row_name_pos):

    csv = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_csv = pd.DataFrame(data=csv)
    transposed_csv = df_csv.T
    transposed_csv.to_csv(file_out, sep=sep_symbol)


def df_take_percentage_by_col(df_in_file, df_out_file):
    df_in = pd.read_csv(df_in_file, sep='\t', header=0, index_col=0)
    df_out = df_in.div(df_in.sum()) * 100
    df_out = df_out.round(3)
    df_out.to_csv(df_out_file, sep='\t')


def get_shared_uniq_elements(list_1, list_2):

    shared_set = set(list_1).intersection(list_2)

    list_1_uniq = []
    for e1 in list_1:
        if e1 not in shared_set:
            list_1_uniq.append(e1)
    list_2_uniq = []
    for e2 in list_2:
        if e2 not in shared_set:
            list_2_uniq.append(e2)

    return shared_set, list_1_uniq, list_2_uniq


def subset_df1(file_in, file_out, rows_to_rm_set, cols_to_rm_set, exclude_rows_cols):

    df_separator    = 'tab'
    column_name_pos = 0
    row_name_pos    = 0

    # setup separator
    if df_separator in ['tab', 'Tab', 'TAB']:
        sep_symbol = '\t'
    elif df_separator in ['comma', 'Comma', 'COMMA']:
        sep_symbol = ','
    else:
        print('Please specify separator as either tab or comma, program exited!')
        exit()

    ###################################### get the id of rows and cols to subset #######################################

    # put all row and col headers in list
    df = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_row_header_list = df.index.values.tolist()
    df_col_header_list = df.columns.values.tolist()

    # read in rows_to_keep_file
    rows_found_set = set()
    rows_missing_set = set()
    for each_r in rows_to_rm_set:
        if rows_to_rm_set is not None:
            if each_r in df_row_header_list:
                rows_found_set.add(each_r)
            else:
                rows_missing_set.add(each_r)

    # read in cols_to_keep_file
    cols_found_set = set()
    cols_missing_set = set()
    if cols_to_rm_set is not None:
        for each_c in cols_to_rm_set:
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

    subset_df.to_csv(file_out, sep=sep_symbol)


def rm_low_sum_cols(otu_table_in,abd_cutoff,otu_table_out):

    otu_table_df = pd.read_csv(otu_table_in, sep='\t', header=0, index_col=0)
    column_sums = otu_table_df.sum()
    column_sum_dict = column_sums.to_dict()
    otu_table_df_filtered = otu_table_df.loc[:, column_sums >= abd_cutoff]
    otu_table_df_filtered.to_csv(otu_table_out, sep='\t')

    # report
    # print('The following samples were removed from the OTU table:')
    # for sample in sorted(column_sum_dict.keys()):
    #     seq_count = column_sum_dict[sample]
    #     if seq_count < abd_cutoff:
    #         print('%s\t%s' % (sample, seq_count))


def subset_df(file_in, file_out, cols_to_keep_set):

    column_name_pos = 0
    row_name_pos    = 0
    sep_symbol      = '\t'

    ###################################### get the id of rows and cols to subset #######################################

    # put all row and col headers in list
    df = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)

    ####################################################################################################################

    # turn sets into lists
    rows_to_keep_list_sorted = set()
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

    subset_df.to_csv(file_out, sep=sep_symbol)


def get_color_list(color_num):

    if color_num <= 8:
        color_list_combined = ['#3787c0', '#39399f', '#ffb939', '#399f39', '#9f399f', '#fb694a', '#9f9f39', '#959595']

    elif 8 < color_num <= 16:
        color_list_combined = ['#2b7bba', '#89bedc', '#2e2e99', '#8a8acc', '#ffa500', '#ffc55c', '#2e992e', '#8acc8a', '#992e99', '#cc8acc', '#d52221', '#fc8161', '#99992e', '#cccc8a', '#5c5c5c', '#adadad']

    else:
        color_num_each = math.ceil(color_num/8) + 2
        color_list_1 = sns.color_palette('Blues', n_colors=color_num_each).as_hex()
        color_list_2 = sns.light_palette('navy',   n_colors=color_num_each).as_hex()
        color_list_3 = sns.light_palette('orange', n_colors=color_num_each).as_hex()
        color_list_4 = sns.light_palette('green',  n_colors=color_num_each).as_hex()
        color_list_5 = sns.light_palette('purple', n_colors=color_num_each).as_hex()
        color_list_6 = sns.color_palette('Reds',  n_colors=color_num_each).as_hex()
        color_list_7 = sns.light_palette('olive',  n_colors=color_num_each).as_hex()
        color_list_8 = sns.color_palette('Greys', n_colors=color_num_each).as_hex()

        color_list_combined = []
        for color_list in [color_list_1, color_list_2, color_list_3, color_list_4, color_list_5, color_list_6, color_list_7, color_list_8]:
            for color in color_list[2:][::-1]:
                color_list_combined.append(color)

    color_list_to_return = random.sample(color_list_combined, color_num)

    color_list_to_return_sorted = []
    for color_to_return in color_list_combined:
        if color_to_return in color_list_to_return:
            color_list_to_return_sorted.append(color_to_return)

    random.shuffle(color_list_to_return_sorted)

    return color_list_to_return_sorted


def nmds(otu_table_txt, otu_classification_txt, min_seq_num, metadata_txt, host_taxon_rank, interested_source, interested_sample_txt, sample_to_exclude_txt, op_dir, op_prefix):

    # define file name
    otu_table_subset                                    = '%s/%s_otu_table_subset1.txt'                             % (op_dir, op_prefix)
    otu_table_subset2_only_classified                   = '%s/%s_otu_table_subset2_only_classified.txt'             % (op_dir, op_prefix)
    otu_table_subset3_only_classified_abd_cutoff        = '%s/%s_otu_table_subset3_only_classified_abd.txt'         % (op_dir, op_prefix)
    otu_table_subset3_only_classified_abd_cutoff_pct    = '%s/%s_otu_table_subset3_only_classified_abd_pct.txt'     % (op_dir, op_prefix)
    otu_table_subset_t                                  = '%s/%s_otu_table_subset4_T.txt'                           % (op_dir, op_prefix)
    otu_table_subset_t_no_0_otu                         = '%s/%s_otu_table_subset5_T_no_0_otu.txt'                  % (op_dir, op_prefix)
    otu_table_subset_t_meta                             = '%s/%s_otu_table_subset6_T_with_group.txt'                % (op_dir, op_prefix)
    otu_table_subset_t_meta_genus_level                 = '%s/%s_otu_table_subset6_T_with_group_genus_level.txt'    % (op_dir, op_prefix)
    output_plot                                         = '%s/%s_nmds.pdf'                                          % (op_dir, op_prefix)
    output_plot_genus_level                             = '%s/%s_nmds_genus_level.pdf'                              % (op_dir, op_prefix)

    otu_table_sample_list = open(otu_table_txt).readline().strip().split('\t')[1:]
    otu_table_sample_list = open(otu_table_txt).readline().strip().split('\t')
    print(otu_table_sample_list)
    # get path to rarefaction_R
    pwd_current_file  = os.path.realpath(__file__)
    current_file_path = '/'.join(pwd_current_file.split('/')[:-1])
    nmds_R     = '%s/NMDS.R' % current_file_path
    if os.path.isfile(nmds_R) is False:
        print('among_host_species_variability.R not found, program exited!')
        exit()

    sample_to_exclude_set = set()
    if sample_to_exclude_txt is not None:
        if os.path.isfile(sample_to_exclude_txt) is True:
            for sample in open(sample_to_exclude_txt):
                sample_to_exclude_set.add(sample.strip().split()[0])

    if interested_source is not None:
        interested_source_set = interested_source.split(',')

    # read in sample_source_txt
    sample_source_dict = dict()
    col_index = dict()
    line_num_index = 0
    for each_line in open(metadata_txt):
        line_num_index += 1
        line_split = each_line.strip().split('\t')
        if line_num_index == 1:
            col_index = {key: i for i, key in enumerate(line_split)}
        else:
            sample_id = line_split[col_index['Sample_ID']]
            sample_source = line_split[col_index['Source']]
            sample_source_dict[sample_id] = sample_source

    # get interested_sample_set
    interested_sample_set = set()
    if (interested_sample_txt is None) and (interested_source is None):
        interested_sample_set = otu_table_sample_list
    elif (interested_sample_txt is not None) and (interested_source is None):
        for each_sample in open(interested_sample_txt):
            interested_sample_set.add(each_sample.strip().split()[0])
    elif (interested_sample_txt is None) and (interested_source is not None):
        interested_source_set = interested_source.split(',')

        for each_sample in sample_source_dict:
            sample_source = sample_source_dict[each_sample]
            if sample_source in interested_source_set:
                interested_sample_set.add(each_sample)

    # get shared and uniq samples
    shared_sample_set, uniq_to_otu_table, uniq_to_interested = get_shared_uniq_elements(otu_table_sample_list, interested_sample_set)
    print('uniq_to_interested')
    print(uniq_to_interested)

    # read in metadata_txt
    sample_group_dict = dict()
    col_index = dict()
    line_num_index = 0
    for each_line in open(metadata_txt):
        line_num_index += 1
        line_split = each_line.strip().split('\t')
        if line_num_index == 1:
            col_index = {key: i for i, key in enumerate(line_split)}
        else:
            sample_id            = line_split[col_index['Sample_ID']]
            sample_source        = line_split[col_index['Source']]
            sample_host_tax_str = line_split[col_index['Host_Taxonomy_NCBI']]
            sample_host_tax_split = sample_host_tax_str.split(';')

            if sample_id not in sample_to_exclude_set:
                if sample_id in shared_sample_set:
                    if sample_source == 'Water':
                        sample_group_dict[sample_id] = 'Water'
                    elif sample_source == 'Sediment':
                        sample_group_dict[sample_id] = 'Sediment'
                    else:
                        needed_tax = '%s__' % host_taxon_rank
                        for each_rank in sample_host_tax_split:
                            if each_rank.startswith(host_taxon_rank):
                                needed_tax = each_rank
                        sample_group_dict[sample_id] = needed_tax

    if len(uniq_to_otu_table) > 0:
        print('Samples uniq to %s:' % otu_table_txt)
        print(','.join(sorted(list(uniq_to_otu_table))))
        print()
    if len(uniq_to_interested) > 0:
        print('Samples uniq to %s:' % metadata_txt)
        print(','.join(sorted(list(uniq_to_interested))))
        print()

    sample_with_grp_info = set()
    sample_without_tax = set()
    for sample in shared_sample_set:
        sample_group = sample_group_dict[sample]
        if sample_group == ('%s__' % host_taxon_rank):
            sample_without_tax.add(sample)
        else:
            sample_with_grp_info.add(sample)

    if len(sample_without_tax) > 0:
        print('Samples with unknown classification at %s level:' % host_taxon_rank)
        print(','.join(sorted(list(sample_without_tax))))
        print()

    # subset OTU table by sample
    if len(sample_with_grp_info) < len(otu_table_sample_list):
        subset_df(otu_table_txt, otu_table_subset, sample_with_grp_info)

    # remove unclassified otu from OTU table
    unclassified_otu_set = set()
    otu_genus_dict = dict()
    for each in open(otu_classification_txt):
        each_split = each.strip().split('\t')
        otu_id = each_split[0]
        otu_tax = each_split[1]
        if otu_tax == 'Unclassified':
            unclassified_otu_set.add(otu_id)
        else:
            otu_tax_split = otu_tax.split(';')
            otu_genus = ''
            for each in otu_tax_split:
                if each.startswith('g__'):
                    otu_genus = each
            if otu_genus != 'g__':
                otu_genus_dict[otu_id] = otu_genus

    subset_df1(otu_table_subset, otu_table_subset2_only_classified, unclassified_otu_set, None, True)

    # rm_low_depth_sample
    rm_low_sum_cols(otu_table_subset2_only_classified, min_seq_num, otu_table_subset3_only_classified_abd_cutoff)

    # take percentage
    df_take_percentage_by_col(otu_table_subset3_only_classified_abd_cutoff, otu_table_subset3_only_classified_abd_cutoff_pct)
    table_to_use = otu_table_subset3_only_classified_abd_cutoff_pct

    # transpose otu table
    transpose_csv(table_to_use, otu_table_subset_t, '\t', 0, 0)
    rm_low_sum_cols(otu_table_subset_t, 1, otu_table_subset_t_no_0_otu)

    # get sample group list
    line_index = 0
    sample_grp_set = set()
    for each_line in open(otu_table_subset_t_no_0_otu):
        line_split = each_line.strip().split('\t')
        if line_index > 0:
            sample_grp = sample_group_dict[line_split[0]]
            sample_grp_set.add(sample_grp)
        line_index += 1

    # define color and shape
    color_num = math.ceil(len(set(sample_grp_set))/3)
    color_list = get_color_list(color_num)*999
    shape_list = ['15', '16', '17']*999

    group_list_test = []
    grp_to_color_dict = dict()
    grp_to_shape_dict = dict()
    grp_index = 0
    for grp in sorted(list(sample_grp_set)):
        if grp in ['Sediment', 'sediment']:
            grp_color = 'black'
            grp_shape = '3'
            group_list_test.append('sediment')
            grp_to_color_dict['sediment'] = grp_color
            grp_to_shape_dict['sediment'] = grp_shape
        elif grp in ['Water', 'water']:
            grp_color = 'black'
            grp_shape = '8'
            group_list_test.append('water')
            grp_to_color_dict['water'] = grp_color
            grp_to_shape_dict['water'] = grp_shape
        else:
            grp_color = color_list[grp_index]
            grp_shape = shape_list[grp_index]
            group_list_test.append(grp)
            grp_index += 1
        grp_to_color_dict[grp] = grp_color
        grp_to_shape_dict[grp] = grp_shape

    shape_list = [grp_to_shape_dict[i] for i in sorted(group_list_test)]
    color_list = [grp_to_color_dict[i] for i in sorted(group_list_test)]
    color_str_for_r = 'c("%s")' % '", "'.join(color_list)
    shape_str_for_r = 'c(%s)' % ', '.join(shape_list)

    #add group info to transposed otu table
    otu_table_subset_t_meta_handle = open(otu_table_subset_t_meta, 'w')
    line_index = 0
    for each_line in open(otu_table_subset_t_no_0_otu):
        line_split = each_line.strip().split('\t')
        if line_index == 0:
            otu_table_subset_t_meta_handle.write('Sample\tGroup\tColor\tShape' + each_line)
        else:
            sample_grp = sample_group_dict[line_split[0]]
            sample_color = grp_to_color_dict[sample_grp]
            sample_shape = grp_to_shape_dict[sample_grp]
            line_split.insert(1, sample_shape)
            line_split.insert(1, sample_color)
            line_split.insert(1, sample_grp)
            str_to_write = '\t'.join(line_split)
            str_to_write = str_to_write.replace('Water', 'water')
            str_to_write = str_to_write.replace('Sediment', 'sediment')
            otu_table_subset_t_meta_handle.write(str_to_write + '\n')
        line_index += 1
    otu_table_subset_t_meta_handle.close()

    # get otu table at genus level
    otu_table_otu_list = open(otu_table_subset_t_meta).readline().strip().split('\t')[4:]
    genus_set = set()
    for otu in otu_table_otu_list:
        otu_genus = otu_genus_dict.get(otu, '')
        if otu_genus != '':
            genus_set.add(otu_genus)

    genus_list_sorted = sorted(list(genus_set))

    otu_table_subset_t_meta_genus_level_handle = open(otu_table_subset_t_meta_genus_level, 'w')
    line_index = 0
    otu_id_list = []
    for each_sample in open(otu_table_subset_t_meta):
        each_sample_split = each_sample.strip().split('\t')
        if line_index == 0:
            otu_id_list = each_sample_split[4:]
            otu_table_subset_t_meta_genus_level_handle.write('%s\t%s\n' % ('\t'.join(each_sample_split[:4]), '\t'.join(genus_list_sorted)))
        else:
            genus_count_dict = dict()
            otu_count_list = each_sample_split[4:]
            for (otu_id, otu_count) in zip(otu_id_list, otu_count_list):
                otu_genus = otu_genus_dict.get(otu_id, '')
                if otu_genus != '':
                    if otu_genus not in genus_count_dict:
                        genus_count_dict[otu_genus] = 0
                    genus_count_dict[otu_genus] += float(otu_count)

            genus_count_list = []
            for genus in genus_list_sorted:
                genus_count = genus_count_dict.get(genus, 0)
                genus_count = float("{0:.4f}".format(genus_count))
                genus_count_list.append(str(genus_count))

            otu_table_subset_t_meta_genus_level_handle.write('%s\t%s\n' % ('\t'.join(each_sample_split[:4]) , '\t'.join(genus_count_list) ))
        line_index += 1
    otu_table_subset_t_meta_genus_level_handle.close()

    # run R script
    nmds_cmd = 'Rscript %s -i %s -o %s' % (nmds_R, otu_table_subset_t_meta, output_plot)
    print(nmds_cmd)
    os.system(nmds_cmd)

    nmds_cmd = 'Rscript %s -i %s -o %s' % (nmds_R, otu_table_subset_t_meta_genus_level, output_plot_genus_level)
    print(nmds_cmd)
    os.system(nmds_cmd)

    print('Done')

    print(pd.read_csv(otu_table_subset, sep='\t', header=0, index_col=0).shape)
    print(pd.read_csv(otu_table_subset2_only_classified, sep='\t', header=0, index_col=0).shape)
    print(pd.read_csv(otu_table_subset3_only_classified_abd_cutoff, sep='\t', header=0, index_col=0).shape)
    print(pd.read_csv(otu_table_subset_t, sep='\t', header=0, index_col=0).shape)
    print(pd.read_csv(otu_table_subset_t_no_0_otu, sep='\t', header=0, index_col=0).shape)


######################################################### Coral ########################################################

# # unclassified OTUs will be ignored from this analysis
# otu_table_txt           = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s07_AllSamples_unoise_otu_table_noEU_min20000_min_pct_0.1.txt'
# sample_metadata_txt     = '/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250511.txt'
# otu_classification_txt  = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s08_AllSamples_unoise_nc.blca.GTDB.2_updated_by_blast_vs_nt.txt'
# interested_sample_txt   = '/Users/songweizhi/Desktop/SMP/00_metadata/Coral_samples_with_barcoding_30_with_Water_Sediment_69.txt'
# minimum_seq_num         = 10000
# interested_source       = None
# sample_to_exclude_txt   = None
# op_dir                  = '/Users/songweizhi/Desktop/SMP/Analysis_3_NMDS'
#
# # group by coral family
# host_taxon_rank         = 'f'
# op_prefix               = 'coral_sample_with_barcoding_30_with_Water_Sediment_69_family_level'
#
# # group by coral genus
# host_taxon_rank         = 'g'
# op_prefix               = 'coral_sample_with_barcoding_30_with_Water_Sediment_69_genus_level'

######################################################## Sponge ########################################################

# unclassified OTUs will be ignored from this analysis
otu_table_txt           = '/Users/songweizhi/Desktop/SMP/02_Usearch16S_20250526_356/s07_AllSamples_unoise_otu_table_noEU_min10000_min0.1pct.txt'
sample_metadata_txt     = '/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250528.txt'
otu_classification_txt  = '/Users/songweizhi/Desktop/SMP/02_Usearch16S_20250526_356/s08_AllSamples_unoise_nc.BLCA.combined.updated.txt'
interested_sample_txt   = '/Users/songweizhi/Desktop/SMP/00_metadata/Sponge_samples_20250524.txt'
minimum_seq_num         = 10000
interested_source       = None
sample_to_exclude_txt   = None
op_dir                  = '/Users/songweizhi/Desktop/SMP/Analysis_3_NMDS'

# group by Sponge order
host_taxon_rank         = 'o'
op_prefix               = 'Sponge_Water_Sediment_113_by_sponge_order'

# group by Sponge family
host_taxon_rank         = 'f'
op_prefix               = 'Sponge_Water_Sediment_113_by_sponge_family'

# group by Sponge genus
host_taxon_rank         = 'g'
op_prefix               = 'Sponge_Water_Sediment_113_by_sponge_genus'

########################################################################################################################

nmds(otu_table_txt, otu_classification_txt, minimum_seq_num, sample_metadata_txt, host_taxon_rank, interested_source, interested_sample_txt, sample_to_exclude_txt, op_dir, op_prefix)
