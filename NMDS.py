import os
import pandas as pd


def transpose_csv(file_in, file_out, sep_symbol, column_name_pos, row_name_pos):

    csv = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_csv = pd.DataFrame(data=csv)
    transposed_csv = df_csv.T
    transposed_csv.to_csv(file_out, sep=sep_symbol)


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


def nmds(otu_table_txt, metadata_txt, tax_rank, interested_group_txt, op_dir):

    # define file name
    otu_table_subset        = '%s/otu_table_subset.txt'                 % op_dir
    otu_table_subset_t      = '%s/otu_table_subset_T.txt'               % op_dir
    otu_table_subset_t_meta = '%s/otu_table_subset_T_with_group.txt'    % op_dir
    output_plot             = '%s/nmds.pdf'                             % op_dir

    # get path to rarefaction_R
    pwd_current_file  = os.path.realpath(__file__)
    current_file_path = '/'.join(pwd_current_file.split('/')[:-1])
    nmds_R     = '%s/NMDS.R' % current_file_path
    if os.path.isfile(nmds_R) is False:
        print('among_host_species_variability.R not found, program exited!')
        exit()

    interested_source_set = set()
    for each_grp in open(interested_group_txt):
        interested_source_set.add(each_grp.strip())

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
            sample_id            = line_split[col_index['Sample_id']]
            sample_source        = line_split[col_index['Source']]
            sample_host_tax_str = line_split[col_index['Host_taxonomy']]
            sample_host_tax_split = sample_host_tax_str.split(';')
            if sample_source in interested_source_set:
                if sample_source == 'Water':
                    sample_group_dict[sample_id] = 'Water'
                elif sample_source == 'Sediment':
                    sample_group_dict[sample_id] = 'Sediment'
                else:
                    needed_tax = '%s__' % tax_rank
                    for each_rank in sample_host_tax_split:
                        if each_rank.startswith(tax_rank):
                            needed_tax = each_rank
                    sample_group_dict[sample_id] = needed_tax

    # get shared and uniq samples
    otu_table_sample_list = open(otu_table_txt).readline().strip().split('\t')[1:]
    shared_sample_set, uniq_to_otu_table, uniq_to_interested = get_shared_uniq_elements(otu_table_sample_list, sample_group_dict.keys())

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
        if sample_group == ('%s__' % tax_rank):
            sample_without_tax.add(sample)
        else:
            sample_with_grp_info.add(sample)

    if len(sample_without_tax) > 0:
        print('Samples with unknown classification at %s level:' % tax_rank)
        print(','.join(sorted(list(sample_without_tax))))
        print()

    # subset OTU table
    otu_table_to_plot = otu_table_txt
    if len(sample_with_grp_info) < len(otu_table_sample_list):
        subset_df(otu_table_txt, otu_table_subset, sample_with_grp_info)
        otu_table_to_plot = otu_table_subset

    # transpose otu table
    transpose_csv(otu_table_subset, otu_table_subset_t, '\t', 0, 0)

    #add group info to transposed otu table
    otu_table_subset_t_meta_handle = open(otu_table_subset_t_meta, 'w')
    line_index = 0
    for each_line in open(otu_table_subset_t):
        line_split = each_line.strip().split('\t')
        if line_index == 0:
            otu_table_subset_t_meta_handle.write('Sample\tGroup' + each_line)
        else:
            sample_grp = sample_group_dict[line_split[0]]
            line_split.insert(1, sample_grp)
            otu_table_subset_t_meta_handle.write('\t'.join(line_split) + '\n')

        line_index += 1

    otu_table_subset_t_meta_handle.close()

    # run R script
    nmds_cmd = 'Rscript %s -i %s -o %s' % (nmds_R, otu_table_to_plot, output_plot)
    print(nmds_cmd)
    os.system(nmds_cmd)

    print('Done')


########################################################################################################################

# file in
otu_table_txt           = '/Users/songweizhi/Desktop/NMDS/s07_AllSamples_unoise_otu_table.txt'
sample_metadata_txt     = '/Users/songweizhi/Desktop/NMDS/metadata.txt'
interested_group_txt    = '/Users/songweizhi/Desktop/NMDS/samples_Sponge_Water_Sediment.txt'
taxon_rank              = 'f'

# file out
op_dir                  = '/Users/songweizhi/Desktop/NMDS'

########################################################################################################################

nmds(otu_table_txt, sample_metadata_txt, taxon_rank, interested_group_txt, op_dir)
