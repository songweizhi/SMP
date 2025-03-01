import os
import pandas as pd


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


def among_host_species_variability(otu_table_txt, metadata_txt, tax_rank, color_code_txt, interested_sample_txt, interested_group_txt, op_prefix, op_dir, default_color):

    # define file name
    otu_table_subset    = '%s/%s_otu_table_subset.txt' % (op_dir, op_prefix)
    grouping_txt        = '%s/%s_grouping.txt'         % (op_dir, op_prefix)
    group_color_txt     = '%s/%s_color.txt'            % (op_dir, op_prefix)
    output_plot         = '%s/%s_plot.pdf'             % (op_dir, op_prefix)
    output_boxplot      = '%s/%s_boxplot.pdf'          % (op_dir, op_prefix)

    # get path to rarefaction_R
    pwd_current_file  = os.path.realpath(__file__)
    current_file_path = '/'.join(pwd_current_file.split('/')[:-1])
    among_host_species_variability_R     = '%s/among_host_species_variability.R' % current_file_path
    if os.path.isfile(among_host_species_variability_R) is False:
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

    sample_group_txt_handle = open(grouping_txt, 'w')
    sample_group_txt_handle.write('SampleID\tSampleGroup\n')
    group_set = set()
    for sample in sorted(list(sample_with_grp_info)):
        sample_grp = sample_group_dict[sample]
        sample_group_txt_handle.write('%s\t%s\n' % (sample, sample_grp))
        group_set.add(sample_grp)
    sample_group_txt_handle.close()

    color_code_dict = dict()
    for each_sample_type in open(color_code_txt):
        sample_type_split = each_sample_type.strip().split('\t')
        sample_type = sample_type_split[0]
        sample_color = sample_type_split[1]
        color_code_dict[sample_type] = sample_color

    group_color_txt_handle = open(group_color_txt, 'w')
    group_color_txt_handle.write('GroupID\tGroupColor\n')
    for grp in sorted(list(group_set)):
        grp_color = color_code_dict.get(grp, default_color)
        group_color_txt_handle.write('%s\t%s\n' % (grp, grp_color))
    group_color_txt_handle.close()

    among_host_species_variability_cmd = 'Rscript %s -i %s -g %s -c %s -b %s -p %s' % (among_host_species_variability_R, otu_table_to_plot, grouping_txt, group_color_txt, output_boxplot, output_plot)
    print(among_host_species_variability_cmd)
    os.system(among_host_species_variability_cmd)

    print('Done')


########################################################################################################################

# file in
otu_table_txt           = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s07_AllSamples_unoise_otu_table_NonEU.txt'
sample_metadata_txt     = '/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250228.txt'
host_taxon_rank         = 'f'
color_code_sample_txt   = '/Users/songweizhi/Desktop/SMP/00_metadata/color_code_sample_type.txt'
interested_sample_txt   = None
interested_group_txt    = '/Users/songweizhi/Desktop/SMP/00_metadata/sample_Sponge_Water_Sediment.txt'
default_color           = '#E67E22'

# file out
op_dir                  = '/Users/songweizhi/Desktop/SMP/Analysis_2_among_host_species_variability'
op_prefix               = 'Sponge_Water_Sediment'

########################################################################################################################

among_host_species_variability(otu_table_txt, sample_metadata_txt, host_taxon_rank, color_code_sample_txt, interested_sample_txt, interested_group_txt, op_prefix, op_dir, default_color)

