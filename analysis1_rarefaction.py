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


def rarefaction(otu_table_txt, metadata_txt, color_code_txt, interested_sample_txt, interested_group_txt, default_color, op_prefix, op_dir):

    # get path to rarefaction_R
    pwd_current_file  = os.path.realpath(__file__)
    current_file_path = '/'.join(pwd_current_file.split('/')[:-1])
    rarefaction_R     = '%s/rarefaction.R' % current_file_path
    if os.path.isfile(rarefaction_R) is False:
        print('rarefaction.R not found, program exited!')
        exit()

    # define file name
    otu_table_subset    = '%s/%s_otu_table_subset.txt' % (op_dir, op_prefix)
    grouping_txt        = '%s/%s_grouping.txt'         % (op_dir, op_prefix)
    group_color_txt     = '%s/%s_color.txt'            % (op_dir, op_prefix)
    output_plot         = '%s/%s_rarefaction.pdf'      % (op_dir, op_prefix)

    otu_table_sample_list = open(otu_table_txt).readline().strip().split('\t')[1:]

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
    if (interested_sample_txt is None) and (interested_group_txt is None):
        interested_sample_set = otu_table_sample_list
    elif (interested_sample_txt is not None) and (interested_group_txt is None):
        for each_sample in open(interested_sample_txt):
            interested_sample_set.add(each_sample.strip())
    elif (interested_sample_txt is None) and (interested_group_txt is not None):

        interested_group_set = set()
        for each_grp in open(interested_group_txt):
            interested_group_set.add(each_grp.strip())

        for each_sample in sample_source_dict:
            sample_source = sample_source_dict[each_sample]
            if sample_source in interested_group_set:
                interested_sample_set.add(each_sample)

    shared_sample_set, uniq_to_otu_table, uniq_to_interested = get_shared_uniq_elements(otu_table_sample_list, interested_sample_set)

    if len(uniq_to_otu_table) > 0:
        print('Samples uniq to %s:' % otu_table_txt.split('/')[-1])
        print(','.join(sorted(list(uniq_to_otu_table))))
    if len(uniq_to_interested) > 0:
        print('Samples uniq to %s:' % metadata_txt.split('/')[-1])
        print(','.join(sorted(list(uniq_to_interested))))

    # subset OTU table
    otu_table_to_plot = otu_table_txt
    if len(shared_sample_set) < len(otu_table_sample_list):
        subset_df(otu_table_txt, otu_table_subset, shared_sample_set)
        otu_table_to_plot = otu_table_subset

    sample_group_txt_handle = open(grouping_txt, 'w')
    sample_group_txt_handle.write('SampleID\tSampleGroup\n')
    group_set = set()
    for sample in sorted(shared_sample_set):
        sample_grp = sample_source_dict.get(sample, 'na')
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

    rarefaction_cmd = 'Rscript %s -i %s -g %s -c %s -o %s' % (rarefaction_R, otu_table_to_plot, grouping_txt, group_color_txt, output_plot)
    print(rarefaction_cmd)
    os.system(rarefaction_cmd)

    print('Done')


########################################################################################################################

# file in
otu_table_txt           = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s07_AllSamples_unoise_otu_table_noEU_n10_pct10_min20000.txt'
sample_metadata_txt     = '/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250403.txt'
color_code_sample_txt   = '/Users/songweizhi/Desktop/SMP/00_metadata/color_code_sample_type.txt'
default_color           = '#999999'
# file out
op_dir                  = '/Users/songweizhi/Desktop/SMP/Analysis_1_Rarefaction'


interested_group_txt    = None
interested_sample_txt   = '/Users/songweizhi/Desktop/4_resequence.txt'
op_prefix               = '4_resequence'


color_code_sample_txt   = '/Users/songweizhi/Desktop/SMP/00_metadata/color_code_sample_type.txt'
interested_sample_txt   = None
interested_group_txt    = '/Users/songweizhi/Desktop/SMP/source_Sponge_Coral_Water_Sediment.txt'
op_prefix               = 'Sponge_Coral_Water_Sediment'


color_code_sample_txt   = '/Users/songweizhi/Desktop/SMP/00_metadata/color_code_sample_type.txt'
interested_sample_txt   = '/Users/songweizhi/Desktop/SMP/coral_sample_with_bcarcoding_35_with_Water_Sediment_74.txt'
interested_group_txt    = None
op_prefix               = 'Coral_Water_Sediment_74'

########################################################################################################################

rarefaction(otu_table_txt, sample_metadata_txt, color_code_sample_txt, interested_sample_txt, interested_group_txt, default_color, op_prefix, op_dir)
