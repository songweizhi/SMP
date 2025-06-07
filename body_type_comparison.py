import os
import pandas as pd


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


def body_type_comparison(metadata_txt, interested_group_txt, otu_table_txt, otu_classification_txt, otu_tax_rank, op_dir, op_prefix):

    # get path to rarefaction_R
    pwd_current_file  = os.path.realpath(__file__)
    current_file_path = '/'.join(pwd_current_file.split('/')[:-1])
    Stacked_bar_plot_R     = '%s/Stacked_bar_plot.R' % current_file_path
    if os.path.isfile(Stacked_bar_plot_R) is False:
        print('Stacked_bar_plot.R not found, program exited!')
        exit()

    ############################################## define output file name #############################################

    otu_table_subset            = '%s/%s_otu_table_subset.txt'      % (op_dir, op_prefix)
    tax_table_txt               = '%s/%s_taxa_table.txt'            % (op_dir, op_prefix)
    tax_table_txt_for_ggplot    = '%s/%s_taxa_table_ggplot.txt'     % (op_dir, op_prefix)
    output_plot                 = '%s/%s_body_type_comparison.pdf'  % (op_dir, op_prefix)

    ################################################ subsample OTU table ###############################################

    otu_table_sample_list = open(otu_table_txt).readline().strip().split('\t')[1:]

    interested_source_set = set()
    for each_grp in open(interested_group_txt):
        interested_source_set.add(each_grp.strip())

    sample_source_dict = dict()
    fauna_to_sample_dict = dict()
    sample_to_fauna_dict = dict()
    sample_to_body_type_dict = dict()
    col_index = dict()
    line_num_index = 0
    for each_line in open(metadata_txt):
        line_num_index += 1
        line_split = each_line.strip().split('\t')
        if line_num_index == 1:
            col_index = {key: i for i, key in enumerate(line_split)}
        else:
            sample_id = line_split[col_index['Sample_id']]
            sample_source = line_split[col_index['Source']]
            sample_host_taxon = line_split[col_index['Host_taxonomy']]
            sample_fauna = line_split[col_index['Dive_fauna']]
            sample_body_type = line_split[col_index['Body_type']]
            sample_to_fauna_dict[sample_id] = sample_fauna
            if sample_source in interested_source_set:
                if sample_fauna not in fauna_to_sample_dict:
                    fauna_to_sample_dict[sample_fauna] = set()
                fauna_to_sample_dict[sample_fauna].add(sample_id)
                sample_to_body_type_dict[sample_id] = sample_body_type
                sample_source_dict[sample_id] = sample_source

    sample_to_group_dict = dict()
    sample_to_label_dict = dict()
    for each_fauna in fauna_to_sample_dict:
        fauna_sample_set = fauna_to_sample_dict[each_fauna]
        if len(fauna_sample_set) > 1:
            for sample in fauna_sample_set:
                sample_source = sample_source_dict[sample]
                sample_body_type = sample_to_body_type_dict[sample]
                sample_label = '%s (%s)' % (sample, sample_body_type)
                sample_to_label_dict[sample] = sample_label
                sample_to_group_dict[sample] = '%s (%s)' % (each_fauna, sample_source)

    # get shared and uniq samples
    shared_sample_set, uniq_to_otu_table, uniq_to_interested = get_shared_uniq_elements(otu_table_sample_list, sample_to_group_dict.keys())

    # check singleton group
    to_plot_fauna_dict = dict()
    for each_sample in shared_sample_set:
        sample_fauna =   sample_to_fauna_dict[each_sample]
        if sample_fauna not in to_plot_fauna_dict:
            to_plot_fauna_dict[sample_fauna] = set()
        to_plot_fauna_dict[sample_fauna].add(each_sample)

    singleton_sample_set = set()
    for each_fauna_to_plot in to_plot_fauna_dict:
        current_sample_set = to_plot_fauna_dict[each_fauna_to_plot]
        if len(current_sample_set) <= 1:
            for each in current_sample_set:
                singleton_sample_set.add(each)

    shared_sample_set_without_singleton = set()
    for each_sample in shared_sample_set:
        if each_sample not in singleton_sample_set:
            shared_sample_set_without_singleton.add(each_sample)

    if len(uniq_to_otu_table) > 0:
        print('Samples uniq to %s (%s):' % (otu_table_txt.split('/')[-1], len(uniq_to_otu_table)))
        print(','.join(sorted(list(uniq_to_otu_table))))
        print()
    if len(uniq_to_interested) > 0:
        print('Samples uniq to %s (%s):' % (metadata_txt.split('/')[-1], len(uniq_to_interested)))
        print(','.join(sorted(list(uniq_to_interested))))
        print()
    if len(singleton_sample_set) > 0:
        print('Singletons (%s):' % len(singleton_sample_set))
        print(','.join(sorted(list(singleton_sample_set))))
        print()

    # subset OTU table
    otu_table_to_plot = otu_table_txt
    if len(shared_sample_set_without_singleton) < len(otu_table_sample_list):
        subset_df(otu_table_txt, otu_table_subset, shared_sample_set_without_singleton)
        otu_table_to_plot = otu_table_subset

    ####################################################################################################################

    otu_tax_dict = dict()
    for otu_tax in open(otu_classification_txt):
        otu_tax_split = otu_tax.strip().split('\t')
        otu_id = otu_tax_split[0]
        tax_str = otu_tax_split[1]
        tax_str_split = tax_str.split(';')
        needed_tax = '%s__' % otu_tax_rank
        for each_rank in tax_str_split:
            if each_rank.startswith(otu_tax_rank):
                needed_tax = each_rank
        otu_tax_dict[otu_id] = needed_tax

    otu_table_df = pd.read_csv(otu_table_to_plot, sep='\t', header=0, index_col=0)
    otu_table_df['Taxon'] = otu_table_df.index.map(otu_tax_dict)
    tax_table_df = otu_table_df.groupby('Taxon', as_index=False).sum()
    tax_table_df_normalized = tax_table_df.copy()
    tax_table_df_normalized.iloc[:, 1:] = tax_table_df.iloc[:, 1:].div(tax_table_df.iloc[:, 1:].sum())
    tax_table_df_normalized_rounded = tax_table_df_normalized.round(3)
    tax_table_df_normalized_rounded.to_csv(tax_table_txt, sep='\t', header=True, index=False)
    tax_table_df_normalized_rounded.set_index('Taxon', inplace=True)

    #################### get stacked bar plot ####################

    tax_table_txt_for_ggplot_handle = open(tax_table_txt_for_ggplot, 'w')
    tax_table_txt_for_ggplot_handle.write('Group\tSample\tTaxa\tAbundance\n')
    header_list = []
    line_index = 0
    for each_line in open(tax_table_txt):
        line_split = each_line.strip().split('\t')
        if line_index == 0:
            header_list = line_split
        else:
            tax_id = line_split[0]
            for (sample_id, taxa_abund) in zip(header_list[1:], line_split[1:]):
                sample_group = sample_to_group_dict[sample_id]
                sample_label = sample_to_label_dict[sample_id]
                tax_table_txt_for_ggplot_handle.write('%s\t%s\t%s\t%s\n' % (sample_group, sample_label, tax_id, taxa_abund))
        line_index += 1
    tax_table_txt_for_ggplot_handle.close()

    # run R script
    nmds_cmd = 'Rscript %s -i %s -o %s' % (Stacked_bar_plot_R, tax_table_txt_for_ggplot, output_plot)
    print(nmds_cmd)
    os.system(nmds_cmd)

    print('Done')

########################################################################################################################

sample_metadata_txt         = '/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250228.txt'
interested_sample_group_txt = '/Users/songweizhi/Desktop/SMP/Analysis_4_Community_composition/sample_All_17.txt'
interested_sample_group_txt = '/Users/songweizhi/Desktop/SMP/Analysis_4_Community_composition/samples_Sponge.txt'
op_dir                      = '/Users/songweizhi/Desktop/SMP/Analysis_5_Body_type'
op_prefix                   = 'Sponge'
otu_table_txt               = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s07_AllSamples_unoise_otu_table_nonEU.txt'
classification_txt          = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s08_AllSamples_unoise_nc.blca.gtdb.2.txt'
otu_tax_rank                = 'd'

########################################################################################################################

body_type_comparison(sample_metadata_txt, interested_sample_group_txt, otu_table_txt, classification_txt, otu_tax_rank, op_dir, op_prefix)
