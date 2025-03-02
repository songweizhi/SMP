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


def community_composition(metadata_txt, host_tax_rank, interested_group_txt, interested_sample_txt, otu_table_txt, classification_txt, otu_tax_rank, op_dir, op_prefix):

    # get path to rarefaction_R
    pwd_current_file  = os.path.realpath(__file__)
    current_file_path = '/'.join(pwd_current_file.split('/')[:-1])
    Stacked_bar_plot_R     = '%s/Stacked_bar_plot.R' % current_file_path
    if os.path.isfile(Stacked_bar_plot_R) is False:
        print('Stacked_bar_plot.R not found, program exited!')
        exit()

    ############################################## define output file name #############################################

    otu_table_subset            = '%s/%s_otu_table_subset.txt'                  % (op_dir, op_prefix)
    tax_table_txt               = '%s/%s_taxa_table.txt'                        % (op_dir, op_prefix)
    tax_table_txt_for_ggplot    = '%s/%s_taxa_table_ggplot.txt'                 % (op_dir, op_prefix)
    output_plot_2               = '%s/%s_community_composition_stacked_bar.pdf' % (op_dir, op_prefix)

    ################################################ get interested_sample_set  ###############################################

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
            sample_id = line_split[col_index['Sample_id']]
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

    ################################################ subsample OTU table ###############################################

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
            sample_host_taxon    = line_split[col_index['Host_taxonomy']]
            if sample_id in interested_sample_set:

                if sample_host_taxon == 'na':
                    sample_group_dict[sample_id] = sample_source
                else:
                    sample_host_taxon_split = sample_host_taxon.split(';')
                    needed_tax = '%s__' % host_tax_rank
                    for each_rank in sample_host_taxon_split:
                        if each_rank.startswith(host_tax_rank):
                            needed_tax = each_rank

                    group_str = '%s__%s' % (sample_source, needed_tax)
                    sample_group_dict[sample_id] = group_str

    # get shared and uniq samples
    shared_sample_set, uniq_to_otu_table, uniq_to_interested = get_shared_uniq_elements(otu_table_sample_list, sample_group_dict.keys())

    if len(uniq_to_otu_table) > 0:
        print('Samples uniq to %s:' % otu_table_txt)
        print(','.join(sorted(list(uniq_to_otu_table))))
        print()
    if len(uniq_to_interested) > 0:
        print('Samples uniq to %s:' % metadata_txt)
        print(','.join(sorted(list(uniq_to_interested))))
        print()

    # subset OTU table
    otu_table_to_plot = otu_table_txt
    if len(shared_sample_set) < len(otu_table_sample_list):
        subset_df(otu_table_txt, otu_table_subset, shared_sample_set)
        otu_table_to_plot = otu_table_subset

    ####################################################################################################################

    otu_tax_dict = dict()
    for otu_tax in open(classification_txt):
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
    tax_table_df_normalized_rounded.set_index('Taxon', inplace=True)
    tax_table_df_normalized_rounded_no_zero  = tax_table_df_normalized_rounded[tax_table_df_normalized_rounded.sum(axis=1) != 0]
    tax_table_df_normalized_rounded_no_zero.to_csv(tax_table_txt, sep='\t', header=True, index=True)

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
                sample_group = sample_group_dict[sample_id]
                tax_table_txt_for_ggplot_handle.write('%s\t%s\t%s\t%s\n' % (sample_group, sample_id, tax_id, taxa_abund))
        line_index += 1
    tax_table_txt_for_ggplot_handle.close()

    # run R script
    nmds_cmd = 'Rscript %s -i %s -o %s' % (Stacked_bar_plot_R, tax_table_txt_for_ggplot, output_plot_2)
    print(nmds_cmd)
    os.system(nmds_cmd)

    print('Done')


########################################################################################################################

# file in
sample_metadata_txt         = '/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250228.txt'
group_host_at_rank          = 'g'  # None
otu_table_txt               = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s07_AllSamples_unoise_otu_table_nonEU.txt'
otu_table_txt               = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s07_AllSamples_unoise_otu_table_nonEU_0.0001.txt'

classification_txt          = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s08_AllSamples_unoise_nc.blca.gtdb.2.txt'
otu_tax_rank                = 'p'
interested_group_txt        = '/Users/songweizhi/Desktop/SMP/Analysis_4_Community_composition/sample_All_17.txt'
op_prefix                   = 'All_17_GTDB'

interested_group_txt        = '/Users/songweizhi/Desktop/SMP/00_metadata/sample_Coral_Water_Sediment.txt'
op_prefix                   = 'Coral_Water_Sediment_GTDB'

interested_group_txt        = '/Users/songweizhi/Desktop/SMP/00_metadata/samples_Coral.txt'
op_prefix                   = 'Coral_GTDB'

interested_sample_txt       = '/Users/songweizhi/Desktop/SMP/sample_Corals_with_abundant_archaea.txt'
interested_group_txt        = None
op_prefix                   = 'sample_coral_with_abundant_archaea_GTDB_0.0001'

op_dir                      = '/Users/songweizhi/Desktop/SMP/Analysis_4_Community_composition'

########################################################################################################################

community_composition(sample_metadata_txt, group_host_at_rank, interested_group_txt, interested_sample_txt, otu_table_txt, classification_txt, otu_tax_rank, op_dir, op_prefix)
