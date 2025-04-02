import os
import argparse
import pandas as pd


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


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


def community_composition(args):

    metadata_txt            = args['m']
    otu_table_txt           = args['otu']
    host_tax_rank           = args['hr']
    interested_source       = args['source']
    interested_sample_txt   = args['sample']
    classification_txt      = args['otu_c']
    otu_tax_rank            = args['mr']
    sample_group_txt        = args['g']
    sample_to_exclude_txt   = args['e']
    output_plot             = args['o']
    plot_width              = args['w']

    f_name, op_dir, f_base, f_ext = sep_path_basename_ext(output_plot)

    if os.path.isdir(op_dir) is False:
        os.mkdir(op_dir)

    # get path to rarefaction_R
    pwd_current_file   = os.path.realpath(__file__)
    current_file_path  = '/'.join(pwd_current_file.split('/')[:-1])
    Stacked_bar_plot_R = '%s/Stacked_bar_plot.R' % current_file_path
    if os.path.isfile(Stacked_bar_plot_R) is False:
        print('Stacked_bar_plot.R not found, program exited!')
        exit()

    ############################################## define output file name #############################################

    otu_table_subset            = '%s/%s_otu_table_subset.txt'                  % (op_dir, f_base)
    tax_table_txt               = '%s/%s_taxa_table.txt'                        % (op_dir, f_base)
    tax_table_txt_t             = '%s/%s_taxa_table_T.txt'                      % (op_dir, f_base)
    tax_table_txt_t_desc        = '%s/%s_taxa_table_T_desc.txt'                 % (op_dir, f_base)
    tax_table_txt_for_ggplot    = '%s/%s_taxa_table_ggplot.txt'                 % (op_dir, f_base)

    ################################################ get interested_sample_set  ###############################################

    sample_to_exclude_set = set()
    if os.path.isfile(sample_to_exclude_txt) is True:
        for sample in open(sample_to_exclude_txt):
            sample_to_exclude_set.add(sample.strip().split()[0])

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
    if (interested_sample_txt is None) and (interested_source is None):
        interested_sample_set = otu_table_sample_list
    elif (interested_sample_txt is not None) and (interested_source is None):
        for each_sample in open(interested_sample_txt):
            interested_sample_set.add(each_sample.strip())
    elif (interested_sample_txt is None) and (interested_source is not None):
        interested_source_set = interested_source.split(',')

        for each_sample in sample_source_dict:
            sample_source = sample_source_dict[each_sample]
            if sample_source in interested_source_set:
                interested_sample_set.add(each_sample)

    ################################################ subsample OTU table ###############################################

    sample_group_dict = dict()
    if sample_group_txt is not None:
        for each_sample in open(sample_group_txt):
            each_sample_split = each_sample.strip().split('\t')
            sample_group_dict[each_sample_split[0]] = each_sample_split[1]
    else:
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
                sample_id         = line_split[col_index['Sample_ID']]
                sample_source     = line_split[col_index['Source']]
                sample_host_taxon = line_split[col_index['Host_Taxonomy_NCBI']]
                if sample_id in interested_sample_set:
                    if sample_host_taxon == 'na':
                        sample_group_dict[sample_id] = sample_source
                    else:
                        sample_host_taxon_split = sample_host_taxon.split(';')
                        needed_tax_str = ''
                        for each_rank in sample_host_taxon_split:
                            rank_abbrev = each_rank.split('__')[0]
                            if rank_abbrev in host_tax_rank:
                                needed_tax_str += '%s;' % each_rank
                        needed_tax_str = needed_tax_str[:-1]
                        group_str = '%s__%s' % (sample_source, needed_tax_str)
                        sample_group_dict[sample_id] = group_str

    # get shared and uniq samples
    shared_sample_set, uniq_to_otu_table, uniq_to_interested = get_shared_uniq_elements(otu_table_sample_list, interested_sample_set)

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

    # transpose_csv
    transpose_csv(tax_table_txt, tax_table_txt_t, '\t', 0, 0)

    #################### get stacked bar plot ####################

    tax_table_txt_for_ggplot_handle = open(tax_table_txt_for_ggplot, 'w')
    tax_table_txt_for_ggplot_handle.write('Group\tSample\tTaxa\tAbundance\n')
    header_list = []
    line_index = 0
    sample_to_group_dict = dict()
    for each_line in open(tax_table_txt):
        line_split = each_line.strip().split('\t')
        if line_index == 0:
            header_list = line_split
        else:
            tax_id = line_split[0]
            for (sample_id, taxa_abund) in zip(header_list[1:], line_split[1:]):

                sample_group = sample_source_dict.get(sample_id, 'na')
                if sample_id in sample_group_dict:
                    sample_group = sample_group_dict[sample_id]
                if sample_id not in sample_to_exclude_set:
                    tax_table_txt_for_ggplot_handle.write('%s\t%s\t%s\t%s\n' % (sample_group, sample_id, tax_id, taxa_abund))
                    sample_to_group_dict[sample_id] = '%s__%s' % (sample_group, sample_id)
        line_index += 1
    tax_table_txt_for_ggplot_handle.close()

    # add desc to tax_table_txt_t
    tax_table_txt_t_desc_handle = open(tax_table_txt_t_desc, 'w')
    line_index = 0
    for each_line in open(tax_table_txt_t):
        if line_index == 0:
            tax_table_txt_t_desc_handle.write(each_line)
        else:
            each_line_split = each_line.strip().split('\t')
            sample_id = each_line_split[0]
            sample_value_list = each_line_split[1:]
            sample_id_with_desc = sample_to_group_dict.get(sample_id, sample_id)
            if sample_id not in sample_to_exclude_set:
                tax_table_txt_t_desc_handle.write('%s\t%s\n' % (sample_id_with_desc, '\t'.join(sample_value_list)))
        line_index += 1
    tax_table_txt_t_desc_handle.close()

    # run R script
    nmds_cmd = 'Rscript %s -i %s -o %s -x %s' % (Stacked_bar_plot_R, tax_table_txt_for_ggplot, output_plot, plot_width)
    print(nmds_cmd)
    os.system(nmds_cmd)

    os.remove(tax_table_txt_t)
    print('Done')


if __name__ == '__main__':

    blast_parser = argparse.ArgumentParser()
    blast_parser.add_argument('-m',      required=True,                          help='metadata file')
    blast_parser.add_argument('-o',      required=True,                          help='output pdf')
    blast_parser.add_argument('-otu',    required=True,                          help='otu table')
    blast_parser.add_argument('-otu_c',  required=True,                          help='otu classification_txt')
    blast_parser.add_argument('-mr',     required=True,                          help='microbiome taxon rank')
    blast_parser.add_argument('-sample', required=False, default=None,           help='interested samples')
    blast_parser.add_argument('-source', required=False, default=None,           help='interested sample source')
    blast_parser.add_argument('-g',      required=False, default=None,           help='sample group txt')
    blast_parser.add_argument('-hr',     required=False, default=None,           help='group sample by host taxonomy, specify taxon rank for grouping')
    blast_parser.add_argument('-e',      required=False, default='',             help='samples to exclude from the output')
    blast_parser.add_argument('-w',      required=False, default=18,type=int,    help='plot width')
    args = vars(blast_parser.parse_args())
    community_composition(args)
