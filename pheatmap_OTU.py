import os
import pandas
import argparse
import numpy as np
import pandas as pd


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def reorder_df_col(file_in, file_out, sep_symbol, column_name_pos, row_name_pos, column_order_list):
    csv = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_csv = pd.DataFrame(data=csv)
    df_csv = df_csv.reindex(column_order_list, axis=1)
    df_csv.to_csv(file_out, sep=sep_symbol)


def pheatmap_OTU(args):

    otu_table_txt           = args['otu']
    metadata_txt            = args['m']
    classification_txt      = args['otu_c']
    otu_label_rank          = args['olr']
    otu_annotation_rank     = args['oar']
    sample_label_rank       = args['slr']
    sample_annotation_rank  = args['sar']
    zero_as_na              = args['na']
    output_plot             = args['o']

    sample_label_rank_list = sample_label_rank.split(',')
    otu_label_rank_list    = otu_label_rank.split(',')

    # get path to rarefaction_R
    pwd_current_file = os.path.realpath(__file__)
    current_file_path = '/'.join(pwd_current_file.split('/')[:-1])
    pheatmap_OTU_R = '%s/pheatmap_OTU.R' % current_file_path
    if os.path.isfile(pheatmap_OTU_R) is False:
        print('pheatmap_OTU.R not found, program exited!')
        exit()

    plot_name, plot_path, plot_base, plot_ext = sep_path_basename_ext(output_plot)

    if os.path.isdir(plot_path) is False:
        os.mkdir(plot_path)

    otu_table_txt_desc_tmp  = '%s/%s_OTU_table_with_desc.tmp.txt'   % (plot_path, plot_base)
    otu_table_txt_desc      = '%s/%s_OTU_table_with_desc.txt'       % (plot_path, plot_base)
    col_annotation_txt      = '%s/%s_col_annotation.txt'            % (plot_path, plot_base)
    row_annotation_txt      = '%s/%s_row_annotation.txt'            % (plot_path, plot_base)

    otu_tax_dict = dict()
    for otu_tax in open(classification_txt):
        otu_tax_split = otu_tax.strip().split('\t')
        otu_id = otu_tax_split[0]
        tax_str = otu_tax_split[1]
        otu_tax_dict[otu_id] = tax_str

    sample_group_dict = dict()
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
            sample_host_taxon = line_split[col_index['Host_Taxonomy_NCBI']]
            sample_group = sample_source
            if sample_host_taxon != 'na':
                sample_group = sample_host_taxon
            sample_group_dict[sample_id] = sample_group

    # write out
    longest_sample_label = 0
    longest_otu_label = 0
    annotation_to_sample_dict = dict()
    label_to_sample_dict = dict()
    sample_num = 0
    otu_num = -1
    otu_table_txt_desc_tmp_handle = open(otu_table_txt_desc_tmp, 'w')
    col_annotation_txt_handle = open(col_annotation_txt, 'w')
    col_annotation_txt_handle.write('\tSample_%s\n' % sample_annotation_rank)
    row_annotation_txt_handle = open(row_annotation_txt, 'w')
    row_annotation_txt_handle.write('\tOTU_%s\n' % otu_annotation_rank)
    line_index = 0
    for each_line in open(otu_table_txt):
        each_line_split = each_line.strip().split('\t')
        sample_num = len(each_line_split)
        if line_index == 0:
            label_list_with_desc = []
            for each_sample in each_line_split:
                sample_group_str = sample_group_dict[each_sample]
                if sample_group_str == 'Water':
                    sample_group_str = 'water'
                if sample_group_str == 'Sediment':
                    sample_group_str = 'sediment'

                sample_label_tax = sample_group_str
                sample_annotation_tax = sample_group_str
                if ';' in sample_group_str:
                    sample_group_str_split = sample_group_str.split(';')
                    sample_label_tax_list = []
                    for each_rank in sample_group_str_split:
                        rank_abbrev = each_rank.split('__')[0]
                        if rank_abbrev in sample_label_rank_list:
                            sample_label_tax_list.append(each_rank)
                        if each_rank.startswith('%s__' % sample_annotation_rank):
                            sample_annotation_tax = each_rank
                    sample_label_tax = ';'.join(sample_label_tax_list)
                sample_new_label = '%s__%s' % (sample_label_tax, each_sample)
                if len(sample_new_label) > longest_sample_label:
                    longest_sample_label = len(sample_new_label)
                col_annotation_txt_handle.write('%s\t%s\n' % (sample_new_label, sample_annotation_tax))

                # get annotation_to_sample_dict
                if sample_annotation_tax not in annotation_to_sample_dict:
                    annotation_to_sample_dict[sample_annotation_tax] = set()
                annotation_to_sample_dict[sample_annotation_tax].add(sample_new_label)

                # label_to_sample_dict
                if sample_label_tax not in label_to_sample_dict:
                    label_to_sample_dict[sample_label_tax] = set()
                label_to_sample_dict[sample_label_tax].add(sample_new_label)

                label_list_with_desc.append(sample_new_label)
            otu_table_txt_desc_tmp_handle.write('\t%s\n' % '\t'.join(label_list_with_desc))
        else:
            otu_id = each_line_split[0]
            otu_count_list = each_line_split[1:]
            otu_tax_str = otu_tax_dict.get(otu_id, 'na')
            otu_tax_str_split = otu_tax_str.split(';')

            # turn 0 to NA
            otu_count_list_0_as_na = []
            for value in otu_count_list:
                if float(value) == 0:
                    otu_count_list_0_as_na.append('NA')
                else:
                    otu_count_list_0_as_na.append(value)
            if zero_as_na is True:
                otu_count_list = otu_count_list_0_as_na

            otu_label_tax = 'na'
            otu_annotation_tax = 'na'
            if otu_tax_str == 'Unclassified':
                otu_label_tax = 'unclassified'
                otu_annotation_tax = 'unclassified'
            else:
                otu_label_tax_list = []
                for each_rank in otu_tax_str_split:
                    rank_abbrev = each_rank.split('__')[0]
                    if rank_abbrev in otu_label_rank_list:
                        otu_label_tax_list.append(each_rank)

                    if each_rank.startswith('%s__' % otu_label_rank):
                        otu_label_tax = each_rank

                    if each_rank.startswith('%s__' % otu_annotation_rank):
                        otu_annotation_tax = each_rank

                otu_label_tax = ';'.join(otu_label_tax_list)

            otu_new_label = '%s__%s' % (otu_label_tax, otu_id)

            if len(otu_new_label) > longest_otu_label:
                longest_otu_label = len(otu_new_label)

            row_annotation_txt_handle.write('%s\t%s\n' % (otu_new_label, otu_annotation_tax))
            otu_table_txt_desc_tmp_handle.write('%s\t%s\n' % (otu_new_label, '\t'.join(otu_count_list)))
        line_index += 1
        otu_num += 1
    otu_table_txt_desc_tmp_handle.close()
    col_annotation_txt_handle.close()
    row_annotation_txt_handle.close()

    #################### reorder df columns and rows ####################

    df = pd.read_csv(otu_table_txt_desc_tmp, sep='\t', header=0, index_col=0)
    df = df[sorted(df.columns)]
    df = df.sort_index(ascending=True)
    df.to_csv(otu_table_txt_desc, sep='\t')

    #################### run R script ####################

    # run R script
    plot_width   = (sample_num*0.2) + (longest_otu_label*0.05)    + 5
    plot_height  = (otu_num*0.35)   + (longest_sample_label*0.03) + 2
    plot_width   = round(plot_width)
    plot_height  = round(plot_height)
    pheatmap_cmd = 'Rscript %s -i %s -o %s -r %s -c %s -x %s -y %s' % (pheatmap_OTU_R, otu_table_txt_desc, output_plot, row_annotation_txt, col_annotation_txt, plot_width, plot_height)
    print(pheatmap_cmd)
    os.system(pheatmap_cmd)

    print('Done!')


if __name__ == '__main__':

    blast_parser = argparse.ArgumentParser()
    blast_parser.add_argument('-m',     required=True,  help='metadata file')
    blast_parser.add_argument('-otu',   required=True,  help='otu table')
    blast_parser.add_argument('-otu_c', required=True,  help='otu classification_txt')
    blast_parser.add_argument('-olr',   required=True,  help='otu_label_rank')
    blast_parser.add_argument('-oar',   required=True,  help='otu_annotation_rank')
    blast_parser.add_argument('-slr',   required=True,  help='sample_label_rank')
    blast_parser.add_argument('-sar',   required=True,  help='sample_annotation_rank')
    blast_parser.add_argument('-na',    required=False, action="store_true",  help='turn 0 to NA')
    blast_parser.add_argument('-o',     required=True,  help='output pdf')
    args = vars(blast_parser.parse_args())
    pheatmap_OTU(args)
