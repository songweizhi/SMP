import os
import argparse


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def pheatmap_OTU(args):

    otu_table_txt           = args['otu']
    metadata_txt            = args['m']
    classification_txt      = args['otu_c']
    otu_label_rank          = args['olr']
    otu_annotation_rank     = args['oar']
    sample_label_rank       = args['slr']
    sample_annotation_rank  = args['sar']
    output_plot             = args['o']

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

    otu_table_txt_desc   = '%s/%s_OTU_table_with_desc.txt'  % (plot_path, plot_base)
    col_annotation_txt   = '%s/%s_col_annotation.txt'       % (plot_path, plot_base)
    row_annotation_txt   = '%s/%s_row_annotation.txt'       % (plot_path, plot_base)

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
    otu_table_txt_desc_handle = open(otu_table_txt_desc, 'w')
    col_annotation_txt_handle = open(col_annotation_txt, 'w')
    col_annotation_txt_handle.write('\tSample_%s\n' % sample_annotation_rank)
    row_annotation_txt_handle = open(row_annotation_txt, 'w')
    row_annotation_txt_handle.write('\tOTU_%s\n' % otu_annotation_rank)
    line_index = 0
    for each_line in open(otu_table_txt):
        each_line_split = each_line.strip().split('\t')
        if line_index == 0:
            label_list_with_desc = []
            for each_sample in each_line_split:
                sample_group_str = sample_group_dict[each_sample]
                sample_group_str_split = sample_group_str.split(';')
                sample_label_tax = 'na'
                sample_annotation_tax = 'na'
                for each_rank in sample_group_str_split:
                    if each_rank.startswith('%s__' % sample_label_rank):
                        sample_label_tax = each_rank
                    if each_rank.startswith('%s__' % sample_annotation_rank):
                        sample_annotation_tax = each_rank
                sample_new_label = '%s__%s' % (each_sample, sample_label_tax)
                col_annotation_txt_handle.write('%s\t%s\n' % (sample_new_label, sample_annotation_tax))
                label_list_with_desc.append(sample_new_label)
            otu_table_txt_desc_handle.write('\t%s\n' % '\t'.join(label_list_with_desc))
        else:
            otu_id = each_line_split[0]
            otu_tax_str = otu_tax_dict.get(otu_id, 'na')
            otu_tax_str_split = otu_tax_str.split(';')

            otu_label_tax = 'na'
            otu_annotation_tax = 'na'
            if otu_tax_str == 'Unclassified':
                otu_label_tax = 'Unclassified'
                otu_annotation_tax = 'Unclassified'
            else:
                for each_rank in otu_tax_str_split:
                    if each_rank.startswith('%s__' % otu_label_rank):
                        otu_label_tax = each_rank
                    if each_rank.startswith('%s__' % otu_annotation_rank):
                        otu_annotation_tax = each_rank

            otu_new_label = '%s__%s' % (otu_id, otu_label_tax)
            row_annotation_txt_handle.write('%s\t%s\n' % (otu_new_label, otu_annotation_tax))
            otu_table_txt_desc_handle.write('%s\t%s\n' % (otu_new_label, '\t'.join(each_line_split[1:])))
        line_index += 1
    otu_table_txt_desc_handle.close()
    col_annotation_txt_handle.close()
    row_annotation_txt_handle.close()

    # run R script
    pheatmap_cmd = 'Rscript %s -i %s -o %s -r %s -c %s' % (pheatmap_OTU_R, otu_table_txt_desc, output_plot, row_annotation_txt, col_annotation_txt)
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
    blast_parser.add_argument('-o',     required=True,  help='output pdf')
    args = vars(blast_parser.parse_args())
    pheatmap_OTU(args)
