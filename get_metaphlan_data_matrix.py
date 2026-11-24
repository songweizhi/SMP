import os
import glob
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


########################################################################################################################

metaphlan_op_dir    = '/Users/songweizhi/Desktop/metaphlan_results'           # the folder where you put metaphlan outputs for all samples
metaphlan_op_ext    = 'tsv'                                                   # metaphlan output file extension, should be tsv
interested_rank     = 't'                                                     # summarize reads number at species level, you can also try t, g, f, etc.
op_data_matrix_txt  = '/Users/songweizhi/Desktop/metaphlan_results_df_t.txt'  # the output data matrix
divide_read_num_by  = 3000                                                    # for debugging/testing, Keep this value as 1 in your formal analysis

########################################################################################################################

op_data_matrix_txt_tmp  = '/Users/songweizhi/Desktop/metaphlan_op_df_tmp.txt'

metaphlan_op_re = '%s/*.%s' % (metaphlan_op_dir, metaphlan_op_ext)
metaphlan_op_list = glob.glob(metaphlan_op_re)

sample_stats_dod = dict()
all_identified_taxon_set = set()
for each_op_file in metaphlan_op_list:
    _, _, sample_id, _ = sep_path_basename_ext(each_op_file)
    current_sample_stats_dict = dict()
    for each_line in open(each_op_file):
        if not each_line.startswith('#'):
            each_line_split = each_line.strip().split('\t')
            taxon_str = each_line_split[0]
            read_num = int(each_line_split[4])
            if each_line.startswith('UNCLASSIFIED'):
                current_sample_stats_dict['UNCLASSIFIED'] = read_num
                all_identified_taxon_set.add('UNCLASSIFIED')
            else:
                interested_rank_taxon = ''
                for each_rank in taxon_str.split('|'):
                    if each_rank.startswith('%s__' % interested_rank):
                        interested_rank_taxon = each_rank
                if interested_rank_taxon != '':
                    if interested_rank_taxon not in current_sample_stats_dict:
                        current_sample_stats_dict[interested_rank_taxon] = 0
                    current_sample_stats_dict[interested_rank_taxon] += read_num
                    all_identified_taxon_set.add(interested_rank_taxon)
    sample_stats_dod[sample_id] = current_sample_stats_dict

all_identified_taxon_list_sorted = sorted(list(all_identified_taxon_set))

# write out to file
op_data_matrix_txt_tmp_handle = open(op_data_matrix_txt_tmp, 'w')
op_data_matrix_txt_tmp_handle.write('\t%s\n' % '\t'.join (all_identified_taxon_list_sorted))
for each_sample in sorted(list(sample_stats_dod.keys())):
    current_count_list = [each_sample]
    for each_taxon in all_identified_taxon_list_sorted:
        current_count_list.append(sample_stats_dod[each_sample].get(each_taxon, 0))
    current_count_list_divided_by_10000 = [int(i/divide_read_num_by//1) for i in current_count_list[1:]]
    current_count_list_divided_by_10000.insert(0, each_sample)
    current_count_list_as_str = [str(i) for i in current_count_list_divided_by_10000]
    op_data_matrix_txt_tmp_handle.write('%s\n' % '\t'.join (current_count_list_as_str))
op_data_matrix_txt_tmp_handle.close()

transpose_csv(op_data_matrix_txt_tmp, op_data_matrix_txt, '\t', 0, 0)

os.system('rm %s' % op_data_matrix_txt_tmp)
print('Done!')
