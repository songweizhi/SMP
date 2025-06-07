import os
import glob
from Bio import SeqIO


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


# file_list1 = [os.path.basename(i) for i in glob.glob('/Users/songweizhi/Desktop/SMP/01_fna_files_DY86II_328/*.fna')]
# file_list2 = [os.path.basename(i) for i in glob.glob('/Users/songweizhi/Desktop/SMP/01_fna_files_DY86II_328_new/*.fna')]
# shared_set, list_1_uniq, list_2_uniq = get_shared_uniq_elements(file_list1, file_list2)
# print(sorted(list_1_uniq))
# print(sorted(list_2_uniq))


# file_to_path_dict = dict()
# for each in open('/Users/songweizhi/Desktop/files.txt'):
#     each_split = each.strip().split('/')
#     file_name  = each_split[-1]
#
#     if file_name not in file_to_path_dict:
#         file_to_path_dict[file_name] = []
#     file_to_path_dict[file_name].append(each.strip())


# print(len(file_to_path_dict))
# for each in file_to_path_dict:
#     if len(file_to_path_dict[each]) > 1:
#         for file in file_to_path_dict[each]:
#             j_num = file.split('/')[0].split('-')[-1]
#             file_name_new = '%s_%s.fna' % (file.split('/')[-1][:-4], j_num)
#             print('cp %s.gz %s.gz' % (file, file_name_new))


fa_1        = '/Users/songweizhi/Desktop/SMP/01_fna_files_duplicated/JL315_B01_4A_J006.fna'
fa_2        = '/Users/songweizhi/Desktop/SMP/01_fna_files_duplicated/JL315_B01_4A_J016.fna'
fa_combined = '/Users/songweizhi/Desktop/SMP/01_fna_files_duplicated/JL315_B01_4A.fna'
prefix      =                                                       'JL315_B01_4A'



fa_combined_handle = open(fa_combined, 'w')
seq_index = 1
for seq in SeqIO.parse(fa_1, 'fasta'):
    fa_combined_handle.write('>%s_%s\n' % (prefix, seq_index))
    fa_combined_handle.write('%s\n' % seq.seq)
    seq_index += 1
for seq in SeqIO.parse(fa_2, 'fasta'):
    fa_combined_handle.write('>%s_%s\n'% (prefix, seq_index))
    fa_combined_handle.write('%s\n'% seq.seq)
    seq_index += 1
fa_combined_handle.close()
