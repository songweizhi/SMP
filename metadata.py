import os
import glob


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


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


########################################################################################################################

metadata_txt    = '/Users/songweizhi/Desktop/SMP/01_metadata/metadata_JL.txt'
fa_dir          = '/Users/songweizhi/Desktop/SMP/00_fa_files'
fa_ext          = 'fna'

########################################################################################################################

# get sample list with sequence file
fa_re        = '%s/*.%s' % (fa_dir, fa_ext)
fa_list      = glob.glob(fa_re)
fa_list_base = [('.'.join(i.split('/')[-1].split('.')[:-1])) for i in fa_list]

# read in metadata_txt
col_index = dict()
line_num_index = 0
for each_line in open(metadata_txt):
    line_num_index += 1
    line_split = each_line.strip().split('\t')
    if line_num_index == 1:
        col_index = {key: i for i, key in enumerate(line_split)}
    else:
        sample_id       = line_split[col_index['Sample']]
        sample_source   = line_split[col_index['Source']]
        sample_dive     = line_split[col_index['Dive']]

        if sample_source in ['Sponge', 'Sediment', 'Water']:
            if sample_id in fa_list_base:
                print('%s\t%s' % (sample_id, sample_source))



