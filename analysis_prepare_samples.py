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

metadata_txt        = '/Users/songweizhi/Desktop/SMP/01_metadata/metadata_JL.txt'
fa_dir              = '/Users/songweizhi/Desktop/SMP/00_fa_files'
fa_ext              = 'fna'
fa_dir_separated    = '/Users/songweizhi/Desktop/SMP/02_filtered_sequences_separated'

########################################################################################################################

metadata_dict = dict()
for each_line in open(metadata_txt):
    line_split = each_line.strip().split('\t')
    gnm_id = line_split[0]
    habitat = line_split[1]
    metadata_dict[gnm_id] = habitat


fa_re = '%s/*.%s' % (fa_dir, fa_ext)
fa_list = glob.glob(fa_re)
fa_file_set = set()
for fa_file in fa_list:
    _, _, f_base, _ = sep_path_basename_ext(fa_file)
    fa_file_set.add(f_base)

shared_sample_list, meta_uniq_sample_list, fa_uniq_sample_list = get_shared_uniq_elements(metadata_dict, fa_file_set)

print('shared_sample_list\t%s\t%s'      % (len(shared_sample_list), shared_sample_list))
print('meta_uniq_sample_list\t%s\t%s'   % (len(meta_uniq_sample_list), meta_uniq_sample_list))
print('fa_uniq_sample_list\t%s\t%s'     % (len(fa_uniq_sample_list), fa_uniq_sample_list))
print()

habitat_to_sample_dict = dict()
for sample in shared_sample_list:
    habitat = metadata_dict[sample]
    if habitat not in habitat_to_sample_dict:
        habitat_to_sample_dict[habitat] = set()
    habitat_to_sample_dict[habitat].add(sample)

for habitat in sorted(list(habitat_to_sample_dict.keys())):

    pwd_subdir = '%s/%s' % (fa_dir_separated, habitat)
    sample_set = habitat_to_sample_dict[habitat]
    print('%s\t%s\t%s' % (habitat, len(sample_set), sample_set))

    if os.path.isdir(pwd_subdir) is False:
        os.system('mkdir %s' % pwd_subdir)

    for sample in sample_set:
        os.system('cp %s/%s.%s %s/' % (fa_dir, sample, fa_ext, pwd_subdir))

