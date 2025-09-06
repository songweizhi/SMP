
sample_update_txt = '/Users/songweizhi/Desktop/sample_update_for_Coral_GRF.txt'


taxa_count_dict_grf = dict()
taxa_count_dict_ms = dict()
line_index = 0
for sample in open(sample_update_txt):
    sample_split = sample.strip().split('\t')
    if line_index > 0:
        coral_f = sample_split[4]
        if coral_f not in taxa_count_dict_grf:
            taxa_count_dict_grf[coral_f] = 0
        taxa_count_dict_grf[coral_f] += 1
        if sample_split[1] == 'yes':
            if coral_f not in taxa_count_dict_ms:
                taxa_count_dict_ms[coral_f] = 0
            taxa_count_dict_ms[coral_f] += 1
    line_index += 1

for each_tax in sorted(list(taxa_count_dict_grf.keys())):
    print('%s\t%s\t%s' % (each_tax, taxa_count_dict_grf[each_tax], taxa_count_dict_ms.get(each_tax, 0) ))


