
########################################################################################################################

metadata_txt_JL         = '/Users/songweizhi/Desktop/SMP/01_metadata/metadata_JL.txt'
metadata_txt_shan       = '/Users/songweizhi/Desktop/SMP/01_metadata/shan/Key_plus_Factor_Metadata_brief_mophospecies.txt'
metadata_txt_JL_final   = '/Users/songweizhi/Desktop/SMP/01_metadata/metadata_JL_final.txt'

########################################################################################################################

sample_info_dict_shan = dict()
for each in open(metadata_txt_shan):
    each_split      = each.strip().split('\t')
    sample_id       = each_split[1]
    sample_info_list = each_split[2:]
    sample_p = each_split[11]
    sample_c = each_split[12]
    sample_o = each_split[13]
    sample_sf = each_split[14]
    sample_f = each_split[15]
    sample_g = each_split[16]
    sample_tax_str = 'p__%s;c__%s;o__%s;sf__%s;f__%s;g__%s' % (sample_p,sample_c,sample_o,sample_sf,sample_f,sample_g)
    if sample_tax_str == 'p__na;c__na;o__na;sf__na;f__na;g__na':
        sample_tax_str = 'na'
    sample_info_list.append(sample_tax_str)
    sample_info_dict_shan[sample_id] = sample_info_list


sample_set = set()
metadata_txt_JL_final_handle = open(metadata_txt_JL_final, 'w')
for sample in open(metadata_txt_JL):
    sample_split = sample.strip().split('\t')
    sample_id = sample_split[0]
    sample_set.add(sample_id)
    if sample_id in sample_info_dict_shan:
        sample_info_list = sample_info_dict_shan[sample_id]
        metadata_txt_JL_final_handle.write('%s\t%s\n' % (sample.strip(), '\t'.join(sample_info_list)))
    else:
        print(sample.strip())
metadata_txt_JL_final_handle.close()













for each in sample_info_dict_shan:
    if each not in sample_set:
        print('%s\t%s' % (each, '\t'.join(sample_info_dict_shan[each])))
