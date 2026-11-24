
itol_connection_txt = '/Users/songweizhi/Desktop/aaa_iTOL.txt'

my_dict = dict()
for each in open('/Users/songweizhi/Desktop/aaa.txt'):
    seq_id = each.strip()
    fauna_id = '_'.join(seq_id.split('_')[:2])
    if fauna_id not in my_dict:
        my_dict[fauna_id] = set()
    my_dict[fauna_id].add(seq_id)

itol_connection_txt_handle = open(itol_connection_txt, 'w')
itol_connection_txt_handle.write('DATASET_CONNECTION\nSEPARATOR TAB\nDATASET_LABEL\tdemo_connections\n')
itol_connection_txt_handle.write('COLOR\t#ff0ff0\nDRAW_ARROWS\t0\nARROW_SIZE\t60\nLOOP_SIZE\t100\n')
itol_connection_txt_handle.write('MAXIMUM_LINE_WIDTH\t3\nCURVE_ANGLE\t45\nCENTER_CURVES\t1\nALIGN_TO_LABELS\t1\nDATA\n')
for each_fauna in my_dict:
    seq_set = my_dict[each_fauna]
    if len(seq_set) > 1:
        seq_list = list(seq_set)
        all_combination_set = set()
        for each_seq1 in seq_list:
            for each_seq2 in seq_list:
                if each_seq1 != each_seq2:
                    current_comb = sorted([each_seq1, each_seq2])
                    current_comb_str = '___|___'.join(current_comb)
                    all_combination_set.add(current_comb_str)

        for each_comb in all_combination_set:
            each_comb_split = each_comb.split('___|___')
            seq_1 = each_comb_split[0]
            seq_2 = each_comb_split[1]
            # NODE1,NODE2,WIDTH,COLOR,STYLE,LABEL
            line_to_write = '%s\t%s\t1\t#0000FF\tnormal\t%s' % (seq_1, seq_2, each_fauna)
            itol_connection_txt_handle.write(line_to_write + '\n')
itol_connection_txt_handle.close()


prior_fauna_set = set()
for each in open('/Users/songweizhi/Desktop/10.txt'):
    fauna_id = '_'.join(each.strip().split('_')[:2])
    prior_fauna_set.add(fauna_id)

for each_fauna in my_dict:
    if each_fauna in prior_fauna_set:
        seq_set = my_dict[each_fauna]
        for each in seq_set:
            print(each)




