
metadata_txt            = '/Users/songweizhi/Desktop/untitled_folder/metadata_20250713.txt'
status_txt              = '/Users/songweizhi/Desktop/untitled_folder/status.txt'
metadata_txt_updated    = '/Users/songweizhi/Desktop/untitled_folder/metadata_20250713_updated.txt'

sequence_status_dict = dict()
for sample in open(status_txt):
    sample_split = sample.strip().split('\t')
    sequence_status_dict[sample_split[0]] = sample_split[1]

metadata_txt_updated_handle = open(metadata_txt_updated, 'w')
col_index = dict()
line_num_index = 0
num = 0
for each_line in open(metadata_txt):
    line_num_index += 1
    line_split = each_line.strip().split('\t')
    if line_num_index == 1:
        col_index = {key: i for i, key in enumerate(line_split)}
        metadata_txt_updated_handle.write(each_line)
    else:
        sample_id = line_split[col_index['Sample_ID']]
        if sample_id not in sequence_status_dict:
            metadata_txt_updated_handle.write('%s\tna\n' % each_line.strip())
        else:
            metadata_txt_updated_handle.write('%s\t%s\n' % (each_line.strip(), sequence_status_dict[sample_id]))
            num += 1
metadata_txt_updated_handle.close()
