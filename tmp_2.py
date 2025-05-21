from Bio import SeqIO


metadata_txt                = '/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250319.txt'
barcoding_by_28S_txt        = '/Users/songweizhi/Desktop/SMP/Host_barcoding/barcoding_by_28S.txt'
barcoding_by_COI_txt        = '/Users/songweizhi/Desktop/SMP/Host_barcoding/barcoding_by_COI.txt'
ncbi_tax_by_barcoding_txt   = '/Users/songweizhi/Desktop/SMP/Host_barcoding/ncbi_tax_by_barcoding.txt'
metadata_txt_new            = '/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250323.txt'


sample_note_dict = dict()
col_index = dict()
line_num_index = 0
for each_line in open(metadata_txt):
    line_num_index += 1
    line_split = each_line.strip().split('\t')
    if line_num_index == 1:
        col_index = {key: i for i, key in enumerate(line_split)}
    else:
        sample_id = line_split[col_index['Sample_ID']]
        sample_note = line_split[col_index['Note']]
        sample_note_dict[sample_id] = sample_note


for each_sample in open('/Users/songweizhi/Desktop/coral54.txt'):
    sample_id = each_sample.strip()
    sample_note = sample_note_dict[sample_id]
    print('%s\t%s' % (sample_id, sample_note))