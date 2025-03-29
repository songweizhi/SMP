from Bio import SeqIO


metadata_txt                = '/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250319.txt'
barcoding_by_28S_txt        = '/Users/songweizhi/Desktop/SMP/Host_barcoding/barcoding_by_28S.txt'
barcoding_by_COI_txt        = '/Users/songweizhi/Desktop/SMP/Host_barcoding/barcoding_by_COI.txt'
ncbi_tax_by_barcoding_txt   = '/Users/songweizhi/Desktop/SMP/Host_barcoding/ncbi_tax_by_barcoding.txt'
metadata_txt_new            = '/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250323.txt'

barcoding_by_both_set = set()
barcoding_by_28S_set = set()
for each in open(barcoding_by_28S_txt):
    barcoding_by_28S_set.add(each.strip())
    barcoding_by_both_set.add(each.strip())

barcoding_by_COI_set = set()
for each in open(barcoding_by_COI_txt):
    barcoding_by_COI_set.add(each.strip())
    barcoding_by_both_set.add(each.strip())

barcoding_dict = dict()
iden_dict = dict()
for each in open(ncbi_tax_by_barcoding_txt):
    each_split = each.strip().split('\t')
    sample_id = each_split[0]
    tax = each_split[1]
    tax_str = tax
    iden = 'na'
    if '(' in tax:
        tax_str = tax.split('(')[0]
        iden = tax.split('(')[1][:-1]
    barcoding_dict[sample_id] = tax_str
    iden_dict[sample_id] = iden


barcoded_sample_set = set()

metadata_txt_new_handle = open(metadata_txt_new, 'w')
# read in metadata_txt
sample_group_dict = dict()
col_index = dict()
line_num_index = 0
for each_line in open(metadata_txt):
    line_num_index += 1
    line_split = each_line.strip().split('\t')
    if line_num_index == 1:
        col_index = {key: i for i, key in enumerate(line_split)}
        metadata_txt_new_handle.write(each_line)
    else:
        sample_id = line_split[col_index['Sample_ID']]

        barcoding_gene = 'na'
        if sample_id in barcoding_by_28S_set:
            barcoding_gene = '28S'
        elif sample_id in barcoding_by_COI_set:
            barcoding_gene = 'COI'

        tax_by_barcoding = barcoding_dict.get(sample_id, 'na')
        identity = iden_dict.get(sample_id, 'na')

        if tax_by_barcoding != 'na':
            barcoded_sample_set.add(sample_id)

        metadata_txt_new_handle.write('%s\t%s\t%s\t%s\n' % (each_line.strip(), barcoding_gene, tax_by_barcoding , identity))

metadata_txt_new_handle.close()




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




shared_set, list_1_uniq, list_2_uniq = get_shared_uniq_elements(barcoding_by_both_set, barcoded_sample_set)


print(shared_set)
print(list_1_uniq)
print(list_2_uniq)
