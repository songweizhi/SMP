from Bio import SeqIO


metadata_txt_in   = '/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250322.txt'
metadata_txt_out  = '/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250324xxx.txt'
ncbi_taxonomy_txt = '/Users/songweizhi/DB/taxdump_20250321/ncbi_taxonomy.txt'


ncbi_taxonomy_dict = dict()
for each in open(ncbi_taxonomy_txt):
    each_split = each.strip().split('\t')
    ncbi_taxonomy_dict[each_split[0]] = each_split[1]

metadata_txt_out_handle = open(metadata_txt_out, 'w')
col_index = dict()
line_num_index = 0
for each_line in open(metadata_txt_in):
    line_num_index += 1
    line_split = each_line.strip().split('\t')
    if line_num_index == 1:
        col_index = {key: i for i, key in enumerate(line_split)}
        metadata_txt_out_handle.write(each_line)
    else:
        sample_id = line_split[col_index['Sample_ID']]
        tax_by_barcoding = line_split[col_index['Taxonomy_by_Barcoding_NCBI']]
        tax_lowest = tax_by_barcoding.split(';')[-1]
        tax_gtdb_format = ncbi_taxonomy_dict.get(tax_lowest, tax_by_barcoding)
        metadata_txt_out_handle.write('%s\t%s\n' % (each_line.strip(), tax_gtdb_format))
metadata_txt_out_handle.close()
