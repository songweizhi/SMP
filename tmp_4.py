import os

sample_metadata_txt     = '/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250511.txt'

sample_name_dict = dict()
for each_line in open('/Users/songweizhi/Desktop/SMP/01_fna_files/for_upload/coral_sample_with_barcoding_30_with_Water_Sediment_69_removed_4nodata_1wrong_metadata_20250508_SZ.txt'):
    each_line_split = each_line.strip().split()
    if not each_line.startswith('Sample_ID'):
        raw_name = each_line_split[0]
        new_name = each_line_split[1]
        sample_name_dict[raw_name] = new_name
print(sample_name_dict)

sample_group_dict = dict()
col_index = dict()
line_num_index = 0
for each_line in open(sample_metadata_txt):
    line_num_index += 1
    line_split = each_line.strip().split('\t')
    if line_num_index == 1:
        col_index = {key: i for i, key in enumerate(line_split)}
    else:
        sample_id = line_split[col_index['Sample_ID']]
        sample_source = line_split[col_index['Source']]
        sample_collect_date = line_split[col_index['Collect_Date']]
        Longitude = float(line_split[col_index['Longitude']])
        Latitude = float(line_split[col_index['Latitude']])
        if Longitude is not 'na':
            Longitude = float("{0:.2f}".format(Longitude))  # 123.45
        if Latitude is not 'na':
            Latitude = float("{0:.2f}".format(Latitude))  # 123.45
        Collect_Depth = line_split[col_index['Collect_Depth']]

        if sample_id  in sample_name_dict:

            sample_name_new = sample_name_dict[sample_id]

            print('%s\t%s\tPRJNA1229863\t%s\t%s\tdeep-sea biome\tdeep sea\tPacific Ocean\t%s N %s E\t%s' % (sample_name_new, sample_name_new, sample_source, sample_collect_date.split()[0], Latitude, Longitude, Collect_Depth))



        # sample_source = line_split[col_index['Source']]
        # sample_host_tax_str = line_split[col_index['Host_Taxonomy_NCBI']]
        # sample_host_tax_split = sample_host_tax_str.split(';')
