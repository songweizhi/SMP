
unclassify_otu_txt = '/Users/songweizhi/Desktop/un.txt'

blast_op_subset    = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB_20250325/s06_AllSamples_unoise_nc_vs_nt_unclassified.txt'

unclassify_otu_set = set()
for each in open(unclassify_otu_txt):
    unclassify_otu_set.add(each.strip())
print(unclassify_otu_set)
print(len(unclassify_otu_set))

blast_op_subset_handle = open(blast_op_subset, 'w')
for each_line in open('/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB_20250325/s06_AllSamples_unoise_nc_vs_nt.txt'):
    q_id = each_line.strip().split('\t')[0]
    if q_id in unclassify_otu_set:
        blast_op_subset_handle.write(each_line)
blast_op_subset_handle.close()
