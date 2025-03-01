import os

########################################################################################################################

blast_op_txt    = '/Users/songweizhi/Desktop/unclassified_against_nt_best_hit.txt'
op_dir          = '/Users/songweizhi/Desktop/unclassified_against_nt'
otu_table       = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s07_AllSamples_unoise_otu_table.txt'
otu_table_nonEU = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s07_AllSamples_unoise_otu_table_nonEU.txt'

########################################################################################################################

ref_to_otu_dict = dict()
matched_ref_set = set()
for each_line in open(blast_op_txt):
    each_split = each_line.strip().split('\t')
    query_id   = each_split[0]
    subject_id = each_split[1]
    matched_ref_set.add(subject_id)
    if subject_id not in ref_to_otu_dict:
        ref_to_otu_dict[subject_id] = set()
    ref_to_otu_dict[subject_id].add(query_id)

empty_file_num = 0
eu_ref_set = set()
for each_id in matched_ref_set:
    op_txt = '%s/%s.txt' % (op_dir, each_id)
    cmd = 'esearch -db nucleotide -query %s | efetch -format docsum | xtract -pattern DocumentSummary -element TaxId | efetch -db taxonomy -format xml | xtract -pattern Taxon -block "*/Taxon" -sep ":" -element Rank,ScientificName > %s' % (each_id, op_txt)
    #os.system(cmd)

    first_line = open(op_txt).readline().strip()

    if len(first_line) == 0:
        print(cmd)
        empty_file_num += 1

    if 'superkingdom:Eukaryota' in first_line:
        eu_ref_set.add(each_id)

eu_otu_set = set()
for each_ref in eu_ref_set:
    matched_otu_set = ref_to_otu_dict.get(each_ref)
    for matched_otu in matched_otu_set:
        eu_otu_set.add(matched_otu)

otu_table_nonEU_handle = open(otu_table_nonEU, 'w')
for each_otu in open(otu_table):
    each_split = each_otu.strip().split('\t')
    otu_id = each_split[0]
    if otu_id not in eu_otu_set:
        otu_table_nonEU_handle.write(each_otu)
otu_table_nonEU_handle.close()


print('empty_file_num: %s' % empty_file_num)
print('eu_otu_set(%s): %s' % (len(eu_otu_set), ','.join(sorted(list(eu_otu_set)))))
