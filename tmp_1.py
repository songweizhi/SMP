from Bio import SeqIO


otu_classification_txt      = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s08_AllSamples_unoise_nc.blca.gtdb.2.txt'
otu_classification_txt      = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s08_AllSamples_unoise_nc.blca.silva.2.txt'
otu_seq_file                = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s06_AllSamples_unoise_nc.fasta'
otu_seq_file_unclassified   = '/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB/s06_AllSamples_unoise_nc_unclassified.fasta'


total_num = 0
unclassified_num = 0
unclassified_otu_set = set()
for otu in open(otu_classification_txt):
    otu_split = otu.strip().split('\t')
    otu_id    = otu_split[0]
    otu_tax   = otu_split[1]
    if otu_tax == 'Unclassified':
        unclassified_num += 1
        unclassified_otu_set.add(otu_id)
    total_num += 1

print(total_num)
print(unclassified_num)
print(unclassified_num*100/total_num)

print(len(unclassified_otu_set))
print(unclassified_otu_set)


otu_seq_file_unclassified_handle = open(otu_seq_file_unclassified, 'w')
for each_seq in SeqIO.parse(otu_seq_file, 'fasta'):
    if each_seq.id in unclassified_otu_set:
        otu_seq_file_unclassified_handle.write('>%s\n' % each_seq.id)
        otu_seq_file_unclassified_handle.write('%s\n' % each_seq.seq)
otu_seq_file_unclassified_handle.close()


