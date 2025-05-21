from Bio import SeqIO

coral_sample_id_txt = '/Users/songweizhi/Desktop/SMP/coral_sample_with_barcoding_34.txt'
combined_28S_raw_fa = '/Users/songweizhi/Desktop/SMP/Host_tree_all/combined_28S_raw.fasta'
combined_COI_raw_fa = '/Users/songweizhi/Desktop/SMP/Host_tree_all/combined_COI_raw.fasta'


coral_sample_set = set()
for coral_sample in open(coral_sample_id_txt):
    coral_sample_set.add(coral_sample.strip())
print(coral_sample_set)
print(len(coral_sample_set))


for seq in SeqIO.parse(combined_COI_raw_fa, 'fasta'):
    seq_id = seq.id
    sample_id = seq_id
    if '_28S_' in seq_id:
        sample_id = seq_id.split('_28S_')[0]
    if '_COI_' in seq_id:
        sample_id = seq_id.split('_COI_')[0]

    if sample_id not in coral_sample_set:
        print('>%s' % seq_id)
        print('%s' % seq.seq)
