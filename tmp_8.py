from Bio import SeqIO


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


set_with_seq = set()
for each_seq in SeqIO.parse('/Users/songweizhi/Desktop/SMP/Host_tree_Sponge/combined_COI_sponge_with_amplicon.fa', 'fasta'):
    seq_id = each_seq.id
    if '_COI_' in seq_id:
        seq_id = seq_id.split('_COI_')[0]
    set_with_seq.add(seq_id)

set_45 = set()
for each in open('/Users/songweizhi/Desktop/45.txt'):
    bin_id = each.strip().split()[0]
    set_45.add(bin_id)
print(len(set_45))

shared_set, set_45_uniq, set_with_seq_uniq = get_shared_uniq_elements(set_45, set_with_seq)
print(len(shared_set), shared_set)
print(len(set_45_uniq), sorted(list(set_45_uniq)))
print(len(set_with_seq_uniq), sorted(list(set_with_seq_uniq)))

