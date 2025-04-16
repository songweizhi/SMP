import os
from Bio import SeqIO
from Bio.Seq import Seq


def best_hit(file_in, file_out):

    file_out_handle = open(file_out, 'w')
    best_hit_line = ''
    best_hit_query_id = ''
    best_hit_score = 0
    for blast_hit in open(file_in):
        blast_hit_split = blast_hit.strip().split('\t')
        query_id = blast_hit_split[0]
        bit_score = float(blast_hit_split[11])

        if best_hit_query_id == '':
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

        elif (query_id == best_hit_query_id) and (bit_score > best_hit_score):
            best_hit_score = bit_score
            best_hit_line = blast_hit

        elif query_id != best_hit_query_id:
            file_out_handle.write(best_hit_line)
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

    file_out_handle.write(best_hit_line)
    file_out_handle.close()


########################################################################################################################

# file in
barcoding_fa        = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/combined_28S_Coral.fa'
ref_seq_fa          = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/28S_ZPVK91YC016-Alignment_reformatted_ref_accession.fasta'

# file out
blast_op            = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/barcoding_vs_ref_blastn.txt'
blast_op_best_hit   = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/barcoding_vs_ref_blastn_best_hit.txt'
fa_out              = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/28S_ZPVK91YC016-Alignment_reformatted_ref_accession_matched_region.fa'

########################################################################################################################

# file in
barcoding_fa        = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/combined_COI_Coral.fa'
ref_seq_fa          = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/COI_ZVF0BP0R016-Alignment_reformatted_ref_accession/accession_sequence.fasta'

# file out
blast_op            = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/COI_ZVF0BP0R016-Alignment_reformatted_ref_accession/barcoding_vs_ref_blastn.txt'
blast_op_best_hit   = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/COI_ZVF0BP0R016-Alignment_reformatted_ref_accession/barcoding_vs_ref_blastn_best_hit.txt'
fa_out              = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/COI_ZVF0BP0R016-Alignment_reformatted_ref_accession/accession_sequence_matched_region.fa'

########################################################################################################################

seq_dict = dict()
for each_seq in SeqIO.parse(ref_seq_fa, 'fasta'):
    seq_dict[each_seq.id] = str(each_seq.seq)
for each_seq in SeqIO.parse(barcoding_fa, 'fasta'):
    seq_dict[each_seq.id] = str(each_seq.seq)

blast_cmd_28s = 'blastn -query %s -subject %s -out %s -evalue 1e-5 -outfmt 6 -task blastn' % (ref_seq_fa, barcoding_fa, blast_op)
os.system(blast_cmd_28s)

best_hit(blast_op, blast_op_best_hit)

fa_out_handle = open(fa_out, 'w')
for line in open(blast_op_best_hit):
    line_split  = line.strip().split('\t')
    q_id        = line_split[0]
    s_id        = line_split[1]
    q_start     = int(line_split[6])
    q_end       = int(line_split[7])
    s_start     = int(line_split[8])
    s_end       = int(line_split[9])

    q_direction = 1
    if q_end < q_start:
        q_direction = -1

    s_direction = 1
    if s_end < s_start:
        s_direction = -1

    q_seq_matched_region = seq_dict[q_id][(q_start-1) : q_end]
    if q_direction != s_direction:
        q_seq_matched_region = Seq(q_seq_matched_region).reverse_complement()

    fa_out_handle.write('>' + q_id + '\n')
    fa_out_handle.write(str(q_seq_matched_region) + '\n')
fa_out_handle.close()
