from Bio import SeqIO


fa_in  = '/Users/songweizhi/Desktop/SMP/Host_tree_Sponge/combined_COI_sponge_with_amplicon_plus_ref.fa'
fa_out = '/Users/songweizhi/Desktop/SMP/Host_tree_Sponge/combined_COI_sponge_with_amplicon_plus_ref.fas'


fa_out_handle = open(fa_out, 'w')
for each in SeqIO.parse(fa_in, 'fasta'):
    fa_out_handle.write('>' + each.id + '\n')
    fa_out_handle.write(str(each.seq) + '\n')
fa_out_handle.close()
