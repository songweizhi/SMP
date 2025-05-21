from Bio import SeqIO


# fa_28s          = '/Users/songweizhi/Desktop/28S.fa'
# fa_28s_renamed  = '/Users/songweizhi/Desktop/28S_renamed.fa'

# fa_coi          = '/Users/songweizhi/Desktop/COI.fa'
# fa_coi_renamed  = '/Users/songweizhi/Desktop/COI_renamed.fa'


# fa_28s_renamed_handle = open(fa_28s_renamed, 'w')
# for each in SeqIO.parse(fa_28s, 'fasta'):
#     seq_id = '_'.join(each.id.split('-')[3:])
#     fa_28s_renamed_handle.write('>%s\n' % seq_id)
#     fa_28s_renamed_handle.write('%s\n' % str(each.seq))
# fa_28s_renamed_handle.close()


# fa_coi_renamed_handle = open(fa_coi_renamed, 'w')
# for each in SeqIO.parse(fa_coi, 'fasta'):
#     seq_id = '_'.join(each.id.split('-')[3:])
#     seq_str = str(each.seq)
#     fa_coi_renamed_handle.write('>%s\n' % seq_id)
#     fa_coi_renamed_handle.write('%s\n' % seq_str)
# fa_coi_renamed_handle.close()


fa_28s          = '/Users/songweizhi/Desktop/SMP/Host_barcoding/28S.fa'
fa_28s_cropped  = '/Users/songweizhi/Desktop/SMP/Host_barcoding/28S_crop.fa'

fa_coi          = '/Users/songweizhi/Desktop/SMP/Host_barcoding/COI_checked_rc.fa'
fa_coi_cropped  = '/Users/songweizhi/Desktop/SMP/Host_barcoding/COI_crop.fa'




fa_coi_cropped_handle = open(fa_coi_cropped, 'w')
for each in SeqIO.parse(fa_coi, 'fasta'):
    seq_str = str(each.seq)
    seq_str_cropped = seq_str[20:-20]
    fa_coi_cropped_handle.write('>%s\n' % each.id)
    fa_coi_cropped_handle.write('%s\n' % seq_str_cropped)
fa_coi_cropped_handle.close()


fa_28s_cropped_handle = open(fa_28s_cropped, 'w')
for each in SeqIO.parse(fa_28s, 'fasta'):
    seq_str = str(each.seq)
    seq_str_cropped = seq_str[20:-20]
    fa_28s_cropped_handle.write('>%s\n' % each.id)
    fa_28s_cropped_handle.write('%s\n' % seq_str_cropped)
fa_28s_cropped_handle.close()

