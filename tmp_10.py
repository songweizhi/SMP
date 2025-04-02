from Bio import SeqIO
from MarkerMAG.matam_16s import str_to_num_list

sample_txt = '/Users/songweizhi/Desktop/coral184.txt'

txt_1       = '/Users/songweizhi/Desktop/SMP/genome_skimming/Genome_skimming_28S_YHRAZKH4013-Alignment_formatted.txt'
txt_1_out   = '/Users/songweizhi/Desktop/SMP/genome_skimming/Genome_skimming_28S_YHRAZKH4013-Alignment_formatted_Coral.txt'

# txt_1       = '/Users/songweizhi/Desktop/SMP/genome_skimming/Genome_skimming_COI_YHRB7PAE013-Alignment_formatted.txt'
# txt_1_out   = '/Users/songweizhi/Desktop/SMP/genome_skimming/Genome_skimming_COI_YHRB7PAE013-Alignment_formatted_Coral.txt'


sample_set = set()
for sample in open(sample_txt):
    sample_set.add(sample.strip())


txt_1_out_handle = open(txt_1_out, 'w')
for each_line in open(txt_1):
    each_line_split = each_line.strip().split('\t')
    sample_id = each_line_split[0].split('_28S_')[0]
    if sample_id not in sample_set:
        txt_1_out_handle.write(each_line)
txt_1_out_handle.close()


# for each_seq in SeqIO.parse('/Users/songweizhi/Desktop/SMP/Host_barcoding/28S_crop_Coral.fa', 'fasta'):
#     if each_seq.id in sample_set:
#         print('>%s' % each_seq.id)
#         print('%s' % each_seq.seq)
#
#
#     # if each_seq.id.split('_COI_')[0] in sample_set:
#     #     print('>%s' % each_seq.id)
#     #     print('%s' % each_seq.seq)


sample_tax_dict = dict()
col_index = dict()
line_num_index = 0
for each_line in open("/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250402.txt"):
    line_num_index += 1
    line_split = each_line.strip().split('\t')
    if line_num_index == 1:
        col_index = {key: i for i, key in enumerate(line_split)}
    else:
        sample_id  = line_split[col_index['Sample_ID']]
        sample_tax = line_split[col_index['Host_Taxonomy_NCBI']]
        sample_tax_dict[sample_id] = sample_tax


longest_gnm_id = 0
longest_tax = 0
for sequence in SeqIO.parse('/Users/songweizhi/Desktop/SMP/Host_tree/combined_COI_Coral_gblocks.aln', 'fasta'):
    seq_id = sequence.id
    gnm_id = seq_id
    if '_COI_' in seq_id:
        gnm_id = seq_id.split('_COI_')[0]
    elif '_28S_' in seq_id:
        gnm_id = seq_id.split('_28S_')[0]
    if len(gnm_id) > longest_gnm_id:
        longest_gnm_id = len(gnm_id)

    gnm_tax = sample_tax_dict[gnm_id]

    if len(gnm_tax) > longest_tax:
        longest_tax = len(gnm_tax)

handle = open('/Users/songweizhi/Desktop/SMP/Host_tree/iTOL_Coral_taxa.txt', 'w')
handle.write('LABELS\nSEPARATOR TAB\nDATA\n')
for sequence in SeqIO.parse('/Users/songweizhi/Desktop/SMP/Host_tree/combined_COI_Coral_gblocks.aln', 'fasta'):
    seq_id = sequence.id
    gnm_id = seq_id
    if '_COI_' in seq_id:
        gnm_id = seq_id.split('_COI_')[0]
    elif '_28S_' in seq_id:
        gnm_id = seq_id.split('_28S_')[0]
    gnm_tax = sample_tax_dict[gnm_id]
    print('%s\t%s%s__%s' % (seq_id, gnm_tax, '_'*(longest_tax-len(gnm_tax)), seq_id))

    str_to_write = '%s\t%s%s__%s' % (seq_id, gnm_tax, '_'*(longest_tax-len(gnm_tax)), seq_id)
    str_to_write = str_to_write.replace('d__Eukaryota;k__Metazoa;p__Cnidaria;c__Anthozoa;', '')
    str_to_write = str_to_write.replace('d__Eukaryota;k__Metazoa;p__Cnidaria;', '')
    handle.write(str_to_write + '\n')
handle.close()













