import os
import glob
from Bio import SeqIO
########################################################################################################################

metadata_txt            = '/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250415.txt'
ref_tax_txt             = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/COI_ZVF0BP0R016-Alignment_reformatted_ref_accession/accession_organism.txt'
voucher_txt             = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/COI_ZVF0BP0R016-Alignment_reformatted_ref_accession/accession_voucher.txt'
tax_ranks               = 'sc,o,f,g,s'

marker_id               = 'COI'
aln_file                = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/combined_%s_Coral_with_refs_gblocks.aln'           % marker_id
itol_colored_label_txt  = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/iTOL_colored_Label_Coral_taxa_%s_with_refs.txt'    % marker_id

########################################################################################################################

metadata_txt            = '/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250415.txt'
ref_tax_txt             = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/28S_ZPVK91YC016-Alignment_reformatted_ref_accession/accession_organism.txt'
voucher_txt             = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/28S_ZPVK91YC016-Alignment_reformatted_ref_accession/accession_voucher.txt'
tax_ranks               = 'sc,o,f,g,s'

marker_id               = '28S'
aln_file                = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/combined_%s_Coral_with_refs_gblocks.aln'           % marker_id
itol_colored_label_txt  = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/iTOL_colored_Label_Coral_taxa_%s_with_refs.txt'    % marker_id

########################################################################################################################

interested_tax_rank_list = tax_ranks.split(',')

jl_sample_set = set()
sample_tax_dict = dict()
col_index = dict()
line_num_index = 0
for each_line in open(metadata_txt):
    line_num_index += 1
    line_split = each_line.strip().split('\t')
    if line_num_index == 1:
        col_index = {key: i for i, key in enumerate(line_split)}
    else:
        sample_id  = line_split[col_index['Sample_ID']]
        sample_tax = line_split[col_index['Host_Taxonomy_NCBI']]

        sample_tax_interested = sample_tax
        if sample_tax not in ['na', 'no_hit', 'no_hit(Bacterium)']:

            sample_tax_split = sample_tax.split(';')
            sample_tax_list = []
            for rank in sample_tax_split:
                rank_abbrev = rank.split('__')[0]
                if rank_abbrev in interested_tax_rank_list:
                    sample_tax_list.append(rank)
            sample_tax_interested = ';'.join(sample_tax_list)

        sample_tax_dict[sample_id] = sample_tax_interested
        jl_sample_set.add(sample_id)

ref_set = set()
for each_line in open(ref_tax_txt):
    line_split = each_line.strip().split('\t')
    sample_id = line_split[0]
    sample_tax = line_split[1]
    sample_tax_split = sample_tax.split(';')
    sample_tax_list = []
    for rank in sample_tax_split:
        rank_abbrev = rank.split('__')[0]
        if rank_abbrev in interested_tax_rank_list:
            sample_tax_list.append(rank)
    sample_tax_interested = ';'.join(sample_tax_list)
    sample_tax_dict[sample_id] = '%s (%s)' % (sample_tax_interested, sample_id)
    #sample_tax_dict[line_split[0]] = line_split[1]
    ref_set.add(sample_id)

voucher_set = set()
for each_line in open(voucher_txt):
    voucher_set.add(each_line.strip().split()[0])

longest_gnm_id = 0
longest_tax = 0
for sequence in SeqIO.parse(aln_file, 'fasta'):
    seq_id = sequence.id
    gnm_id = seq_id
    if '_COI_' in seq_id:
        gnm_id = seq_id.split('_COI_')[0]
    elif '_28S_' in seq_id:
        gnm_id = seq_id.split('_28S_')[0]
    if len(gnm_id) > longest_gnm_id:
        longest_gnm_id = len(gnm_id)

    gnm_tax = sample_tax_dict.get(gnm_id, gnm_id)
    if len(gnm_tax) > longest_tax:
        longest_tax = len(gnm_tax)

itol_colored_label_txt_handle = open(itol_colored_label_txt, 'w')
itol_colored_label_txt_handle.write('DATASET_TEXT\nSEPARATOR TAB\nDATASET_LABEL\tColoredLabel\nDATA\n')
for sequence in SeqIO.parse(aln_file, 'fasta'):
    seq_id = sequence.id
    gnm_id = seq_id
    if '_COI_' in seq_id:
        gnm_id = seq_id.split('_COI_')[0]
    elif '_28S_' in seq_id:
        gnm_id = seq_id.split('_28S_')[0]
    gnm_tax = sample_tax_dict.get(gnm_id, gnm_id)
    sample_color = '#000000'
    if seq_id in voucher_set:
        sample_color = '#FF0000'
    elif gnm_id in jl_sample_set:
        sample_color = '#00FF00'

    if seq_id in ref_set:
        str_to_write = '%s\t%s' % (seq_id, gnm_tax)
    else:
        str_to_write = '%s\t%s%s (%s)' % (seq_id, gnm_tax, ''*(longest_tax-len(gnm_tax)), seq_id.split('_%s_' % marker_id)[0])
    #itol_colored_label_txt_handle.write(str_to_write + '\n')
    itol_colored_label_txt_handle.write('%s\t-1\t%s\tnormal\t1\t0\n' % (str_to_write, sample_color))
    gnm_tax_split = gnm_tax.split(';')
itol_colored_label_txt_handle.close()
