import os
import glob
from Bio import SeqIO
from MarkerMAG.matam_16s import str_to_num_list

########################################################################################################################

metadata_txt            = '/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250408.txt'
samples_to_ignore_txt   = '/Users/songweizhi/Desktop/SMP/Coral_samples_to_ignore_51_Shan.txt'
samples_to_ignore_txt   = '/Users/songweizhi/Desktop/SMP/Coral_samples_to_ignore_66.txt'
sample_txt              = '/Users/songweizhi/Desktop/SMP/coral_sample_with_barcoding_30_with_Water_Sediment_69.txt'
amplicon_dir            = '/Users/songweizhi/Desktop/SMP/01_fna_files'
amplicon_ext            = 'fna'
ref_tax_txt             = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/28S_ZPVK91YC016-Alignment_reformatted_ref_taxon.txt'
ref_tax_txt             = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/COI_ZVF0BP0R016-Alignment_reformatted_ref_accession/accession_organism.txt'

marker_id               = '28S'
marker_id               = 'COI'
aln_file                = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/combined_%s_Coral_gblocks.aln'                     % marker_id
aln_file                = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/combined_%s_Coral_with_refs_gblocks.aln'           % marker_id
itol_label_txt          = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/iTOL_Label_Coral_taxa_%s_with_refs.txt'            % marker_id
itol_colored_label_txt  = '/Users/songweizhi/Desktop/SMP/Host_tree_Coral/iTOL_colored_Label_Coral_taxa_%s_with_refs.txt'    % marker_id

########################################################################################################################

amplicon_file_re = '%s/*.%s' % (amplicon_dir, amplicon_ext)
amplicon_file_list = [os.path.basename(i) for i in glob.glob(amplicon_file_re)]
amplicon_id_list = ['.'.join(i.split('.')[:-1]) for  i in amplicon_file_list]

to_ignore_sample_set = set()
for sample in open(samples_to_ignore_txt):
    to_ignore_sample_set.add(sample.strip().split()[0])

sample_set = set()
for sample in open(sample_txt):
    sample_set.add(sample.strip())

sample_depth_dict = dict()
sample_tax_dict = dict()
sample_source_dict = dict()
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
        sample_tax_dict[sample_id] = sample_tax
        sample_source = line_split[col_index['Source']]
        sample_depth = line_split[col_index['Collect_Depth']]
        sample_source_dict[sample_id] = sample_source
        sample_depth_dict[sample_id] = sample_depth

for each_line in open(ref_tax_txt):
    line_split = each_line.strip().split('\t')
    sample_tax_dict[line_split[0]] = line_split[1]


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
    print(gnm_tax)
    if len(gnm_tax) > longest_tax:
        longest_tax = len(gnm_tax)



handle = open(itol_label_txt, 'w')
handle.write('LABELS\nSEPARATOR TAB\nDATA\n')
for sequence in SeqIO.parse(aln_file, 'fasta'):
    seq_id = sequence.id
    gnm_id = seq_id
    if '_COI_' in seq_id:
        gnm_id = seq_id.split('_COI_')[0]
    elif '_28S_' in seq_id:
        gnm_id = seq_id.split('_28S_')[0]
    gnm_tax = sample_tax_dict.get(gnm_id, gnm_id)
    #print('%s\t%s%s__%s' % (seq_id, gnm_tax, '_'*(longest_tax-len(gnm_tax)), seq_id))

    gnm_source = sample_source_dict.get(gnm_id, gnm_id)

    str_to_write = '%s\t%s%s__%s' % (seq_id, gnm_tax, '_'*(longest_tax-len(gnm_tax)), seq_id.split('_%s_' % marker_id)[0])
    str_to_write = str_to_write.replace('d__Eukaryota;k__Metazoa;p__Cnidaria;c__Anthozoa;', '')
    str_to_write = str_to_write.replace('d__Eukaryota;k__Metazoa;p__Cnidaria;c__Hydrozoa;o__;f__;g__;s__', 'c__Hydrozoa')
    handle.write(str_to_write + '\n')

    gnm_tax_split = gnm_tax.split(';')
    #print('%s\t%s' % (seq_id, '\t'.join(gnm_tax_split)))

handle.close()



