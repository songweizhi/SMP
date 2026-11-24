import os
from Bio import SeqIO


def get_hit_seq_seq(sample_id, marker, assembly_fa, ref_seq_file, blast_op, extracted_seq_file):

    assembly_seq_dict = dict()
    for each_seq in SeqIO.parse(assembly_fa, 'fasta'):
        assembly_seq_dict[each_seq.id] = str(each_seq.seq)

    seq_len_dict_28s = dict()
    for each_seq in SeqIO.parse(ref_seq_file, 'fasta'):
        seq_len_dict_28s[each_seq.id] = len(each_seq.seq)

    ctg_to_28s_dict = dict()
    for each_hit in open(blast_op):
        each_hit_split = each_hit.strip().split('\t')
        query_id = each_hit_split[0]
        subject_id = each_hit_split[1]
        query_len = seq_len_dict_28s[query_id]
        iden = float(each_hit_split[2])
        aln_len = int(each_hit_split[3])
        aln_cov = float("{0:.3f}".format(aln_len * 100 / query_len))
        sstart = int(each_hit_split[8])
        send = int(each_hit_split[9])
        if (iden >= 80) and (aln_cov >= 80):
            if subject_id not in ctg_to_28s_dict:
                ctg_to_28s_dict[subject_id] = set()
            ctg_to_28s_dict[subject_id].add('%s_%s' % (sorted([sstart, send])[0], sorted([sstart, send])[1]))

    extracted_seq_file_handle = open(extracted_seq_file, 'a')
    seq_index = 1
    for each_ctg in ctg_to_28s_dict:
        hit_list = sorted(list(ctg_to_28s_dict[each_ctg]))
        left_end_list = []
        right_end_list = []
        for seg in hit_list:
            seg_l = int(seg.split('_')[0])
            seg_r = int(seg.split('_')[1])
            left_end_list.append(seg_l)
            right_end_list.append(seg_r)

        if (max(left_end_list) - min(left_end_list) <= 50) and (max(right_end_list) - min(right_end_list) <= 50):
            ctg_seq = assembly_seq_dict[each_ctg]
            hit_seg_seq = ctg_seq[min(left_end_list):max(right_end_list)]
            extracted_seq_file_handle.write('>%s_%s_%s_%s\n' % (sample_id, marker, seq_index, each_ctg))
            extracted_seq_file_handle.write(hit_seg_seq + '\n')
            extracted_seq_file_handle.flush()
            seq_index += 1
        else:
            print('%s\t%s\t%s' % (sample_id, each_ctg, '\t'.join(hit_list)))
    extracted_seq_file_handle.close()


########################################################################################################################

# inputs
ref_seq_file_18s            = '/project/spongeholobiont/Sponge_r226/05_Sponge_tree/RefSeqs_with_AOA_18S.fasta'
ref_seq_file_28s            = '/project/spongeholobiont/Sponge_r226/05_Sponge_tree/RefSeqs_with_AOA_28S.fasta'
ref_seq_file_coi            = '/project/spongeholobiont/Sponge_r226/05_Sponge_tree/RefSeqs_with_AOA_COI_iden95.ffn'
ref_seq_file_18s            = '/Users/songweizhi/Desktop/Sponge_r226/07_Sponge_tree/000/17_sponge_samples_with_AOA/RefSeqs_with_AOA_18S.fasta'
ref_seq_file_28s            = '/Users/songweizhi/Desktop/Sponge_r226/07_Sponge_tree/000/17_sponge_samples_with_AOA/RefSeqs_with_AOA_28S.fasta'
ref_seq_file_coi            = '/Users/songweizhi/Desktop/Sponge_r226/07_Sponge_tree/000/17_sponge_samples_with_AOA/RefSeqs_with_AOA_COI_iden95.ffn'

# 17 sponge samples with AOA
sample_id_txt               = '/Users/songweizhi/Desktop/Sponge_r226/07_Sponge_tree/000/17_sponge_samples_with_AOA/JL_samples_17.txt'
assembly_dir                = '/Users/songweizhi/Desktop/Sponge_r226/07_Sponge_tree/000/17_sponge_samples_with_AOA/JL_sample_17_spades_assembly'
blast_op_dir                = '/Users/songweizhi/Desktop/Sponge_r226/07_Sponge_tree/000/17_sponge_samples_with_AOA/JL_sample_17_spades_assembly_blast_op'
genome_skimming_18s_fa      = '/Users/songweizhi/Desktop/Sponge_r226/07_Sponge_tree/000/17_sponge_samples_with_AOA/18S_17.fasta'
genome_skimming_28s_fa      = '/Users/songweizhi/Desktop/Sponge_r226/07_Sponge_tree/000/17_sponge_samples_with_AOA/28S_17.fasta'
genome_skimming_coi_fa      = '/Users/songweizhi/Desktop/Sponge_r226/07_Sponge_tree/000/17_sponge_samples_with_AOA/COI_17.fasta'
extract_barcoding_gene_seq  = True

########################################################################################################################

if extract_barcoding_gene_seq is True:
    if os.path.isfile(genome_skimming_18s_fa):
        print('%s already exists, program exited!' % genome_skimming_18s_fa)
        exit()
    if os.path.isfile(genome_skimming_28s_fa):
        print('%s already exists, program exited!' % genome_skimming_28s_fa)
        exit()
    if os.path.isfile(genome_skimming_coi_fa):
        print('%s already exists, program exited!' % genome_skimming_coi_fa)
        exit()

file_index = 1
for each in open(sample_id_txt):
    sample_id = each.strip()

    # mkblastdb
    mkblastdb_cmd = 'makeblastdb -in %s.fasta -dbtype nucl -parse_seqids' % sample_id
    mkblastdb_cmd = 'BioSAK hpc4 -q amd -a marmolecol -wt 00:59:59 -t 36 -conda mybase2 -n mkblastdb_%s -c "%s"' % (file_index, mkblastdb_cmd)
    #print(mkblastdb_cmd)

    # blast barcoding gene against assemblies
    blast_cmd_18s = 'blastn -query %s -db %s.fasta -out %s_18S.txt -evalue 1e-5 -outfmt 6 -task blastn'         % (ref_seq_file_18s, sample_id, sample_id)
    blast_cmd_28s = 'blastn -query %s -db %s.fasta -out %s_28S.txt -evalue 1e-5 -outfmt 6 -task blastn'         % (ref_seq_file_28s, sample_id, sample_id)
    blast_cmd_coi = 'blastn -query %s -db %s.fasta -out %s_COI.txt -evalue 1e-5 -outfmt 6 -task blastn'         % (ref_seq_file_coi, sample_id, sample_id)
    blast_cmd_both = '%s; %s; %s'                                                                               % (blast_cmd_18s, blast_cmd_28s, blast_cmd_coi)
    blast_cmd_both = blast_cmd_28s
    blast_cmd_both = 'BioSAK hpc4 -q amd -a marmolecol -wt 00:29:59 -t 36 -conda mybase2 -n blast_%s -c "%s"'   % (file_index, blast_cmd_both)
    print(blast_cmd_both)

    # extract sequences
    assembly_fa  = '%s/%s.fasta'   % (assembly_dir, sample_id)
    blast_op_18s = '%s/%s_18S.txt' % (blast_op_dir, sample_id)
    blast_op_28s = '%s/%s_28S.txt' % (blast_op_dir, sample_id)
    blast_op_coi = '%s/%s_COI.txt' % (blast_op_dir, sample_id)
    if extract_barcoding_gene_seq is True:
        if os.path.isfile(blast_op_18s):
            get_hit_seq_seq(sample_id, '18S', assembly_fa, ref_seq_file_18s, blast_op_18s, genome_skimming_18s_fa)
        if os.path.isfile(blast_op_28s):
            get_hit_seq_seq(sample_id, '28S', assembly_fa, ref_seq_file_28s, blast_op_28s, genome_skimming_28s_fa)
        if os.path.isfile(blast_op_coi):
            get_hit_seq_seq(sample_id, 'COI', assembly_fa, ref_seq_file_coi, blast_op_coi, genome_skimming_coi_fa)

    # then run online blast for the extracted sequences and download the results (select the "Text" format)
    # then parse the downloaded blast results with "BioSAK blast"

    ####################################################################################################################

    file_index += 1


'''

# blast extracted barcoding gene sequences and parse the results
cd /Users/songweizhi/Desktop/SMP/Genome_Skimming_and_Metagenome
BioSAK blast -i Genome_skimming_28S_unclassified_7_FCMHVNK9014-Alignment.txt -o Genome_skimming_28S_unclassified_7_FCMHVNK9014-Alignment_formatted.txt -n 20 -iden 80 -tax /Users/songweizhi/DB/taxdump_20250321/ncbi_taxonomy.txt
BioSAK blast -i Genome_skimming_COI_unclassified_7_FCMJKGMR015-Alignment.txt -o Genome_skimming_COI_unclassified_7_FCMJKGMR015-Alignment_formatted.txt -n 20 -iden 80 -tax /Users/songweizhi/DB/taxdump_20250321/ncbi_taxonomy.txt

cd /Users/songweizhi/Desktop/SMP/Genome_Skimming_and_Metagenome
BioSAK blast -i Genome_skimming_28S_batch4_MW498N6V016-Alignment.txt -o Genome_skimming_28S_batch4_MW498N6V016-Alignment_formatted.txt -n 20 -iden 80 -tax /Users/songweizhi/DB/taxdump_20250321/ncbi_taxonomy.txt
BioSAK blast -i Genome_skimming_COI_batch4_MW49JA26016-Alignment.txt -o Genome_skimming_COI_batch4_MW49JA26016-Alignment_formatted.txt -n 20 -iden 80 -tax /Users/songweizhi/DB/taxdump_20250321/ncbi_taxonomy.txt

cd /Users/songweizhi/Desktop/888
BioSAK blast -i MWFTM9N2014-Alignment.txt -o MWFTM9N2014-Alignment_formatted.txt -n 5 -iden 80 -tax /Users/songweizhi/DB/taxdump_20250321/ncbi_taxonomy.txt
BioSAK blast -i MWFTZADS014-Alignment.txt -o MWFTZADS014-Alignment_formatted.txt -n 5 -iden 80 -tax /Users/songweizhi/DB/taxdump_20250321/ncbi_taxonomy.txt

cd /Users/songweizhi/Desktop/Sponge_r226/07_Sponge_tree/17_sponge_samples_with_AOA
BioSAK blast -i 18S_17-UGEKG5XG014-Alignment.txt -o 18S_17-UGEKG5XG014-Alignment_formatted.txt -n 10 -iden 75 -tax /Users/songweizhi/DB/taxdump_20250321/ncbi_taxonomy.txt
BioSAK blast -i COI_17-UGEECZSE014-Alignment.txt -o COI_17-UGEECZSE014-Alignment_formatted.txt -n 10 -iden 75 -tax /Users/songweizhi/DB/taxdump_20250321/ncbi_taxonomy.txt

cd /Users/songweizhi/Desktop/Sponge_r226/07_Sponge_tree/000/17_sponge_samples_with_AOA
BioSAK blast -i 28S_17_V082P3DE016-Alignment.txt -o 28S_17_V082P3DE016-Alignment_formatted.txt -n 10 -iden 75 -tax /Users/songweizhi/DB/taxdump_20250321/ncbi_taxonomy.txt

'''
