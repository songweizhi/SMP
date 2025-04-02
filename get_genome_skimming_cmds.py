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
sample_id_txt           = '/Users/songweizhi/Desktop/SMP/genome_skimming/sample_id_39.txt'
ref_seq_file_28s        = '/Users/songweizhi/Desktop/SMP/Host_barcoding/28S_crop.fa'
ref_seq_file_coi        = '/Users/songweizhi/Desktop/SMP/Host_barcoding/COI_crop.fa'
assembly_dir            = '/Users/songweizhi/Desktop/SMP/genome_skimming/spades_k55-127_scaffolds'
blast_op_dir            = '/Users/songweizhi/Desktop/SMP/genome_skimming/blast'

# outputs
genome_skimming_28s_fa = '/Users/songweizhi/Desktop/SMP/genome_skimming/Genome_skimming_28S.fasta'
genome_skimming_coi_fa = '/Users/songweizhi/Desktop/SMP/genome_skimming/Genome_skimming_COI.fasta'

########################################################################################################################

if os.path.isfile(genome_skimming_28s_fa):
    os.remove(genome_skimming_28s_fa)
if os.path.isfile(genome_skimming_coi_fa):
    os.remove(genome_skimming_coi_fa)

file_index = 1
for each in open(sample_id_txt):

    sample_id = each.strip()
    fastqc_cmd = 'fastqc %s_1.fq.gz %s_2.fq.gz' % (sample_id, sample_id)
    # print(fastqc_cmd)

    trim_cmd = 'trimmomatic PE -threads 30 -phred33 %s_1.fq.gz %s_2.fq.gz %s_1_P.fq.gz %s_1_UP.fq.gz %s_2_P.fq.gz %s_2_UP.fq.gz ILLUMINACLIP:/scratch/PI/ocessongwz/DB/trimmomatic_adapters/TruSeq3-PE-2.fa:2:30:10 CROP:145 HEADCROP:5 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:25 MINLEN:35' % (sample_id, sample_id, sample_id, sample_id, sample_id, sample_id)
    trim_cmd = 'BioSAK hpc3 -q cpu -a oces -wt 11:59:59 -t 36 -conda mybase2 -n trim%s -c "%s"' % (file_index, trim_cmd)
    # print(trim_cmd)

    fastqc_cmd = 'BioSAK hpc3 -q cpu -a oces -wt 11:59:59 -t 36 -conda mybase2 -n qc_%s -c "fastqc %s_1_P.fq.gz %s_2_P.fq.gz"' % (file_index, sample_id, sample_id)
    # print(fastqc_cmd)

    # combine unpaired reads
    combine_unpaired_reads_cmd = 'gunzip %s_1_UP.fq.gz; gunzip %s_2_UP.fq.gz; cat %s_1_UP.fq %s_2_UP.fq > %s_UP.fq; gzip %s_UP.fq' % (sample_id, sample_id, sample_id, sample_id, sample_id, sample_id)
    # print(combine_unpaired_reads_cmd)

    # spades
    spades_cmd = 'spades.py --meta -t 36 -k 55,77,99,127 -1 %s_1_P.fq.gz -2 %s_2_P.fq.gz -s %s_UP.fq.gz -o %s_spades_k55-127' % (sample_id, sample_id, sample_id, sample_id)
    spades_cmd = 'BioSAK hpc3 -q cpu -a oces -wt 71:59:59 -t 36 -conda mybase2 -n spades_%s -c "%s"' % (file_index, spades_cmd)
    # print(spades_cmd)

    # copy assemblies
    copy_assemblies_cmd = 'cp %s_spades_k55-127/scaffolds.fasta spades_k55-127_scaffolds/%s.fasta' % (sample_id, sample_id)
    # print(copy_assemblies_cmd)

    # compress assemblies
    gzip_cmd = 'gzip %s.fasta' % sample_id
    gzip_cmd = 'BioSAK hpc3 -q cpu-share -wt 11:59:59 -t 36 -conda mybase2 -n gzip%s -c "%s"' % (file_index, gzip_cmd)
    # print(gzip_cmd)

    # mkblastdb
    mkblastdb_cmd = 'makeblastdb -in %s.fasta -dbtype nucl -parse_seqids' % sample_id
    mkblastdb_cmd = 'BioSAK hpc3 -q cpu-share -wt 11:59:59 -t 36 -conda mybase2 -n mkblastdb_%s -c "%s"' % (file_index, mkblastdb_cmd)
    # print(mkblastdb_cmd)

    # blast barcoding gene against assemblies
    blast_cmd_28s = 'blastn -query 28S_crop.fa -db /scratch/PI/ocessongwz/SMP/genome_skimming/blastdb/%s.fasta -out %s_28S.txt -evalue 1e-5 -outfmt 6 -task blastn' % (sample_id, sample_id)
    blast_cmd_coi = 'blastn -query COI_crop.fa -db /scratch/PI/ocessongwz/SMP/genome_skimming/blastdb/%s.fasta -out %s_COI.txt -evalue 1e-5 -outfmt 6 -task blastn' % (sample_id, sample_id)
    blast_cmd_both = '%s; %s' % (blast_cmd_28s, blast_cmd_coi)
    blast_cmd_both = 'BioSAK hpc3 -q cpu-share -wt 11:59:59 -t 36 -conda mybase2 -n blast_%s -c "%s"' % (file_index, blast_cmd_both)
    # print(blast_cmd_both)

    ############################################ extract_barcoding_gene_seq ############################################

    extract_barcoding_gene_seq = True

    assembly_fa  = '%s/%s.fasta'   % (assembly_dir, sample_id)
    blast_op_28s = '%s/%s_28S.txt' % (blast_op_dir, sample_id)
    blast_op_coi = '%s/%s_COI.txt' % (blast_op_dir, sample_id)
    if extract_barcoding_gene_seq is True:
        if os.path.isfile(blast_op_28s):
            get_hit_seq_seq(sample_id, '28S', assembly_fa, ref_seq_file_28s, blast_op_28s, genome_skimming_28s_fa)
        if os.path.isfile(blast_op_coi):
            get_hit_seq_seq(sample_id, 'COI', assembly_fa, ref_seq_file_coi, blast_op_coi, genome_skimming_coi_fa)

    ####################################################################################################################

    file_index += 1


'''

# blast extracted barcoding gene sequences and parse the results
cd /Users/songweizhi/Desktop/SMP/genome_skimming
BioSAK blast -i Genome_skimming_COI_YHRB7PAE013-Alignment.txt -o Genome_skimming_COI_YHRB7PAE013-Alignment_formatted.txt -n 20 -iden 80 -tax /Users/songweizhi/DB/taxdump_20250321/ncbi_taxonomy.txt
BioSAK blast -i Genome_skimming_28S_YHRAZKH4013-Alignment.txt -o Genome_skimming_28S_YHRAZKH4013-Alignment_formatted.txt -n 20 -iden 80 -tax /Users/songweizhi/DB/taxdump_20250321/ncbi_taxonomy.txt


Genome_skimming_COI_YHRB7PAE013-Alignment.txt
Genome_skimming_28S_YHRAZKH4013-Alignment.txt

'''
