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
ref_seq_file_28s        = '/Users/songweizhi/Desktop/SMP/Host_barcoding_Sanger/28S_crop.fa'
ref_seq_file_coi        = '/Users/songweizhi/Desktop/SMP/Host_barcoding_Sanger/COI_crop.fa'


# 39
sample_id_txt           = '/Users/songweizhi/Desktop/SMP/genome_skimming/sample_id_39.txt'
assembly_dir            = '/Users/songweizhi/Desktop/SMP/genome_skimming/spades_k55-127_scaffolds'
blast_op_dir            = '/Users/songweizhi/Desktop/SMP/genome_skimming/spades_k55-127_scaffolds_blast_op'
genome_skimming_28s_fa  = '/Users/songweizhi/Desktop/SMP/genome_skimming/Genome_skimming_28S_39.fasta'
genome_skimming_coi_fa  = '/Users/songweizhi/Desktop/SMP/genome_skimming/Genome_skimming_COI_39.fasta'


# Maeva 11
sample_id_txt          = '/Users/songweizhi/Desktop/SMP/genome_skimming/sample_id_maeva11.txt'
assembly_dir            = '/Users/songweizhi/Desktop/SMP/genome_skimming/Mito_contigs_Maeva'
blast_op_dir            = '/Users/songweizhi/Desktop/SMP/genome_skimming/Mito_contigs_Maeva_blast_op'
genome_skimming_28s_fa  = '/Users/songweizhi/Desktop/SMP/genome_skimming/Genome_skimming_28S_maeva11.fasta'
genome_skimming_coi_fa  = '/Users/songweizhi/Desktop/SMP/genome_skimming/Genome_skimming_COI_maeva11.fasta'


# batch2 13
sample_id_txt               = '/Users/songweizhi/Desktop/SMP/Genome_Skimming/spades_k55-127_scaffolds_batch2_13.txt'
assembly_dir                = '/Users/songweizhi/Desktop/SMP/Genome_Skimming/spades_k55-127_scaffolds_batch2_13'
blast_op_dir                = '/Users/songweizhi/Desktop/SMP/genome_skimming/spades_k55-127_scaffolds_batch2_13_blast_op'
genome_skimming_28s_fa      = '/Users/songweizhi/Desktop/SMP/genome_skimming/Genome_skimming_28S_batch2_13.fasta'
genome_skimming_coi_fa      = '/Users/songweizhi/Desktop/SMP/genome_skimming/Genome_skimming_COI_batch2_13.fasta'
extract_barcoding_gene_seq  = True



sample_id_txt          = '/Users/songweizhi/Desktop/SMP/Genome_Skimming/sample_id_13_2025_07_30.txt'

# Changyu Zhu
sample_id_txt          = '/Users/songweizhi/Desktop/Ophiuroidea/Ophiuroidea_samples.txt'
'''


mv 海蛇尾6gut_1.fq.gz Ophiuroidea_gut6_1.fq.gz
mv 海蛇尾6gut_2.fq.gz Ophiuroidea_gut6_2.fq.gz

'''

sample_id_txt           = '/Users/songweizhi/Desktop/SMP/Genome_Skimming/spades_k55-127_scaffolds_batch2_13.txt'


########################################################################################################################

if os.path.isfile(genome_skimming_28s_fa):
    os.remove(genome_skimming_28s_fa)
if os.path.isfile(genome_skimming_coi_fa):
    os.remove(genome_skimming_coi_fa)

file_index = 1
for each in open(sample_id_txt):

    sample_id = each.strip()
    fastqc_cmd = 'fastqc %s_1.fq.gz %s_2.fq.gz' % (sample_id, sample_id)
    fastqc_cmd = 'BioSAK hpc3 -q cpu -a oces -wt 11:59:59 -t 36 -conda mybase2 -n fastqc%s -c "%s"' % (file_index, fastqc_cmd)
    #print(fastqc_cmd)

    trim_cmd = 'trimmomatic PE -threads 30 -phred33 %s_1.fq.gz %s_2.fq.gz %s_1_P.fq.gz %s_1_UP.fq.gz %s_2_P.fq.gz %s_2_UP.fq.gz ILLUMINACLIP:/scratch/PI/ocessongwz/DB/trimmomatic_adapters/TruSeq3-PE-2.fa:2:30:10 CROP:145 HEADCROP:5 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:25 MINLEN:35' % (sample_id, sample_id, sample_id, sample_id, sample_id, sample_id)
    trim_cmd = 'BioSAK hpc3 -q cpu -a boqianpy -wt 11:59:59 -t 36 -conda mybase2 -n trim%s -c "%s"' % (file_index, trim_cmd)
    #print(trim_cmd)

    fastqc_cmd = 'BioSAK hpc3 -q cpu -a oces -wt 11:59:59 -t 36 -conda mybase2 -n qc_%s -c "fastqc %s_1_P.fq.gz %s_2_P.fq.gz"' % (file_index, sample_id, sample_id)
    #print(fastqc_cmd)

    # combine unpaired reads
    combine_unpaired_reads_cmd = 'gunzip %s_1_UP.fq.gz; gunzip %s_2_UP.fq.gz; cat %s_1_UP.fq %s_2_UP.fq > %s_UP.fq; gzip %s_UP.fq' % (sample_id, sample_id, sample_id, sample_id, sample_id, sample_id)
    combine_unpaired_reads_cmd = 'BioSAK hpc3 -q cpu -a oces -wt 11:59:59 -t 36 -conda mybase2 -n cat_up_%s -c "%s"' % (file_index, combine_unpaired_reads_cmd)
    #print(combine_unpaired_reads_cmd)

    # spades
    spades_cmd = 'spades.py --meta -t 36 -k 55,77,99,127 -1 %s_1_P.fq.gz -2 %s_2_P.fq.gz -s %s_UP.fq.gz -o %s_spades_k55-127' % (sample_id, sample_id, sample_id, sample_id)
    spades_cmd = 'BioSAK hpc3 -q cpu -a oces -wt 71:59:59 -t 36 -conda mybase2 -n spades_%s -c "%s"' % (file_index, spades_cmd)
    #print(spades_cmd)

    # copy and rename spades assemblies
    copy_assemblies_cmd = 'cp %s_spades_k55-127/scaffolds.fasta spades_k55-127_scaffolds/%s.fasta' % (sample_id, sample_id)
    #print(copy_assemblies_cmd)

    # compress spades_wd
    tar_cmd = 'tar -czvf %s_spades_k55-127.tar.gz %s_spades_k55-127'                            % (sample_id, sample_id)
    tar_cmd = 'BioSAK hpc3 -q cpu -a oces -wt 11:59:59 -t 1 -conda mybase2 -n tar%s -c "%s"'    % (file_index, tar_cmd)
    #print(tar_cmd)

    # compress assemblies
    gzip_cmd = 'gzip %s.fasta' % sample_id
    gzip_cmd = 'BioSAK hpc3 -q cpu-share -wt 11:59:59 -t 36 -conda mybase2 -n gzip%s -c "%s"' % (file_index, gzip_cmd)
    # print(gzip_cmd)

    # mkblastdb
    mkblastdb_cmd = 'makeblastdb -in %s.fasta -dbtype nucl -parse_seqids' % sample_id
    mkblastdb_cmd = 'BioSAK hpc3 -q cpu -a oces -wt 11:59:59 -t 36 -conda mybase2 -n mkblastdb_%s -c "%s"' % (file_index, mkblastdb_cmd)
    # print(mkblastdb_cmd)

    # blast barcoding gene against assemblies
    blast_cmd_28s = 'blastn -query /scratch/PI/ocessongwz/SMP/genome_skimming/blast/28S_crop.fa -db %s.fasta -out %s_28S.txt -evalue 1e-5 -outfmt 6 -task blastn' % (sample_id, sample_id)
    blast_cmd_coi = 'blastn -query /scratch/PI/ocessongwz/SMP/genome_skimming/blast/COI_crop.fa -db %s.fasta -out %s_COI.txt -evalue 1e-5 -outfmt 6 -task blastn' % (sample_id, sample_id)
    blast_cmd_both = '%s; %s' % (blast_cmd_28s, blast_cmd_coi)
    blast_cmd_both = 'BioSAK hpc3 -q cpu-share -wt 11:59:59 -t 36 -conda mybase2 -n blast_%s -c "%s"' % (file_index, blast_cmd_both)
    # print(blast_cmd_both)

    ############################################ extract_barcoding_gene_seq ############################################

    assembly_fa  = '%s/%s.fasta'   % (assembly_dir, sample_id)
    blast_op_28s = '%s/%s_28S.txt' % (blast_op_dir, sample_id)
    blast_op_coi = '%s/%s_COI.txt' % (blast_op_dir, sample_id)
    if extract_barcoding_gene_seq is True:
        if os.path.isfile(blast_op_28s):
            get_hit_seq_seq(sample_id, '28S', assembly_fa, ref_seq_file_28s, blast_op_28s, genome_skimming_28s_fa)
        if os.path.isfile(blast_op_coi):
            get_hit_seq_seq(sample_id, 'COI', assembly_fa, ref_seq_file_coi, blast_op_coi, genome_skimming_coi_fa)

    ##################################################### binning ######################################################

    # decompress
    gunzip_cmd = 'gunzip %s_1_P.fq.gz; gunzip %s_2_P.fq.gz' % (sample_id, sample_id)
    gunzip_cmd = 'BioSAK hpc3 -q cpu-share -wt 11:59:59 -t 36 -conda mybase2 -n gunzip%s -c "%s"' % (file_index, gunzip_cmd)
    # print(gunzip_cmd)

    # rename fq files for metaWRAP
    rename_cmd = 'mv %s_1_P.fq %s_1.fastq; mv %s_2_P.fq %s_2.fastq' % (sample_id, sample_id, sample_id, sample_id)
    #print(rename_cmd)

    # binning
    metawrap_cmd = 'metaWRAP binning -t 36 --metabat2 --maxbin2 --concoct -a /scratch/PI/boqianpy/songweizhi/genome_skimming_2/spades_k55-127_scaffolds/%s.fasta -o %s_metaWRAP_wd /scratch/PI/boqianpy/songweizhi/genome_skimming_2/reads_after_qc/%s_1.fastq /scratch/PI/boqianpy/songweizhi/genome_skimming_2/reads_after_qc/%s_2.fastq' % (sample_id, sample_id, sample_id, sample_id)
    metawrap_cmd = 'BioSAK hpc3 -q cpu -a boqianpy -wt 71:59:59 -t 36 -conda metawrap-env -n metawrap%s -c "%s"' % (file_index, metawrap_cmd)
    # print(metawrap_cmd)

    # refine
    refine_cmd = 'metawrap bin_refinement -o %s_refine_wd -t 36 -c 50 -x 5 -A %s_metaWRAP_wd/metabat2_bins -B %s_metaWRAP_wd/maxbin2_bins -C %s_metaWRAP_wd/concoct_bins' % (sample_id, sample_id, sample_id, sample_id)
    refine_cmd = 'BioSAK hpc3 -q cpu -a oces -wt 23:59:59 -t 36 -conda metawrap-env -n refine%s -c "%s"' % (file_index, refine_cmd)
    # print(refine_cmd)

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
