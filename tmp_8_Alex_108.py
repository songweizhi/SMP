import os

sample_id_txt   = '/Users/songweizhi/Desktop/SMP/01_Alex_108/PRJNA930637_108.txt'
fastq_dir       = '/Users/songweizhi/Desktop/SMP/01_Alex_108/PRJNA930637_108'


for each_sample in open(sample_id_txt):
    sample_id = each_sample.split()[0]
    fq_r1 = '%s/%s_1.fastq' % (fastq_dir, sample_id)
    fq_r2 = '%s/%s_2.fastq' % (fastq_dir, sample_id)

    # trim
    trim_cmd = 'trimmomatic PE -threads 10 -phred33 %s_1.fastq %s_2.fastq %s_1_P.fastq %s_1_UP.fastq %s_2_P.fastq %s_2_UP.fastq ILLUMINACLIP:/Users/songweizhi/DB/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 CROP:145 HEADCROP:5 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:100' % (sample_id, sample_id, sample_id, sample_id, sample_id, sample_id)
    trim_cmd = 'trimmomatic PE -threads 10 -phred33 %s_1.fastq %s_2.fastq %s_1_P.fastq %s_1_UP.fastq %s_2_P.fastq %s_2_UP.fastq ILLUMINACLIP:/Users/songweizhi/DB/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:100' % (sample_id, sample_id, sample_id, sample_id, sample_id, sample_id)
    # print(trim_cmd)
    # os.system(trim_cmd)

    # fastqc
    fastqc_cmd_raw = 'fastqc %s_1.fastq %s_2.fastq'     % (sample_id, sample_id)
    fastqc_cmd_qc  = 'fastqc %s_1_P.fastq %s_2_P.fastq' % (sample_id, sample_id)
    # print(fastqc_cmd_raw)
    # print(fastqc_cmd_qc)
    # os.system(fastqc_cmd_raw)
    # os.system(fastqc_cmd_qc)

    # merge
    fastq_mergepairs_cmd = 'usearch -fastq_mergepairs %s_1.fastq -reverse %s_2.fastq -fastqout %s_merged.fastq -fastq_minlen 250 -fastq_minovlen 8' % (sample_id, sample_id, sample_id)
    fastq_mergepairs_cmd = 'usearch -fastq_mergepairs %s_1.fastq -reverse %s_2.fastq -fastqout %s_merged.fastq -relabel @ -fastq_minovlen 8 -fastq_maxdiffs 5 -fastq_pctid 80 -fastq_minmergelen 245 -fastq_maxmergelen 275 -sample %s' % (sample_id, sample_id, sample_id, sample_id)
    print(fastq_mergepairs_cmd)


