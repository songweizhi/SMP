import os
from Bio import SeqIO


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def run_mmseqs_linclust(pwd_combined_faa, num_threads, iden_cutoff, cov_cutoff, force_overwrite, output_folder):

    # create output folder
    if os.path.isdir(output_folder) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % output_folder)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % output_folder)

    # define fine name
    faa_name, faa_path, faa_base, faa_ext = sep_path_basename_ext(pwd_combined_faa)
    mmseqs_db   = '%s.db'                                           % pwd_combined_faa
    mmseqs_clu  = '%s/%s.mmseqs.iden%s.cov%s'                       % (output_folder, faa_base, iden_cutoff, cov_cutoff)
    mmseqs_tmp  = '%s/%s.mmseqs.iden%s.cov%s.tmp'                   % (output_folder, faa_base, iden_cutoff, cov_cutoff)
    mmseqs_tsv  = '%s/%s.mmseqs.iden%s.cov%s.tsv'                   % (output_folder, faa_base, iden_cutoff, cov_cutoff)
    rep_seq_txt = '%s/%s.mmseqs.iden%s.cov%s.representatives.txt'   % (output_folder, faa_base, iden_cutoff, cov_cutoff)
    rep_seq_faa = '%s/%s.mmseqs.iden%s.cov%s.representatives.faa'   % (output_folder, faa_base, iden_cutoff, cov_cutoff)

    mmseqs_createdb_cmd  = 'mmseqs createdb %s %s > /dev/null' % (pwd_combined_faa, mmseqs_db)
    # print(mmseqs_createdb_cmd)
    os.system(mmseqs_createdb_cmd)

    mmseqs_cluster_cmd = 'mmseqs linclust %s %s %s --threads %s --min-seq-id %s --cov-mode 1 -c %s > /dev/null' % (mmseqs_db, mmseqs_clu, mmseqs_tmp, num_threads, iden_cutoff, cov_cutoff)
    # print(mmseqs_cluster_cmd)
    os.system(mmseqs_cluster_cmd)

    mmseqs_createtsv_cmd = 'mmseqs createtsv %s %s %s %s > /dev/null' % (mmseqs_db, mmseqs_db, mmseqs_clu, mmseqs_tsv)
    # print(mmseqs_createtsv_cmd)
    os.system(mmseqs_createtsv_cmd)

    # get representative sequences
    rep_seq_set = set()
    for each_line in open(mmseqs_tsv):
        each_line_split = each_line.strip().split('\t')
        rep_seq_set.add(each_line_split[0])

    with open(rep_seq_txt, 'w') as f:
        f.write('\n'.join(sorted(list(rep_seq_set))))

    rep_seq_faa_handle = open(rep_seq_faa, 'w')
    for each_seq in SeqIO.parse(pwd_combined_faa, 'fasta'):
        if each_seq.id in rep_seq_set:
            rep_seq_faa_handle.write('>' + each_seq.id + '\n')
            rep_seq_faa_handle.write(str(each_seq.seq) + '\n')
    rep_seq_faa_handle.close()

    print('cov%s iden%s: %s' % (cov_cutoff, iden_cutoff, len(rep_seq_set)))


pwd_combined_faa    = '/Users/songweizhi/Desktop/Jianwei_Maggie/All_PeptidaseS9-RiPPs_dna100_realRiPPS9100_addref.faa'
num_threads         = '10'
force_overwrite     = True

cov_cutoff_list     = [0.75, 0.8, 0.85]
iden_cutoff_list    = [0.9, 0.95, 0.97, 0.99]

cov_cutoff_list     = [0.85]
iden_cutoff_list    = [0.35, 0.5, 0.6, 0.7, 0.8]

for cov_cutoff in cov_cutoff_list:
    for iden_cutoff in iden_cutoff_list:
        op_dir              = '/Users/songweizhi/Desktop/Jianwei_Maggie/mmseqs_cov%s_iden%s' % (cov_cutoff, iden_cutoff)
        run_mmseqs_linclust(pwd_combined_faa, num_threads, iden_cutoff, cov_cutoff, force_overwrite, op_dir)


'''

cov0.85 iden0.35: 667
cov0.85 iden0.5: 904
cov0.85 iden0.6: 1121
cov0.85 iden0.7: 1341
cov0.85 iden0.8: 1624

'''
