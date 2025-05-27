import os
import glob
import numpy as np
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def plot_seq_len(fna_file, png_file):

    len_list = []
    len_set = set()
    for seq in SeqIO.parse(fna_file, 'fasta'):
        len_list.append(len(seq.seq))
        len_set.add(len(seq.seq))
    plt.hist(len_list)
    plt.savefig(png_file, dpi=300)
    plt.close()



file_dir = '/Users/songweizhi/Desktop/SMP/01_fna_files_combined_DY86II_Haima_Hongkong_357'
file_ext = 'fna'
file_re = '%s/*.%s' % (file_dir, file_ext)
file_list = [os.path.basename(i)[:-4] for i in glob.glob(file_re)]


for file in file_list:
    pwd_fna = '/Users/songweizhi/Desktop/SMP/01_fna_files_combined_DY86II_Haima_Hongkong_357/%s.fna' % (file)
    pwd_png = '/Users/songweizhi/Desktop/plots/%s.png' % (file)
    plot_seq_len(pwd_fna, pwd_png)
