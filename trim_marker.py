import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


def best_hit(file_in, file_out):

    file_out_handle = open(file_out, 'w')
    best_hit_line = ''
    best_hit_query_id = ''
    best_hit_score = 0
    for blast_hit in open(file_in):
        blast_hit_split = blast_hit.strip().split('\t')
        query_id = blast_hit_split[0]
        bit_score = float(blast_hit_split[11])

        if best_hit_query_id == '':
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

        elif (query_id == best_hit_query_id) and (bit_score > best_hit_score):
            best_hit_score = bit_score
            best_hit_line = blast_hit

        elif query_id != best_hit_query_id:
            file_out_handle.write(best_hit_line)
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

    file_out_handle.write(best_hit_line)
    file_out_handle.close()


def trim_marker(args):

    ref_fa            = args['ref']
    to_trim_fa        = args['trim']
    op_dir            = args['o']
    force_overwrite   = args['f']

    blast_op          = '%s/blastn.txt'           % op_dir
    blast_op_best_hit = '%s/blastn_best_hit.txt'  % op_dir
    fa_out            = '%s/matched_region.fa'    % op_dir

    if os.path.isdir(op_dir):
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('%s already exists, program exited!' % op_dir)
            exit()
    os.mkdir(op_dir)

    seq_dict = dict()
    for each_seq in SeqIO.parse(to_trim_fa, 'fasta'):
        seq_dict[each_seq.id] = str(each_seq.seq)
    for each_seq in SeqIO.parse(ref_fa, 'fasta'):
        seq_dict[each_seq.id] = str(each_seq.seq)

    blast_cmd_28s = 'blastn -query %s -subject %s -out %s -evalue 1e-5 -outfmt 6 -task blastn' % (to_trim_fa, ref_fa, blast_op)
    os.system(blast_cmd_28s)

    best_hit(blast_op, blast_op_best_hit)

    fa_out_handle = open(fa_out, 'w')
    for line in open(blast_op_best_hit):
        line_split  = line.strip().split('\t')
        q_id        = line_split[0]
        s_id        = line_split[1]
        q_start     = int(line_split[6])
        q_end       = int(line_split[7])
        s_start     = int(line_split[8])
        s_end       = int(line_split[9])

        q_direction = 1
        if q_end < q_start:
            q_direction = -1

        s_direction = 1
        if s_end < s_start:
            s_direction = -1

        q_seq_matched_region = seq_dict[q_id][(q_start-1) : q_end]
        if q_direction != s_direction:
            q_seq_matched_region = Seq(q_seq_matched_region).reverse_complement()

        fa_out_handle.write('>' + q_id + '\n')
        fa_out_handle.write(str(q_seq_matched_region) + '\n')
    fa_out_handle.close()


if __name__ == '__main__':

    trim_marker_parser = argparse.ArgumentParser()
    trim_marker_parser.add_argument('-ref',  required=True,                       help='reference sequences')
    trim_marker_parser.add_argument('-trim', required=True,                       help='sequences to trim')
    trim_marker_parser.add_argument('-o',    required=True,                       help='output directory')
    trim_marker_parser.add_argument('-f',    required=False, action="store_true", help='force overwrite')
    args = vars(trim_marker_parser.parse_args())
    trim_marker(args)
