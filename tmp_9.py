
file_index = 1
for each_sample in open('/Users/songweizhi/Desktop/Sponge_r226/07_Sponge_tree/JL_samples_17.txt'):
    sample_id      = each_sample.strip()
    blast_cmd_18s  = 'blastn -query RefSeqs_with_AOA_18S.fasta -subject JL_sample_17_spades_assembly/%s.fasta -out %s_18S.txt -evalue 1e-5 -outfmt 6 -task blastn'        % (sample_id, sample_id)
    blast_cmd_coi  = 'blastn -query RefSeqs_with_AOA_COI_iden95.ffn -subject JL_sample_17_spades_assembly/%s.fasta -out %s_COI.txt -evalue 1e-5 -outfmt 6 -task blastn'   % (sample_id, sample_id)
    blast_cmd_both = '%s; %s' % (blast_cmd_18s, blast_cmd_coi)
    blast_cmd_both = 'BioSAK hpc4 -q amd -a marmolecol -wt 00:59:59 -t 36 -conda mybase2 -n blast_%s -c "%s"' % (file_index, blast_cmd_both)
    #print(blast_cmd_both)
    file_index += 1

    # get best hit
    blast_op_18s            = '%s_18S.txt'          % sample_id
    blast_op_coi            = '%s_COI.txt'          % sample_id
    blast_op_18s_best_hit   = '%s_18S_best_hit.txt' % sample_id
    blast_op_coi_best_hit   = '%s_COI_best_hit.txt' % sample_id
    best_hit_cmd_18s = 'BioSAK BestHit -i %s -o %s' % (blast_op_18s, blast_op_18s_best_hit)
    best_hit_cmd_coi = 'BioSAK BestHit -i %s -o %s' % (blast_op_coi, blast_op_coi_best_hit)
    print(best_hit_cmd_18s)
    print(best_hit_cmd_coi)
