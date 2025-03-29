
blast_op            = '/Users/songweizhi/Desktop/s06_AllSamples_unoise_nc_vs_nt.txt'
blast_op_formatted  = '/Users/songweizhi/Desktop/s06_AllSamples_unoise_nc_vs_nt_formatted.txt'
line_num_to_keep    = 99999

blast_op            = 's06_AllSamples_unoise_nc_vs_nt.txt'
blast_op_formatted  = 's06_AllSamples_unoise_nc_vs_nt_formatted.txt'
line_num_to_keep    = 99999

blast_op_formatted_handle = open(blast_op_formatted, 'w')
wrote_line_num = 0
current_query = ''
count_line = 0
for each_line in open(blast_op):
    if each_line.startswith('Query= '):
        current_query = each_line.strip()[6:].strip()
        wrote_line_num = 0
    elif each_line.startswith('Sequences producing significant alignments:'):
        count_line = 1
    elif each_line.startswith('>'):
        count_line = 0
    if (count_line == 1) and ('Sequences producing significant alignments:' not in each_line) and (len(each_line.strip()) > 0) and (wrote_line_num < line_num_to_keep):
        blast_op_formatted_handle.write('%s\t%s\n' % (current_query, each_line.strip()))
        wrote_line_num += 1
blast_op_formatted_handle.close()


'''

Zotu13917 KX624041.1 Uncultured bacterium clone 16S(V3-V4)-1310 16S ribosom...  390     5e-104

'''
