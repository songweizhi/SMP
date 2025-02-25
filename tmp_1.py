
def rarefaction(otu_table_txt, sample_group_txt, group_color_txt, output_plot):

    pwd_current_file    = os.path.realpath(__file__)
    current_file_path   = '/'.join(pwd_current_file.split('/')[:-1])
    rarefaction_R       = '%s/rarefaction_R' % current_file_path
    if os.path.isfile(rarefaction_R) is False:
        print('rarefaction.R not found, program exited!')
        exit()

    # define file name
    grouping_txt    = '%s/grouping.txt'     % op_dir
    group_color_txt = '%s/color.txt'        % op_dir
    output_plot     = '%s/rarefaction.pdf'  % op_dir

    sample_list = open(otu_table_txt).readline().strip().split('\t')[1:]

    metadata_dict = dict()
    for sample in open(metadata_txt):
        sample_split = sample.strip().split('\t')
        sample_id = sample_split[0]
        sample_grp = sample_split[1]
        metadata_dict[sample_id] = sample_grp

    sample_group_txt_handle = open(grouping_txt, 'w')
    sample_group_txt_handle.write('SampleID\tSampleGroup\n')
    group_set = set()
    for sample in sorted(sample_list):
        sample_grp = metadata_dict.get(sample, 'na')
        sample_group_txt_handle.write('%s\t%s\n' % (sample, sample_grp))
        group_set.add(sample_grp)
    sample_group_txt_handle.close()

    color_code_dict = dict()
    for each_sample_type in open(color_code_txt):
        sample_type_split = each_sample_type.strip().split('\t')
        sample_type = sample_type_split[0]
        sample_color = sample_type_split[1]
        color_code_dict[sample_type] = sample_color

    group_color_txt_handle = open(group_color_txt, 'w')
    group_color_txt_handle.write('GroupID\tGroupColor\n')
    for grp in sorted(list(group_set)):
        grp_color = color_code_dict[grp]
        group_color_txt_handle.write('%s\t%s\n' % (grp, grp_color))
    group_color_txt_handle.close()

    rarefaction_cmd = 'Rscript %s -i %s -g %s -c %s -o %s' % \
    (rarefaction_R, otu_table_txt, grouping_txt, group_color_txt, output_plot)
    print(rarefaction_cmd)
    os.system(rarefaction_cmd)

    print('Done')


