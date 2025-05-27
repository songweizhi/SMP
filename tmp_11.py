import os
import glob


coral_microbiome_project_sample_set = set()
for line in open('/Users/songweizhi/Desktop/SMP/00_Manuscript_Coral/Coral microbiome manuscript v12 - Supplementary file.txt'):
    if not line.startswith('Sample_ID\t'):
        coral_microbiome_project_sample_set.add(line.strip().split()[0])
print(len(coral_microbiome_project_sample_set))

sponge_microbiome_project_sample_set = set()
for line in open('/Users/songweizhi/Desktop/SMP/00_metadata/Sponge_samples_20250524.txt'):
    sponge_microbiome_project_sample_set.add(line.strip().split()[0])
print(len(sponge_microbiome_project_sample_set))

ignored_sponge_microbiome_project_sample_set = set()
for line in open('/Users/songweizhi/Desktop/SMP/00_metadata/Sponge_samples_20250524_ignored.txt'):
    ignored_sponge_microbiome_project_sample_set.add(line.strip().split()[0])
print(len(ignored_sponge_microbiome_project_sample_set))


metadata_txt = '/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250526_SZ_WS.txt'
metadata_new = '/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250526_SZ_WS_updated.txt'


metadata_new_handle = open(metadata_new, 'w')
for line in open(metadata_txt):
    if line.startswith('0_Sample_Index'):
        metadata_new_handle.write('%s\tCoral_Microbiome_Project\tSponge_Microbiome_Project\n' % line.strip())
    else:
        line_split = line.strip().split('\t')
        sample_id = line_split[1]

        in_cmp = 'no'
        if sample_id in coral_microbiome_project_sample_set:
            in_cmp = 'yes'

        if sample_id in sponge_microbiome_project_sample_set:
            in_smp = 'yes'
        elif sample_id in ignored_sponge_microbiome_project_sample_set:
            in_smp = 'no'
        else:
            in_smp = 'to_be_determined'
        metadata_new_handle.write('%s\t%s\t%s\n' % (line.strip(), in_cmp, in_smp))
metadata_new_handle.close()


# file_dir = '/Users/songweizhi/Desktop/SMP/01_fna_files_combined_DY86II_Haima_Hongkong_357'
# file_ext = 'fna'
# file_re = '%s/*.%s' % (file_dir, file_ext)
# file_list = [os.path.basename(i)[:-4] for i in glob.glob(file_re)]
# print(file_list)
# print(len(file_list))
#
#
# for file in file_list:
#     if file not in samples_with_metadata:
#         print(file)

