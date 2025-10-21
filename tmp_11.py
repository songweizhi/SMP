
sample_set = set()
id_txt = '/Users/songweizhi/Desktop/555_uniq.txt'
for each in open(id_txt):
    sample_set.add(each.strip())
print(sample_set)
print(len(sample_set))

handle = open('/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250905_updated.txt', 'w')

found_set = set()
for each in open('/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250905.txt'):
    sample_id = each.strip().split()[1]
    if sample_id in sample_set:
        handle.write('%s\tyes\n' % (each.strip()))
        found_set.add(sample_id)
    else:
        handle.write('%s\tno\n' % (each.strip()))
handle.close()


print(len(found_set))


for each in sample_set:
    if each not in found_set:
        print(each)

