
s_set = set()
for each in open('/Users/songweizhi/Desktop/13.txt'):
    s_set.add(each.strip())


for each in open('/Users/songweizhi/Desktop/SMP/00_metadata/metadata_20250905.txt'):
    each_split = each.strip().split('\t')
    if each_split[1] in s_set:
        print(each.strip())

