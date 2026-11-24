
txt = '/Users/songweizhi/Desktop/Book1.txt'

sample_run_dict = dict()
for line in open(txt, 'r'):
    line_split = line.strip().split('\t')
    sample = line_split[1]
    run = line_split[0]
    if sample not in sample_run_dict:
        sample_run_dict[sample] = []
    sample_run_dict[sample].append(run)




m = 0
for each in sample_run_dict:
    if len(sample_run_dict[each]) > 1:
        print(each, sample_run_dict[each][0])
        print(each, sample_run_dict[each][1])

    elif len(sample_run_dict[each]) == 1:
        m += 1

print(m)

'''

JL302_B05_3 X101SC23061594-Z05-J020	X101SC25013054-Z01-J015
JL313_B35_4 X101SC23061594-Z05-J020	X101SC25013054-Z01-J006
JL313_B24_4 X101SC23061594-Z05-J020	X101SC25013054-Z01-J006
JL306_B12_1 X101SC23061594-Z05-J020	X101SC25013054-Z01-J006
JL305_B19_1 X101SC23061594-Z05-J020	X101SC25013054-Z01-J006
JL306_B05_2 X101SC23061594-Z05-J020	X101SC25013054-Z01-J006

'''