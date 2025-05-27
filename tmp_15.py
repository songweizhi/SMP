
for each in open('/Users/songweizhi/Desktop/file.txt'):
    bin_id = each.strip()
    blast_cmd = 'blastn -query /Users/songweizhi/Desktop/SMP/genome_skimming/refined_MAGs/%s.fa -subject Zotu24_51_64_455_458_464.fasta -out %s.txt -outfmt 6' % (bin_id, bin_id)
    print(blast_cmd)




