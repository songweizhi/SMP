import os
import argparse


get_unclassified_ASV_domain_usage = '''
================================= get_unclassified_ASV_domain example commands =================================

python3 /Users/songweizhi/PycharmProjects/SMP/get_unclassified_ASV_domain.py -o /Users/songweizhi/Desktop/aaa.txt -n 5 -i /Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB_20250513_357/s06_AllSamples_unoise_nc_BLCA_GTDB_r226_Unclassified_3818_vs_core_nt.txt -node /Users/songweizhi/DB/taxdump_20250321/nodes.dmp -db /Users/songweizhi/DB/NCBI/RefSeq_taxonomy.txt -id /Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB_20250513_357/s06_AllSamples_unoise_nc_BLCA_GTDB_r226_Unclassified_1750.txt


python3 /Users/songweizhi/PycharmProjects/SMP/get_unclassified_ASV_domain.py -o /Users/songweizhi/Desktop/aaa.txt -n 5 -i /Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB_20250513_357/s06_AllSamples_unoise_nc_BLCA_GTDB_r226_Unclassified_3818_vs_core_nt.txt -node /Users/songweizhi/DB/taxdump_20250321/nodes.dmp -db /Users/songweizhi/DB/NCBI/RefSeq_taxonomy.txt -id /Users/songweizhi/Desktop/SMP/02_Usearch16S_20250525_356/
s08_AllSamples_unoise_nc_Unclassified_1539.txt
/Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB_20250513_357/s06_AllSamples_unoise_nc_BLCA_GTDB_r226_Unclassified_1750.txt

===============================================================================================
'''


def group_tax_by_domain(nodes_dmp):
    child_to_parent_dict = dict()
    for each_node in open(nodes_dmp):
        each_node_split = each_node.strip().split('|')
        tax_id = int(each_node_split[0].strip())
        parent_tax_id = int(each_node_split[1].strip())
        child_to_parent_dict[tax_id] = parent_tax_id

    tax_id_set_bac = set()
    tax_id_set_ar = set()
    tax_id_set_eu = set()
    tax_id_set_virus = set()
    for tax_id in child_to_parent_dict.keys():
        # get tax_lineage_list
        tax_lineage_list = [tax_id]
        while tax_id != 1:
            tax_id = child_to_parent_dict[tax_id]
            tax_lineage_list.append(tax_id)
        # Bacteria
        if 2 in tax_lineage_list:
            for each_id in tax_lineage_list[:-3]:
                tax_id_set_bac.add(each_id)
        # Archaea
        elif 2157 in tax_lineage_list:
            for each_id in tax_lineage_list[:-3]:
                tax_id_set_ar.add(each_id)
        # Eukaryota
        elif 2759 in tax_lineage_list:
            for each_id in tax_lineage_list[:-3]:
                tax_id_set_eu.add(each_id)
        # Virus
        elif 10239 in tax_lineage_list:
            for each_id in tax_lineage_list[:-3]:
                tax_id_set_virus.add(each_id)
    return tax_id_set_bac, tax_id_set_ar, tax_id_set_eu, tax_id_set_virus


def get_unclassified_ASV_domain(args):

    blast_op_txt            = args['i']
    interested_asv_txt      = args['id']
    asv_classification_txt  = args['o']
    nodes_dmp               = args['node']
    refseq_tax_txt          = args['db']
    hit_num_to_keep         = args['n']

    tax_set_bac, tax_set_ar, tax_set_eu, tax_set_virus = group_tax_by_domain(nodes_dmp)

    interested_asv_set = set()
    if os.path.isfile(interested_asv_txt):
        with open(interested_asv_txt) as f:
            for each_line in f:
                interested_asv_set.add(each_line.strip().split()[0])

    refseq_tax_dict = dict()
    for each_ref in open(refseq_tax_txt):
        each_ref_split = each_ref.strip().split('\t')
        if len(each_ref_split) == 2:
            refseq_tax_dict[each_ref_split[0]] = each_ref_split[1]

    query_to_hit_dict = dict()
    for each_line in open(blast_op_txt):
        each_line_split = each_line.strip().split('\t')
        query_id = each_line_split[0]
        subject_id = each_line_split[1]

        to_process = False
        if len(interested_asv_set) == 0:
            to_process = True
        else:
            if query_id in interested_asv_set:
                to_process = True

        # add to query_to_hit_dict
        if to_process is True:
            if query_id not in query_to_hit_dict:
                query_to_hit_dict[query_id] = set()
            if len(query_to_hit_dict[query_id]) < hit_num_to_keep:
                query_to_hit_dict[query_id].add(subject_id)

    unclassified_refseq_set = set()
    eu_otu_set = set()
    bac_otu_set = set()
    ar_otu_set = set()
    unclassified_otu_set = set()
    for query in query_to_hit_dict:

        hit_tax_list = []
        for each_hit in query_to_hit_dict[query]:
            hit_tax = refseq_tax_dict.get(each_hit, 'na')
            if hit_tax == 'na':
                hit_tax_list.append(hit_tax)
                unclassified_refseq_set.add(each_hit)
            else:
                hit_tax_list.append(int(hit_tax))

        # get stat
        hit_num_bac = 0
        hit_num_ar = 0
        hit_num_eu = 0
        hit_num_virus = 0
        hit_num_na = 0
        for hit_tax in hit_tax_list:
            if hit_tax in tax_set_bac:
                hit_num_bac += 1
            elif hit_tax in tax_set_ar:
                hit_num_ar += 1
            elif hit_tax in tax_set_eu:
                hit_num_eu += 1
            elif hit_tax in tax_set_virus:
                hit_num_virus += 1
            elif hit_tax == 'na':
                hit_num_na += 1

        # get domain level assignment
        if hit_num_bac == hit_num_to_keep:
            bac_otu_set.add(query)
        elif hit_num_ar == hit_num_to_keep:
            ar_otu_set.add(query)
        else:
            unclassified_otu_set.add(query)

    # decide ASV set to report
    asv_set_to_report = query_to_hit_dict.keys
    if len(interested_asv_set) > 0:
        asv_set_to_report = interested_asv_set

    # write out
    asv_classification_txt_handle = open(asv_classification_txt, 'w')
    for each_asv in sorted(list(asv_set_to_report)):
        if each_asv in ar_otu_set:
            asv_classification_txt_handle.write('%s\td__Archaea;p__;c__;o__;f__;g__;s__\n' % each_asv)
        elif each_asv in bac_otu_set:
            asv_classification_txt_handle.write('%s\td__Bacteria;p__;c__;o__;f__;g__;s__\n' % each_asv)
        else:
            asv_classification_txt_handle.write('%s\tUnclassified\n' % each_asv)
    asv_classification_txt_handle.close()


if __name__ == '__main__':

    get_unclassified_ASV_domain_parser = argparse.ArgumentParser(usage=get_unclassified_ASV_domain_usage)
    get_unclassified_ASV_domain_parser.add_argument('-i',    required=True,                              help='blast result, outfmt need to be 6')
    get_unclassified_ASV_domain_parser.add_argument('-id',   required=False,                             help='ASV file')
    get_unclassified_ASV_domain_parser.add_argument('-o',    required=True,                              help='output of identified eukaryotic OTUs')
    get_unclassified_ASV_domain_parser.add_argument('-node', required=True,                              help='nodes.dmp ')
    get_unclassified_ASV_domain_parser.add_argument('-db',   required=True,                              help='RefSeq taxonomy file')
    get_unclassified_ASV_domain_parser.add_argument('-n',    required=False, type=int, default=10,       help='number of top hits to consider, default is 10')
    args = vars(get_unclassified_ASV_domain_parser.parse_args())
    get_unclassified_ASV_domain(args)
