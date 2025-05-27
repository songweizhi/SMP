import argparse


def cat_blca_tax(args):

    pc  = args['pc']
    sc  = args['sc']
    out = args['o']

    seq_set = set()
    pc_dict = dict()
    for each_pc in open(pc):
        each_pc_split = each_pc.strip().split('\t')
        pc_dict[each_pc_split[0]] = each_pc_split[1]
        seq_set.add(each_pc_split[0])

    sc_dict = dict()
    for each_sc in open(sc):
        each_sc_split = each_sc.strip().split('\t')
        sc_dict[each_sc_split[0]] = each_sc_split[1]
        seq_set.add(each_sc_split[0])

    # combine classifications
    combined_dict = dict()
    unclassified_set = set()
    for each_query in seq_set:
        pc_tax = pc_dict.get(each_query, 'Unclassified')
        sc_tax = sc_dict.get(each_query, 'Unclassified')
        if pc_tax != 'Unclassified':
            combined_dict[each_query] = pc_tax
        else:
            if sc_tax != 'Unclassified':
                combined_dict[each_query] = sc_tax
            else:
                unclassified_set.add(each_query)

    # write out
    out_handle = open(out, 'w')
    for each_q in sorted(list(unclassified_set)):
        out_handle.write('%s\tUnclassified\n' % each_q)
    for each_q in sorted(list(combined_dict.keys())):
        q_tax = combined_dict[each_q]
        out_handle.write('%s\t%s\n' % (each_q, q_tax))
    out_handle.close()


if __name__ == '__main__':

    cat_blca_tax_parser = argparse.ArgumentParser()
    cat_blca_tax_parser.add_argument('-pc',    required=True,  help='primary classifications')
    cat_blca_tax_parser.add_argument('-sc',    required=True,  help='secondary classifications')
    cat_blca_tax_parser.add_argument('-o',     required=True,  help='combined classifications')
    args = vars(cat_blca_tax_parser.parse_args())
    cat_blca_tax(args)
