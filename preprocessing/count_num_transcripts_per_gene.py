
"""
*********************
        Main
*********************
"""
if __name__ == "__main__":

    rpkm_filtered_path = '../../../../Documents/Stanford/CS341_Data/transcript_rpkm_top_10000_var.txt'
    rpkm_filtered = open(rpkm_filtered_path)

    gene_counts = {}
    firstLine = True
    for line in rpkm_filtered:
        if firstLine:
            firstLine = False
            continue
        first_tab_index = line.index('\t')
        second_tab_index = line.index('\t', first_tab_index+1)
        gene_id = line[first_tab_index+1:second_tab_index]
        gene_id = gene_id[0:gene_id.index('.')]
        if gene_id in gene_counts:
            gene_counts[gene_id] += 1
        else:
            gene_counts[gene_id] = 1

    gene_counts_file = open('../data/num_transcripts_per_gene.txt', 'w')
    gene_counts_file.write('Gene_Symbol\tNum_TargetIDs\n')

    for gene_id in gene_counts.keys():
        gene_counts_file.write(gene_id + '\t' + str(gene_counts[gene_id]) + '\n')