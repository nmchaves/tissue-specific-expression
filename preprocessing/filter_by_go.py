"""
    This file performs the 1st step of preprocessing of the GTEx RPKM file:
    "GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt".

    It goes through each line of the above rpkm file. If the current transcript
    corresponds to a gene in Gene Ontology, then the current line is written out
    to the file "transcript_rpkm_in_go.txt". Otherwise, the current line is written
    out to the file "transcript_rpkm_not_in_go.txt".

"""

"""
*********************
        Main
*********************
"""
if __name__ == "__main__":

    path_to_rpkm_file = '../../../Downloads/GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt'
    path_to_go_file = '../../../Downloads/trans_gene_name_filtered_by_GO.txt'

    # Read in GO file and generate dictionary containing transcripts.
    go_file = open(path_to_go_file)
    transcriptIdsInGO = {}
    firstLine = True
    for line in go_file:
        if firstLine:
            firstLine = False
            continue
        transcriptId = line[0:line.index('\t')]
        transcriptIdsInGO[transcriptId] = True
    go_file.close()

    rpkm_file = open(path_to_rpkm_file)
    rpkm_file_retained = open('transcript_rpkm_in_go.txt', 'w')
    rpkm_file_removed = open('transcript_rpkm_not_in_go.txt', 'w')
    firstLine = True
    for line in rpkm_file:
        if firstLine:
            firstLine = False
            # Write the header row into both files
            rpkm_file_retained.write(line)
            rpkm_file_removed.write(line)
            continue
        else:
            transcriptId = line[0:line.index('\t')]
            if '.' in transcriptId:
                transcriptId = transcriptId[0:transcriptId.index('.')]

            if transcriptId in transcriptIdsInGO:
                rpkm_file_retained.write(line)
            else:
                rpkm_file_removed.write(line)

    rpkm_file_retained.close()
    rpkm_file_removed.close()


