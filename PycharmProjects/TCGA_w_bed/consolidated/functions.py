def create_mutation_bed():
    with open("/Users/radhikarawat/PycharmProjects/CosmicMutantExport.tsv", "U") as mutation_data_file:
        with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/mutation_bed.bed", "w") as f:
            for row in mutation_data_file:
                row = row.split("\t")
                exonCount = row[8]
                if row[22] == "38":
                    if "breast" in row:
                        mutchrom, bps = row[23].split(":")
                        mutchrom = "chr" + mutchrom
                        mutstart, mutend = bps.split("-")
                        name = row[1]
                        gene_name = row[0]
                        TCGA_info = row[4]
                        # score=row[11]

                        f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (mutchrom, mutstart, mutend, name, gene_name, TCGA_info))

def create_exon_bed(bracket=1):
      #
    ##format row
    with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/transcripts.txt", "U") as exon_data_file:
        with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/exon_splice_sites_bed.bed", "w") as f:
            for row in exon_data_file:
                row = row.split("\t")

                if "e" not in row[8]:
                    exonCount = int(row[8])
                    if exonCount > 1:
                        chrom = row[2]
                        name = row[1]
                        name2 = row[12]
                        score = row[11]
                        exonStarts = row[9].split(",")
                        exonEnds = row[10].split(",")
                        numexons = len(exonEnds)
                        for i in range(0, numexons):
                            if exonStarts[i] != "":
                                if exonEnds[i] != "":
                                    SpliceSiteStart1 = str(int(exonStarts[i]) - bracket)
                                    SpliceSiteEnd1 = str(int(exonStarts[i]) + bracket)
                                    if int(SpliceSiteStart1) > 0:
                                        f.write("%s\t%s\t%s\t%s\n" % (
                                        chrom, SpliceSiteStart1, SpliceSiteEnd1, name + '_' + str(i)))
                                    SpliceSiteStart2 = str(int(exonEnds[i]) - bracket)
                                    SpliceSiteEnd2 = str(int(exonEnds[i]) + bracket)
                                    if int(SpliceSiteStart2) > 0:
                                        f.write("%s\t%s\t%s\t%s\t\n" % (chrom, SpliceSiteStart2, SpliceSiteEnd2,
                                                                        name + '_' + str(i)))
    return bracket

def intersect_bed_files(bracket=1):
    import os
    output_path="/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/bed_intersect_"+str(bracket)
    os.system("bedtools intersect -a /Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/mutation_bed.bed -b /Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/exon_splice_sites_bed.bed -wo > "+output_path)
    return output_path

def process_splice_mutations(bed_intersect_path):
    f = open("/Users/radhikarawat/PycharmProjects/TCGA_BRCA_exp_HiSeqV2_exon-2015-02-24/genomicMatrix")
    samples = f.readline()
    f.close()
    length = len(samples)

    with open(bed_intersect_path, "U") as f:

        with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/TCGA_splice_site_mutations.txt", "w") as output:
            for row in f:
                if "TCGA" in row:
                    splitrow = row.split("\t")
                    if splitrow[5] in samples:
                        output.write("%i\t%s\n" % (int(samples.index(splitrow[5])), splitrow[-2]))
    print "splice mutations list created: index of sample in genomix matrix \t exon_name"


def create_patient_sample_header():
    with open("/Users/radhikarawat/PycharmProjects/TCGA_BRCA_exp_HiSeqV2_exon-2015-02-24/genomicMatrix","U") as samples:
        with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/patient_sample#_header.txt","w") as patient_sample_header:
            patient_sample_header.write(samples.readline())
    print "patient_sample_header created"

def create_gene_dict():
    genes_dict = {}  # contains gene coordinates (chr:start-end:strand) and gene name ("name")
    with open("/Users/radhikarawat/PycharmProjects/TCGA_BRCA_exp_HiSeqV2_exon-2015-02-24/probe", "U") as gene_probe:
        for gene_row in gene_probe:
            gene_row = gene_row.split('\t')
            genes_dict[gene_row[0]] = gene_row[1]
    return genes_dict

def add_gene_names_to_matrix(gene_dictionary):
    with open("/Users/radhikarawat/PycharmProjects/TCGA_BRCA_exp_HiSeqV2_exon-2015-02-24/genomicMatrix",
              "U") as genomic_matrix:
        with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/TCGA_splice_site_mutations.txt") as mutations:
            with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/results.txt") as results:
                for mutation in mutations:
                    mutation = mutation.split('\t')
                    index, exon_name = mutation[0], mutation[10]
                    for row in genomic_matrix:
                        row = row.split('\t')
                        if row[0] in gene_dictionary.keys:
                            row[0] = gene_dictionary[row[0]]
                        results.write(exon_name + '\t' +'')



