from datetime import datetime
import numpy

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
    output_path = "/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/bed_intersect_" + str(bracket)
    os.system(
        "bedtools intersect -a /Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/mutation_bed.bed -b /Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/exon_splice_sites_bed.bed -wo > " + output_path)
    return output_path


def process_splice_mutations(bed_intersect_path):
    f = open("/Users/radhikarawat/PycharmProjects/TCGA_BRCA_exp_HiSeqV2_exon-2015-02-24/genomicMatrix")
    samples = f.readline()
    f.close()
    split_samples = samples.split('\t')

    with open(bed_intersect_path, "U") as f:

        with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/TCGA_splice_site_mutations.txt",
                  "w") as output:

            for row in f:

                if "TCGA" in row:
                    splitrow = row.split("\t")
                    if splitrow[5] in split_samples:
                        index = int(split_samples.index(splitrow[5]))
                        print len(split_samples), index,
                        output.write("%i\t%s\n" % (index, splitrow[-2]))
    print "splice mutations list created: index of sample in genomix matrix \t exon_name"


def create_patient_sample_header():
    with open("/Users/radhikarawat/PycharmProjects/TCGA_BRCA_exp_HiSeqV2_exon-2015-02-24/genomicMatrix",
              "U") as samples:
        with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/patient_sample#_header.txt",
                  "w") as patient_sample_header:
            patient_sample_header.write(samples.readline())
    print "patient_sample_header created"


def create_gene_dict():
    genes_dict = {}  # contains gene coordinates (chr:start-end:strand) and gene name ("name")
    with open("/Users/radhikarawat/PycharmProjects/TCGA_BRCA_exp_HiSeqV2_exon-2015-02-24/probe", "U") as gene_probe:
        for gene_row in gene_probe:
            gene_row = gene_row.split('\t')
            genes_dict[gene_row[0]] = gene_row[1]
    return genes_dict


# def exon_dict():
#     exon_dict={} #creates dict from splice mutations file
#     with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/TCGA_splice_site_mutations.txt") as \
#             mutations:
#         for row in mutations:
#     return exon_dict

def trim_matrix():
    start = datetime.now()
    # get matrix columns only for indexes in mutations list:
    dict = {}
    with open(
            "/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/TCGA_splice_site_mutations.txt") as mutations:
        for row in mutations:
            row = row.split('\t')
            dict[row[0]] = ''

    with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/SplitMatrix.txt",
              "w") as output_matrix:
        with open("/Users/radhikarawat/PycharmProjects/TCGA_BRCA_exp_HiSeqV2_exon-2015-02-24/genomicMatrix",
                  "U") as genomic_matrix:
            j = 0
            for row in genomic_matrix:
                if j % 10000 == 0:
                    print j, start, datetime.now() - start
                j += 1
                row = row.split('\t')
                output_matrix.write(row[0])

                for i in dict.keys():
                    # print len(row), int(i), len(row)>i
                    output_matrix.write('\t' + row[int(i)])
                output_matrix.write('\n')
    print 'elapsed time:', datetime.now() - start


def exons_to_matrix():
    with open(
            "/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/TCGA_splice_site_mutations.txt") as mutations:
        with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/SplitMatrix.txt", "U") as genomic_matrix:
            genomic_matrix.readline()

            with open(
                    "/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/exon_matrix.txt", "w") as exon_matrix:

                with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/gene_bp_list.txt",
                          'w') as bp_list:
                    j = 0
                    for row in mutations:
                        print j

                        index, exon_name = row.split('\t')
                        exon_data = ''
                        gene_bp_list = ''
                        i = 0
                        for gene_row in genomic_matrix:
                            gene_row1 = gene_row.split('\t')
                            if i == 0:
                                i += 1
                                # bp_list.write(gene_row[0]+'\t')

                            addition = gene_row1[int(index)] + '\t'
                            exon_data += addition
                            if j > 0: print addition
                        j += 1
                        print exon_data
                        exon_matrix.write(exon_name.split('\n')[0] + '\t' + exon_data + '\n')

# def add_gene_names_to_matrix(gene_dictionary):
#     with open("/Users/radhikarawat/PycharmProjects/TCGA_BRCA_exp_HiSeqV2_exon-2015-02-24/genomicMatrix",
#               "U") as genomic_matrix:
#         with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/TCGA_splice_site_mutations.txt") as mutations:
#             with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/results.txt") as results:
#                 for mutation in mutations:
#                     mutation = mutation.split('\t')
#                     index, exon_name = mutation[0], mutation[10]
#                     for row in genomic_matrix:
#                         row = row.split('\t')
#                         if row[0] in gene_dictionary.keys:
#                             row[0] = gene_dictionary[row[0]]
#                         results.write(exon_name + '\t' +'')


def exons_to_matrix_using_np():

    start=datetime.now()
    with open('/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/TCGA_splice_site_mutations.txt','U') as mutation_file:
        mutation_array= numpy.loadtxt(mutation_file,dtype=str)

    print mutation_array, mutation_array.shape, datetime.now()-start



    with open("/Users/radhikarawat/PycharmProjects/TCGA_BRCA_exp_HiSeqV2_exon-2015-02-24/genomicMatrix", "U") as matrix_file:

        split_matrix_array=numpy.loadtxt(matrix_file,dtype=str)
    print split_matrix_array, split_matrix_array.shape,datetime.now()-start

    new_array=split_matrix_array

    for j in range (0,17796-1): #every index in mutation list
        index,exon_name=mutation_array[j,:]
        new_array[0,int(index)]=exon_name


    print new_array, new_array.shape, datetime.now()-start

    delete_list=[]
    for i in range(0,new_array.shape[1]):
        if "TCGA" in str(new_array[0,i]):
            delete_list.append(i)
    print delete_list, datetime.now()-start
    output_array=numpy.delete(new_array,delete_list,1)

    with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/results_array2.txt", "w+") as output_file:
        numpy.savetxt(output_file, output_array, delimiter="\t", newline="\n",fmt="%s")

    print datetime.now()-start

def add_genes_to_matrix_np():
    start=datetime.now()
    with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/results_array2.txt", "r+") as array_file:
        array=numpy.loadtxt(array_file, dtype=str)

    gene_dict=create_gene_dict()
    for i in range(0,array.shape[0]):
        gene_bps = str(array[i, 0])
        #print i
        try:
            #if gene_bps in gene_dict.keys():
            array[i,0]=gene_dict[gene_bps]
            gene_dict[gene_bps]=str(gene_dict[gene_bps])+'*'
        except KeyError: array[i,0]=gene_bps
    with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/exon_v_genes_results_unique.txt",
              "w") as output_file:
        numpy.savetxt(output_file, array, delimiter="\t", newline="\n", fmt="%s")
    print datetime.now()-start, start

def create_plot():
    start = datetime.now()
    import pandas as pd
    df = pd.read_csv("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/exon_v_genes_results.txt", sep='\t')

    info = numpy.array(df)
    # print info[:5,:4]

    print datetime.now()-start, "file loaded"
    # raise SystemError
    # #
    #
    #start=datetime.now()
    import matplotlib.pyplot as plt
    # with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/consolidated/exon_v_genes_results.txt",
    #           "r+") as array_file:
    #     print "A"
    #     info=numpy.genfromtxt(array_file, dtype='str', filling_values=0.0)
    #     # info = numpy.loadtxt(, dtype=str)
    #print "loaded", datetime.now()-start
    # #for gene 0
    y = info[1,1:]

    names = info[0,1:]
    width = 1 / 1.5
    print names, y

    import matplotlib.pyplot as plt
    import numpy as np


    x = np.arange(len(y))

    plt.bar(x, y)
    plt.xticks(x + 0.5, names, rotation=90)
    try: plt.show()
    except: plt.savefig()
