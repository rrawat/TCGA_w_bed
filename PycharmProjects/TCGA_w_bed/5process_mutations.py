#119	chr11	100193393	100193529	ENST00000528682	CNTN5	TCGA-A8-A08L-01	chr11	100193484	100193490
# ENST00000527185.2_14	6


with open("/Users/radhikarawat/PycharmProjects/TCGA_BRCA_exp_HiSeqV2_exon-2015-02-24/genomicMatrix","U") as samples:
    print "genomic matrix rows:" , sum(1 for line in samples),
    for line in samples.readline():print len(line)

with open("/Users/radhikarawat/PycharmProjects/TCGA_BRCA_exp_HiSeqV2_exon-2015-02-24/probe") as probe:
    print "probe rows", sum(1 for row in probe)

for line in open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/4out_gene_base_pairs.txt"):
    print "gene bp entries in list:", sum (1 for element in line)

for line in open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/4out_genes.txt"):
    print "gene entries in list:", sum (1 for element in line)

raise SystemError(0)

with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/TCGA_splice_site_mutations.txt",
          "U") as TCGA_splice_site_mutations:
    for TCGA_row in TCGA_splice_site_mutations:
        row=row.split("\t")
        index=row[0]



with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/5out_exon_list.txt", 'w') as exon_list:
    exon=TCGA_row[10]
    exon_list.write(exon)


with open ("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/5out_results.txt",'w') as results:
    results.write(data)

