f=open("/Users/radhikarawat/PycharmProjects/TCGA_BRCA_exp_HiSeqV2_exon-2015-02-24/genomicMatrix")
samples= f.readline()
f.close()
length=len(samples)

with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/bed_intersect_b1","U") as f:

    with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/TCGA_splice_site_mutations.txt","w") as output:
        for row in f:
            if "TCGA" in row:
                splitrow=row.split("\t")
                if splitrow[5] in samples:
                    output.write("%i\t%s\n" %(int(samples.index(splitrow[5])),row))



