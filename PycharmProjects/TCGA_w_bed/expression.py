




























raise SystemError(0)


#from 6/1/16

with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/TCGA_splice_site_mutations.txt", "U") as mutations:


    with open("/Users/radhikarawat/PycharmProjects/TCGA_BRCA_exp_HiSeqV2_exon-2015-02-24/genomicMatrix",
              "U") as samples:
        with open("/Users/radhikarawat/PycharmProjects/TCGA_BRCA_exp_HiSeqV2_exon-2015-02-24/probe",
                  "U") as gene_matchup:


            samples.readline()
            bps=[]
            for sample_row in samples:
                sample_row = sample_row.split('\t')
                bps.append(sample_row[0])
            genes=""
            for gene_row in gene_matchup:
                gene_row=gene_row.split('\t')
                for bp_set in bps:
                    if str(bp_set) in str(gene_row[0]):
                        genes+=str(gene_row[1])+'\t'

            genes+='\n'




        with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/results.txt","w") as output:
            output.write(genes)

            for row in mutations:
                row=row.split('\t')
                print row
                if row[0] !="":
                    index=int(row[0])
                    exonID=row[10]
                    data=''
                    for sample_row in samples:
                        sample_row.split('\t')
                        data+=str(sample_row[index])+'\t'
                    print data
                        #mutation in exon-> decrease/increase of gene a, b, c, etc.

                    output.write('%s\t%s'  %(exonID, data))



