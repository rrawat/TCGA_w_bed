#with open("/Users/radhikarawat/PycharmProjects/TCGA_BRCA_exp_HiSeqV2_exon-2015-02-24/genomicMatrix","U") as samples:
    # with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/4out_patient_sample#_header.txt","w") as patient_sample_header:
    #     patient_sample_header.write(samples.readline())


#     with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/4out_gene_base_pairs.txt","w") as bp_list:
#         bps = ""
#         i=0
#         for sample_row in samples:
#             sample_row = sample_row.split('\t')
#             bp_list.write(sample_row[0]+'\t')
#             i+=1
# print "created bps list"
# raise SystemError(0)

#with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/4out_genes.txt", "w") as gene_output_list


genes_dict={}
with open("/Users/radhikarawat/PycharmProjects/TCGA_BRCA_exp_HiSeqV2_exon-2015-02-24/probe","U") as gene_probe:
    for gene_row in gene_probe:
        gene_row=gene_row.split('\t')
        genes_dict[gene_row[0]]=gene_row[1]
bp_list=open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/4out_gene_base_pairs.txt", "U")
bp_row=[]
for info in bp_list:
    bp_row=info.split('\t')
bp_list.close()

with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/4out_genes.txt", "w") as gene_output_list:
            for i in range(0,len(bp_row)):
                if i % 100 == 0: print i
                bps=bp_row[i]
                #print bps
                #if bps in genes_dict.keys():
                try:
                    gene_name=genes_dict[bps]
                    if gene_name=="":
                        gene_output_list.write(bps + '\t')
                    else:gene_output_list.write(gene_name+'\t')
                except KeyError:
                    gene_output_list.write(bps + '\t')



raise SystemError(0)


#     print gene_probe.readline() #remove header row
#     with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/4out_gene_base_pairs.txt","U") as bp_list:
#         #j=0
#         with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/4out_genes.txt", "w") as gene_output_list:
#             for gene_row in gene_probe:
#                 print gene_row
#                 gene_row = gene_row.split('\t')
#                 for bp_list_row in bp_list:
#                     for bps in bp_list_row.split('\t'):
#                         # if i%1000 ==0:
#                         #     #print bps, gene_row[0]
#                         # i+=1
#                         if bps in gene_row:
#                             gene_output_list.write(gene_row[1] + '\t')
#                             print "YES"
#                             #genes += gene_row[1] + '\t'
#                             # j+=1
#                             # if j%100==0:
#                             #     print j, i
#                         else:gene_output_list.write(bps + '\t')
# #             genes += '\n'
# # with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/4out_genes.txt","w") as gene_file:
# #      gene_file.write(genes)