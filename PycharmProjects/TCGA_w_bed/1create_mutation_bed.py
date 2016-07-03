#ene name	Accession Number	Gene CDS length	HGNC ID	Sample name	ID_sample	ID_tumour	Primary site	Site subtype 1	Site subtype 2	Site subtype 3	Primary histology	Histology subtype 1	Histology subtype 2	Histology subtype 3	Genome-wide screen	Mutation ID	Mutation CDS	Mutation AA	Mutation Description	Mutation zygosity	LOH	GRCh	Mutation genome position	Mutation strand	SNP	Resistance Mutation	FATHMM prediction	FATHMM score	Mutation somatic status	Pubmed_PMID	ID_STUDY	Sample source	Tumour origin	Age


with open("/Users/radhikarawat/PycharmProjects/CosmicMutantExport.tsv","U") as mutation_data_file:
    with open("mutation_bed.bed","w") as f:
        for row in mutation_data_file:
            row=row.split("\t")
            exonCount=row[8]
            if row[22]=="38":
                if "breast" in row:

                    mutchrom,bps = row[23].split(":")
                    mutchrom="chr"+mutchrom
                    mutstart,mutend = bps.split("-")
                    name=row[1]
                    gene_name=row[0]
                    TCGA_info=row[4]
                    #score=row[11]

                    f.write ("%s\t%s\t%s\t%s\t%s\t%s\n" % (mutchrom, mutstart, mutend,name,gene_name,TCGA_info))