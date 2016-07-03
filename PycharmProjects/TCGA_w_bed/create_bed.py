#create bed file from exon data

##import exon data file


##bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames


##format row
with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/transcripts.txt","U") as exon_data_file:
    with open("exon_bed.bed","w") as f:
        for row in exon_data_file:
            row=row.split("\t")
            exonCount=row[8]
            if exonCount>1:
                chrom = row[2]
                name=row[1]
                name2=row[12]
                score=row[11]
                exonStarts=row[9].split(",")
                exonEnds=row[10].split(",")
                numexons=len(exonEnds)
                for i in range(0,numexons):
                    f.write ("%s\t%s\t%s\t%s\n" % (chrom, exonStarts[i], exonEnds[i],name))


# The first three required BED fields are:
#
# chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
# chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
# chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.


##write row


##close exon data file, bed file



#create bed file from mutation data

##import mutation data file

##format row

##write row


##close mutation data file, bed file
