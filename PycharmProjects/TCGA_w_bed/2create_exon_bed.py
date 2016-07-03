#create bed file from exon data

##import exon data file


##bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames

bracket=3#
##format row
with open("/Users/radhikarawat/PycharmProjects/TCGA_w_bed/transcripts.txt","U") as exon_data_file:
    with open("exon_splice_sites_bed.bed","w") as f:
        for row in exon_data_file:
            row=row.split("\t")

            if "e" not in row[8]:
                exonCount=int(row[8])
                if exonCount>1:
                    chrom = row[2]
                    name=row[1]
                    name2=row[12]
                    score=row[11]
                    exonStarts=row[9].split(",")
                    exonEnds=row[10].split(",")
                    numexons=len(exonEnds)
                    for i in range(0,numexons):
                        if exonStarts[i] !="":
                            if exonEnds[i]!="":
                                SpliceSiteStart1=str(int(exonStarts[i])-bracket)
                                SpliceSiteEnd1=str(int(exonStarts[i])+bracket)
                                if int(SpliceSiteStart1) >0:

                                    f.write ("%s\t%s\t%s\t%s\n" % (chrom, SpliceSiteStart1, SpliceSiteEnd1,name+'_'+str(i)))
                                SpliceSiteStart2 = str(int(exonEnds[i]) - bracket)
                                SpliceSiteEnd2 = str(int(exonEnds[i]) + bracket)
                                if int(SpliceSiteStart2) > 0:
                                    f.write("%s\t%s\t%s\t%s\t\n" % (chrom, SpliceSiteStart2, SpliceSiteEnd2,
                                                                    name+'_'+str(i)))



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
