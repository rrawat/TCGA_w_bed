from functions import *
create_mutation_bed()
create_exon_bed(bracket=1)
bed_intersect_path=intersect_bed_files(bracket=1)
process_splice_mutations(bed_intersect_path)
#create_patient_sample_header()
#gene_dictionary=create_gene_dict()
#add_gene_names_to_matrix(gene_dictionary)



