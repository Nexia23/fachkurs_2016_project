#read_in

#open the fasta file
chr1_fasta = open("fsa_sequences/S288C_Chromosome I.fsa")

#create a list where each element is one line of the fasta file
chr1_list = chr1_fasta.read().splitlines()

#put the header of the fasta file in an extra string
chr1_header = chr1_list[0]

#delete the header
chr1_list[0] = ""

#create a string with the whole sequence
chr1_seq = "".join(chr1_list)
