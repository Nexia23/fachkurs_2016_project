chr1_fasta = open("fsa_sequences/S288C_Chromosome I.fsa")
chr1_list = chr1_fasta.read().splitlines()
chr1_header = chr1_list[0]
chr1_list[0] = ""
chr1_seq = "".join(chr1_list)
print(chr1_header)
print(chr1_seq)