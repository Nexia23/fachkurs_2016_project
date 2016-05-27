import itertools

class Gene(BioMoleculeCount):

    def __init__(self, mid, name, location, chr, sequence, strand, count=0):
        super().__init__(mid, name, count)
        self.__location = location
        self.__chr  = chr
        self.__sequence = sequence 
        self.sequence_binding=[0]*len(sequence)


    ###### COMMENT for DATA GROUP #######
    #feel free to replace 'sequence'-information by start-, end-positions and strand (+/-)

    ###### COMMENT FOR REPLICATION_GROUP #######
    #count: 1 for unreplicated gene, 2 for copied gene 

    @property
    def location(self):
        return self.__location

    @property
    def chr(self):
        return self.__chr

    @property
    def sequence(self):
        return self.__sequence
        
    @sequence.setter
    def sequence(self, value):
        if not isinstance(value, str):
            raise Exception("sequence must be a string")
            # TODO: check for valid nucleotides here
        self.__sequence = value.upper()



def creategene():

	with open("fsa_sequences/orf_coding.fasta") as orf_fasta:		#open the file, read it and create a list, where each element is a gene with header+sequence
		orf_list = orf_fasta.read().replace("i>", "").replace("sub>", "").replace("->", "").split(">")
	

    orf_splitlist = [""]*len(orf_list)	#initialise the new list[gen][header=0 or gen=1]
    orf_list = orf_list[1:] #entfernen des ersten nicht vorhandenden elements
    
    for i in range(len(orf_list)):	#Trennen von header und sequenz
        orf_splitlist[i] = orf_list[i].split("\n", 1)	
    
    for i in range(len(orf_list)):	#replace "\n" 
        orf_splitlist[i][1] = orf_splitlist[i][1].replace("\n", "")

    orf_splitlist_unnested = [x for i in orf_splitlist for x in i]	#unnesten

    gene_seq = orf_splitlist_unnested[1::2]
    

    #Gen ID

    header_list = orf_splitlist_unnested[0::2]
    header_split = [""]*len(header_list)

    for i in range(len(header_list)):
        header_split[i] = header_list[i].split(" ", 1)

    header_split_unnested = [x for i in header_split for x in i]
    gene_id = header_split_unnested[0::2]

    return gene_seq, gene_id, 


