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

    #open the file
    orf_fasta = open("fsa_sequences/orf_coding.fasta")
    #read the file in and create a list, where each element is a gen with header+sequence
    orf_list = orf_fasta.read()
    orf_list = orf_list.replace("i>", "").replace("sub>", "").replace("->", "")
    orf_list = orf_list.split(">")
    #create a list of length 2 with the first element being the header and the second being the sequence (still with line breaks (use join later))
    #combine with loop to create len(orf_list) lists
    #initialise the new list[gen][header=0 or gen=1]
    orf_splitlist = [""]*len(orf_list)

    #Trennen von header und sequenz
    for i in range(0,len(orf_list)):
        orf_splitlist[i] = orf_list[i].split("\n", 1)
    
    #zusammenfügen der Sequenz 
    for i in range(1, len(orf_list)):
        orf_splitlist[i][1] = "".join(orf_splitlist[i][1].rsplit())

    #entfernen des ersten nicht vorhandenden elements
    orf_splitlist = orf_splitlist[1:]

    #unnesten
    orf_splitlist_unnested = [x for i in orf_splitlist for x in i]

  
    gene_seq = orf_splitlist_unnested[1::2]
    

    #Gen ID

    header_list = orf_splitlist_unnested[0::2]

    header_split = [""]*len(header_list)

    for i in range(0,len(header_list)):
        header_split[i] = header_list[i].split(" ", 1)

    header_split_unnested = [x for i in header_split for x in i]
    gen_id = header_split_unnested[0::2]

    return gene_seq, gen_id, 


