import itertools
import re


class Gene(object):

    def __init__(self, mid, chr, sequence, count=0):
        #self.__location = location
        self.__mid = mid
        self.__chr  = chr
        self.__sequence = sequence 
        self.sequence_binding=[0]*len(sequence)


    ###### COMMENT for DATA GROUP #######
    #feel free to replace 'sequence'-information by start-, end-positions and strand (+/-)

    ###### COMMENT FOR REPLICATION_GROUP #######
    #count: 1 for unreplicated gene, 2 for copied gene 

    #@property
    #def location(self):
    #    return self.__location

    @property
    def mid(self):
        return self.__mid

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



def creategenes():

    with open("fsa_sequences/orf_coding.fasta") as orf_fasta:       #open the file, read it and create a list, where each element is a gene with header+sequence
        orf_list = orf_fasta.read().replace("i>", "").replace("sub>", "").replace("->", "").split(">")
    
    orf_list = orf_list[1:] #entfernen des ersten nicht vorhandenden elements
    orf_splitlist = [""]*len(orf_list)  #initialise the new list[gen][header=0 or gen=1]
    
    for i in range(len(orf_list)):  #Trennen von header und sequenz
        orf_splitlist[i] = orf_list[i].split("\n", 1)   

    """
    Gen Sequenz
    """

    for i in range(len(orf_list)):  #replace "\n" 
        orf_splitlist[i][1] = orf_splitlist[i][1].replace("\n", "")


    gene_seq = [x[1] for x in orf_splitlist]
    

    """
    Gen ID
    """

    header_list = [x[0] for x in orf_splitlist]
    header_split = [""]*len(header_list)

    for i in range(len(header_list)):
        header_split[i] = header_list[i].split(" ", 1)


    gene_id = [x[0] for x in header_split]

    """
    Which Chr
    """

    which_chr_list = re.findall(",\s(Chr\s([IVX]+|Mito)|2-micron plasmid)", "".join(header_list)) # find out which chr with regular expressions and put them in a nested list

    which_chr = [x[0] for x in which_chr_list]

    converter_dictionary = {"Chr I" : 1, "Chr II": 2, "Chr III": 3, "Chr IV": 4, "Chr V": 5, "Chr VI": 6, "Chr VII": 7, "Chr VIII": 8, "Chr IX": 9, "Chr X": 10, "Chr XI": 11, "Chr XII": 12, "Chr XIII": 13, "Chr XIV": 14, "Chr XV": 15, "Chr XVI": 16, "Chr Mito": 17, "2-micron plasmid": 18}

    for i in range(len(which_chr)):
        which_chr[i] = converter_dictionary[which_chr[i]]

    gene = {}
    for i in range(len(gene_id)):
        gene[gene_id[i]] = Gene(gene_id[i], which_chr[i], gene_seq[i])

    return gene

    """
    Location
    """

    loc_list = re.findall("([^-][\d]+?)-([\d]+?),", "".join(header_list))

gene = creategenes()
