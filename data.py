import re
import sys
import pandas as pd
import numpy 


class Chromosome:
    """
    An object for the chromosome with the attributes:

    id -> gives out the chromosome in integer form 1-16, 17 = mitochondrial chromosome
    fastaname -> is only important for the initiation of the object and used by the function createchromosomes
    sequence -> sequence of the whole chromosome
    revsequence -> the reverse sequence of the whole chromosome

    possible operations are:
    gene.id -> gives the chromosome number 1-16, 17 = mitochondrial chromosome
    gene.sequence -> gives the sequence of the chromosome
    gene.revsequence -> gives out the reverse sequence of the chromosome
    """
    def __init__(self, id,  arf, fastaname):
        self._id=id
        self._fastaname=fastaname
        self._arf = arf
        self.binding_molecules=[[],[]]  #list with tuples of start and end positions of occupied regions in [0] and the binding molecule in [1]
        self.replication_ori_bound = False
        
        #read in the file, delete the headers, concatenate the single lines of the sequence and store them as a single string in sequence
        with open(self._fastaname) as chr_fasta:
            chr_list = chr_fasta.read().splitlines() 
        chr_list[0] = ""
        #self._sequence = "".join(chr_list)
        


        sequence_str = "".join(chr_list)
        self._sequence = list(sequence_str)
        


        #generate the reverse sequence

        converter_dictionary = {"A" : "T", "T": "A", "C": "G", "G": "C"}

        sequence_str = list(sequence_str)
        for i in range(len(sequence_str)):
            sequence_str[i] = converter_dictionary[sequence_str[i]]
        self._revsequence = sequence_str

    #add befehl   
    def __add__(self,chromosome):
        if not isinstance(chromosome,Chromosome):
            raise TypeError
        self.sequence=self.sequence+ chromosome.sequence
        self.revsequence=self.revsequence+ chromosome.revsequence
        return self
    
    #getter für id, sequence & revsequence

    @property
    def id(self):
        return self._id
    @id.setter
    def id(self, value):
        if not isinstance(value, int):
            raise TypeError("ID must be an Integer.")
        self._id = value
        
    @property
    def sequence(self):
        return self._sequence

    @property
    def arf(self):
        return self._arf
    
        
    @property
    def revsequence(self):
        return self._revsequence
    @revsequence.setter
    def revsequence(self, value):
        if not isinstance(value, str):
            raise TypeError("RevSequence must be a String.")
        self._revsequence = value
        
    @property
    def fastaname(self):
        return self._fastaname


#####this method will check for bound regions in the chromosome######

    def chromosome_bound(self, range):

        if isinstance(range, list):
            #range[0]: startposition of test, range[1]: endposition of test
            iter=0
            bound_stuff=[]
            bound_found=False
            for tuples in self.binding_molecules[0]:

                if bound_found==False:
                    if tuples[0]>range[1] and iter==0:
                        return None
                    elif tuples[0]>range[0] and tuples[0]>=range[1]:
                        return None
                    elif tuples[1]<range[0] and tuples[1]<range[0]:
                        pass
                    elif tuples[0]<=range[0] and tuples[1]>=range[1]:
                        bound_stuff.append(self.binding_molecules[1][self.binding_molecules[0].index(tuples)])
                        return bound_stuff
                    elif tuples[1]<range[1]:
                        bound_stuff.append(self.binding_molecules[1][self.binding_molecules[0].index(tuples)])
                        bound_found=True
                    
                else:
                    if tuples[0]>range[1]:
                        return bound_stuff
                    elif tuples[1]>range[1]:
                        bound_stuff.append(self.binding_molecules[1][self.binding_molecules[0].index(tuples)])
                        return bound_stuff
                    else:
                        bound_stuff.append(self.binding_molecules[1][self.binding_molecules[0].index(tuples)])
                iter+=1
            return bound_stuff



        elif isinstance(range, int):
            for tuples in self.binding_molecules[0]:
                if tuples[0]<=range and tuples[1]>=range:
                    return self.binding_molecules[1][self.binding_molecules[0].index(tuples)]
                elif tuples[0]>range:
                    break
           
        else:
            print("Argument type not expected (list or int)")

    #### method needed which stores a tuple of start&end and bound molecule in the ordered (!!!) list bindnig_molecules
    def bind_to_chrom(start, end):
        pass

class Gene:
    """
    an object for the gene sequence with the attributes:
    id -> gene id
    name -> gene name
    chr -> which chromosome the gene is on in integer form, 1-15 und 17 = mitochondrial chromosome and 18 = 2-micron plasmid
    sequence -> sequence of the dna which is transcribed
    count ->
    sequence_binding -> records if the gene is bound (=1) or currently not bound (=0), important for transcription vs replikation (i think)
    transrate -> float of the transkription rate, if none was found a median of all values was used, set to per second
    halflive -> float of the halflive time, used for decay of the mRNA, if none was found, a median of all values was used set to per second

    possible operations are:
    gene.mid -> provides the gene name
    gene.chr -> provides the chromosome on which the gene is on
    gene.sequence -> provides the sequence of the gene
    gene.name -> provides the name of the gene
    gene.translate -> provides the transcription rate
    gene.halflive -> provides the halflive of the corresponding mRNA
    """

    def __init__(self, mid, name, chr, sequence, location, transrate, halflive, count=0):
        
        self.__mid = mid
        self.__name = name
        self.__location = location
        self.__chr  = chr
        self.__sequence = sequence
        self.__transrate = transrate
        self.__halflive = halflive 
        self.sequence_binding=[0]*len(sequence)
        self.rnas_transcribed=0                             #number of transcribed RNAs during one transcription-process
        self.pol_on_gen = []                                #nucleotides transcribed by polymerases on the gene
        self.rate=32
                        
        
        if numpy.isnan(self.__transrate):
            self.__transrate = 0.00123056
        if numpy.isnan(self.__halflive):
            self.__halflive = 0.262
    ###### COMMENT for DATA GROUP #######
    #feel free to replace 'sequence'-information by start-, end-positions and strand (+/-)

    ###### COMMENT FOR REPLICATION_GROUP #######
    #count: 1 for unreplicated gene, 2 for copied gene 

    @property
    def mid(self):
        return self.__mid

    @property
    def chr(self):
        return self.__chr

    @property
    def sequence(self):
        return self.__sequence
        
    @property
    def name(self):
        return self.__name

    @property
    def location(self):
        return self.__location

    @property
    def transrate(self):
        return self.__transrate
    
    @property
    def halflive(self):
        return self.__halflive


        
    @sequence.setter
    def sequence(self, value):
        if not isinstance(value, str):
            raise Exception("sequence must be a string")
            # TODO: check for valid nucleotides here
        self.__sequence = value.upper()


def createwholegenome(chr_list):

    whole_genome = ""
    for i in range(len(chr_list)):
        whole_genome += whole_genome + chr_list[i].sequence

    return whole_genome

def createchromosomes():
    """
    arf = [[]]

    ori_raw = pd.read_excel("fsa_sequences/ORI.xls")
    
    db_chr1 = ori_raw[ori_raw.ix[:,0]=="chr1"]
    arf_chr1 = 
    """


    chr1=mol.Chromosome("Chr 1","fsa_sequences/S288C_Chromosome I.fsa")
    chr2=mol.Chromosome("Chr 2","fsa_sequences/S288C_Chromosome II.fsa")
    chr3=mol.Chromosome("Chr 3","fsa_sequences/S288C_Chromosome III.fsa")
    chr4=mol.Chromosome("Chr 4","fsa_sequences/S288C_Chromosome IV.fsa")
    chr5=mol.Chromosome("Chr 5","fsa_sequences/S288C_Chromosome V.fsa")
    chr6=mol.Chromosome("Chr 6","fsa_sequences/S288C_Chromosome VI.fsa")
    chr7=mol.Chromosome("Chr 7","fsa_sequences/S288C_Chromosome VII.fsa")
    chr8=mol.Chromosome("Chr 8","fsa_sequences/S288C_Chromosome VIII.fsa")
    chr9=mol.Chromosome("Chr 9","fsa_sequences/S288C_Chromosome IX.fsa")
    chr10=mol.Chromosome("Chr 10","fsa_sequences/S288C_Chromosome X.fsa")
    chr11=mol.Chromosome("Chr 11","fsa_sequences/S288C_Chromosome XI.fsa")
    chr12=mol.Chromosome("Chr 12","fsa_sequences/S288C_Chromosome XII.fsa")
    chr13=mol.Chromosome("Chr 13","fsa_sequences/S288C_Chromosome XIII.fsa")
    chr14=mol.Chromosome("Chr 14","fsa_sequences/S288C_Chromosome XIV.fsa")
    chr15=mol.Chromosome("Chr 15","fsa_sequences/S288C_Chromosome XV.fsa")
    chr16=mol.Chromosome("Chr 16","fsa_sequences/S288C_Chromosome XVI.fsa")
    chrmito=mol.Chromosome("Chr 17","fsa_sequences/S288C_Chromosome Mito.fsa")

    chr_list=[chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chrmito]
        
    return chr_list


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


    """
    gene_name
    """

    gene_name_list = [""]*len(header_list)

    for i in range(len(header_split)):
        gene_name_list[i] = header_split[i][1].split(" ", 1)

    gene_name = [x[0] for x in gene_name_list]



    """
    Location
    """

    loc_list = [0]*len(header_list)

    for i in range(len(header_list)):   #creates a list where loc_list[0] = all the tuples of locations for gene 0 in str form
        loc_list[i] = re.findall("([^-][\d]+?)-([\d]+?),", header_list[i])

    
    """
    transrate
    """
    transcription_raw = pd.read_excel("fsa_sequences/msb2010112-s1.xls")
    
    
    transcription_raw[transcription_raw.ix[:,1]<0]=0    #negative werte auf 0 setzen
    transcription_raw.ix[:,1]=transcription_raw.ix[:,1].divide(9000)    #in sekunden umrechnen

    transcription = transcription_raw.ix[gene_id, 1]

    """
    halflive
    """
    halflive_raw = pd.read_excel("fsa_sequences/msb2010112-s2.xls")

    halflive_raw.ix[:,1] = halflive_raw.ix[:,1].divide(60)

    halflive = halflive_raw.ix[gene_id, 1]

    """
    return
    """
    gene = {}
    for i in range(len(gene_id)):
        gene[gene_id[i]] = mol.Gene(gene_id[i], gene_name[i], which_chr[i], gene_seq[i], loc_list[i], transcription[i], halflive[i])

    return gene