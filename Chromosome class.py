import sys
import random
import pandas as pd
import numpy as np


class Chromosome(object):
    
    def __init__(self, id ,fastaname):
        self._id=id
        self._revsequence= None
        self._fastaname=fastaname
        
        #einlesen des files, löschen des headers, zusammenfügen der einzelnen zeilen als einzelnen string der als sequence abgespeichert wird
        with open(self._fastaname) as chr_fasta:
            chr_list = chr_fasta.read().splitlines() 

        chr_list[0] = ""
        self._sequence = "".join(chr_list)

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
    @sequence.setter
    def sequence(self, value):
        if not isinstance(value, str):
            raise TypeError("Sequence must be a String.")
        self._sequence = value
        
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



# Objekte erstellen
#objekte immer klein schrieben und vlt auch abkürzen -> Chromosome1 = chr1

def createchromosomes():
    
    chr1=Chromosome(1,"fsa_sequences/S288C_Chromosome I.fsa")
    chr2=Chromosome(2,"fsa_sequences/S288C_Chromosome II.fsa")
    chr3=Chromosome(3,"fsa_sequences/S288C_Chromosome III.fsa")
    chr4=Chromosome(4,"fsa_sequences/S288C_Chromosome IV.fsa")
    chr5=Chromosome(5,"fsa_sequences/S288C_Chromosome V.fsa")
    chr6=Chromosome(6,"fsa_sequences/S288C_Chromosome VI.fsa")
    chr7=Chromosome(7,"fsa_sequences/S288C_Chromosome VII.fsa")
    chr8=Chromosome(8,"fsa_sequences/S288C_Chromosome VIII.fsa")
    chr9=Chromosome(9,"fsa_sequences/S288C_Chromosome IX.fsa")
    chr10=Chromosome(10,"fsa_sequences/S288C_Chromosome X.fsa")
    chr11=Chromosome(11,"fsa_sequences/S288C_Chromosome XI.fsa")
    chr12=Chromosome(12,"fsa_sequences/S288C_Chromosome XII.fsa")
    chr13=Chromosome(13,"fsa_sequences/S288C_Chromosome XIII.fsa")
    chr14=Chromosome(14,"fsa_sequences/S288C_Chromosome XIV.fsa")
    chr15=Chromosome(15,"fsa_sequences/S288C_Chromosome XV.fsa")
    chr16=Chromosome(16,"fsa_sequences/S288C_Chromosome XVI.fsa")
    chrmito=Chromosome(17,"fsa_sequences/S288C_Chromosome Mito.fsa")

    chr_list=[chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chrmito]
    
    return chr_list

def createwholegenome(chr_list):

    whole_genome = ""
    for i in range(len(chr_list)):
        whole_genome += whole_genome + chr_list[i].sequence

    return whole_genome


chr_list = createchromosomes()
createwholegenome(chr_list)



