import sys
import random
import pandas as pd
import numpy as np


class Chromosome(object):
    
    def __init__(self, id , length, sequence, revsequence, fastaname):
        self._id=id
        self._length=length
        self._sequence=sequence
        self._revsequence=revsequence
        self._fastname=fastaname
    
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
        
# Objekte erstellen
#objekte immer klein schrieben und vlt auch abkürzen -> Chromosome1 = chr1

def createchromosome():
    
    chr1=Chromosome(1,0,"deinema","","fsa_sequences/S288C_Chromosome I.fsa")
    chr2=Chromosome(2,0,"","","fsa_sequences/S288C_Chromosome II.fsa")
    chr3=Chromosome(3,0,"","","fsa_sequences/S288C_Chromosome III.fsa")
    chr4=Chromosome(4,0,"","","fsa_sequences/S288C_Chromosome IV.fsa")
    chr5=Chromosome(5,0,"","","fsa_sequences/S288C_Chromosome V.fsa")
    chr6=Chromosome(6,0,"","","fsa_sequences/S288C_Chromosome VI.fsa")
    chr7=Chromosome(7,0,"","","fsa_sequences/S288C_Chromosome VII.fsa")
    chr8=Chromosome(8,0,"","","fsa_sequences/S288C_Chromosome VIII.fsa")
    chr9=Chromosome(9,0,"","","fsa_sequences/S288C_Chromosome IX.fsa")
    chr10=Chromosome(10,0,"","","fsa_sequences/S288C_Chromosome X.fsa")
    chr11=Chromosome(11,0,"","","fsa_sequences/S288C_Chromosome XI.fsa")
    chr12=Chromosome(12,0,"","","fsa_sequences/S288C_Chromosome XII.fsa")
    chr13=Chromosome(13,0,"","","fsa_sequences/S288C_Chromosome XIII.fsa")
    chr14=Chromosome(14,0,"","","fsa_sequences/S288C_Chromosome XIV.fsa")
    chr15=Chromosome(15,0,"","","fsa_sequences/S288C_Chromosome XV.fsa")
    chr16=Chromosome(16,0,"","","fsa_sequences/S288C_Chromosome XVI.fsa")
    chrmito=Chromosome(17,0,"","","fsa_sequences/S288C_Chromosome Mito.fsa")

    chr_list=[chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chrmito]
    
    return chr_list

    
getchr=createchromosome()
print(getchr[0].sequence)
    



