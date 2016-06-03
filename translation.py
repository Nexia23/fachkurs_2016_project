import processes
import molecules
import numpy as np
#import random
import copy


class Translation(processes.Process):
    """
    Translation is instantiated in the Cell to produce proteins.

    Defines Translation process. It iterates over all ribosomes and decides what
    they should do. They either bind to mRNA or elongate/terminate a protein if
    they are already bound.

    """
    code = dict([('UCA', 'S'), ('UCG', 'S'), ('UCC', 'S'), ('UCU', 'S'),
                 ('UUU', 'F'), ('UUC', 'F'), ('UUA', 'L'), ('UUG', 'L'),
                 ('UAU', 'Y'), ('UAC', 'Y'), ('UAA', '*'), ('UAG', '*'),
                 ('UGU', 'C'), ('UGC', 'C'), ('UGA', '*'), ('UGG', 'W'),
                 ('CUA', 'L'), ('CUG', 'L'), ('CUC', 'L'), ('CUU', 'L'),
                 ('CCA', 'P'), ('CCG', 'P'), ('CCC', 'P'), ('CCU', 'P'),
                 ('CAU', 'H'), ('CAC', 'H'), ('CAA', 'Q'), ('CAG', 'Q'),
                 ('CGA', 'R'), ('CGG', 'R'), ('CGC', 'R'), ('CGU', 'R'),
                 ('AUU', 'I'), ('AUC', 'I'), ('AUA', 'I'), ('AUG', 'M'),
                 ('ACA', 'T'), ('ACG', 'T'), ('ACC', 'T'), ('ACU', 'T'),
                 ('AAU', 'N'), ('AAC', 'N'), ('AAA', 'K'), ('AAG', 'K'),
                 ('AGU', 'S'), ('AGC', 'S'), ('AGA', 'R'), ('AGG', 'R'),
                 ('GUA', 'V'), ('GUG', 'V'), ('GUC', 'V'), ('GUU', 'V'),
                 ('GCA', 'A'), ('GCG', 'A'), ('GCC', 'A'), ('GCU', 'A'),
                 ('GAU', 'D'), ('GAC', 'D'), ('GAA', 'E'), ('GAG', 'E'),
                 ('GGA', 'G'), ('GGG', 'G'), ('GGC', 'G'), ('GGU', 'G')])

    def __init__(self, id, name, debug=True):
        # call the constructor of the base class (processes.Process in this case)
        super().__init__(id, name)
        self.__initiate_ribosomes()     #initiating one Ribosome
        self.debug=debug

    def __repr__(self):
        # representation of object: should be unambiguous!
        return ','.join([self.name, str(type(self))])

    def __str__(self):
        # representation of object: to be nice and printable
        return ','.join([self.name, str(type(self))])


    def __initiate_ribosomes(self):
        self.ribosomes = molecules.Ribosome('Ribosomes', 'Ribosome', 1)


    def update(self, model):
        """
        Update all mrnas and translate proteins.
        """
        # enzymes_id -> the only one initialized ribosome
       
        self.ribosomes = model.states[list(self.enzyme_ids)[0]] # call in the dictionary states for the Ribsosome-object

        mRNA_list = [] # list of current mRNA objects
        for state in model.states:

            if isinstance(model.states[state], list): # check whether model.states is a list
            # if yes --> iterate over all elements and check whether they are mRNA
                for elem in model.states[state]:
                    if isinstance (elem, molecules.MRNA):
                        # check whether first start codon is in frame or not. If not: append 1 or 2 nucleotides to 5'-end of mRNA. Currently not needed, since mRNAs always start with AUG.
                        bp_before_AUG = elem.sequence.find('AUG')
                        if (bp_before_AUG != -1) and (bp_before_AUG%3 != 0): 
                            if bp_before_AUG%3 == 1:
                                elem.sequence = "UU" + elem.sequence
                            if bp_before_AUG%3 ==2:
                                elem.sequence = "U" + elem.sequence

                        mRNA_list.append(elem) # store current mRNAs for translation
                        
                        
        #print(mRNA_list)
        for i in range (10):
                          
            np.random.shuffle(mRNA_list)         # List wird durch shuffle gemischt
            for mrna in mRNA_list:           # substrate should be a list for ids of all mrnas
                                                          
                # ribosoms work: bind, move, initialise, elongate
                if mrna.sequence_triplet_binding[0] == 0:  # check if 1st codon is empty
                   
                    self.bind(mrna)  # bind to first position on mRNA
                
                prot = self.move(mrna)

                if isinstance(prot, molecules.Protein):     # storing the protein in the states-dictionary
                    #print('creation of protein')
                    if prot.name in model.states:
                        model.states[prot.name].append(prot)
                    else:
                        model.states[prot.name] = [prot]
                    
    def bind(self, mrna):
        """
        Bind to 5'-end of mRNA --> initiate / move without protein synthesis
        
        """
        
        if np.random.poisson(self.ribosomes.count) > 1:  # check if ribosome can bind
            mrna.sequence_triplet_binding[0] = 'R'          # !Jens! codon wird 'R', Teil des Ribosomes binden an 0. Stelle binden!

            if mrna.sequence[0:3]=='AUG':                   # if first codon is START codon -> initiate
                self.initiate(mrna, 0)
                
            self.ribosomes.count -= 1  # remove bound ribosome from list of free ribosomes

    def move(self, mrna):
        """
        move not initiated ribosomes and also other ribosomes

        """

        ribo_pos = [j for j, val in enumerate(mrna.sequence_triplet_binding) if val=='R'] #liste die R-positionen speichert also den index in triplet_binding 
        prot_pos = [j for j, val in enumerate(mrna.sequence_triplet_binding) if isinstance(val, molecules.Protein)] #liste der Protein-Positionen
        
        pos = ribo_pos+prot_pos # Liste aller R- und Proteinpositionen
    
        np.random.shuffle(pos) # shuffles the list so random order is created
        
        for i in pos:  # iterate through all ribosome/protein positions  
            
            if i == len(mrna.sequence_triplet_binding)-1: # wenn das Ende der mRNA erreicht ist (letztes codon)
                self.setfree(mrna, i)
                self.ribosomes.count += 1   # und erhöht die Menge freier Ribosomen um 1
                continue                                                           
            
            if mrna.sequence_triplet_binding[i+1] == 1:  # besetzt also muss es warten
                
                continue 
        
            # falls nächste Position nicht besetzt, entweder initiate und/oder elongate
            if mrna.sequence_triplet_binding[i] == 'R': # falls R noch nicht initiiert...
                if mrna[i * 3:i * 3 + 3] == 'AUG': #...und das nächste Codon das Startcodon ist
                    self.initiate(mrna, i)
                else: # Bewegung ohne Initiation
                    self.occupy(mrna, i+1)
                    continue
            elif not isinstance(mrna.sequence_triplet_binding[i], molecules.Protein):
                raise Exception("Problem: we should have a protein at this position")

            return self.elongate( mrna, i)
        return 0

    def initiate(self, mrna, i):
        """
        Initiate translation, create protein object.

        @type mrna: MRNA
        """
        # !Jens! Da muss jemand kraeftig umbauen
        mrna.sequence_triplet_binding[i] = molecules.Protein("Protein_{}".format(mrna.mid),
                                                             "Protein_{0}".format(mrna.name.split("_")[-1]),
                                                             "",)
        
    def elongate(self, mrna, i):
        """
        Elongate the new protein by the correct amino acid. Check if an
        MRNA is bound and if ribosome can move to next codon.
        Terminate if the ribosome reaches a STOP codon.

        @type return: Protein or False
        """
        print('elongation')
        if isinstance(mrna.sequence_triplet_binding[i], molecules.Protein):  

            codon = mrna[i * 3:i * 3 + 3]
            aa = self.code[codon]

            if aa == "*":  # terminate at stop codon
                return self.terminate(mrna, i)

            if i + 1 >= len(mrna.sequence_triplet_binding):
                return self.terminate(mrna, i)  # terminate if mrna ends

            if mrna.sequence_triplet_binding[i + 1] == 0:  # if the next rna position is free

                mrna.sequence_triplet_binding[i] += aa  
                print(mrna.sequence_triplet_binding[i].sequence)   
                self.occupy(mrna, i+1)                      # nächste Stelle besetzten durch ocupyfunction
              
        return 0
    def terminate(self, mrna, i):
        """
        Splits the ribosome/MRNA complex and returns a protein.
        """
      
        protein = mrna.sequence_triplet_binding[i]  # bound mRNA
        self.setfree(mrna,i)
        self.ribosomes.count += 1
        # put proteins in lists (state)
        return protein  

    def setfree(self, mrna, i):      
        if i < 10:
            k=0
            while k <= i:                             #entkoppelt alles von anfang an
                mrna.sequence_triplet_binding[k] = 0
                k += 1 

        elif i >=10:

            k=0
            while k <=10:                              # Zehnstellen davor werden entkoppelt
                mrna.sequence_triplet_binding[i-k] = 0
                k+=1

    def occupy(self, mrna, i): 
        """
        Movement of Ribosome OR Protein to next codon. 
        !!! Moves from position i-1 to position i TODO:Compatibility to elongate (i-1 --> i)
        Occupies the last 9 codons and frees up the 10th (backwards)
        """
       
        if isinstance(mrna.sequence_triplet_binding[i-1], molecules.Protein):
            mrna.sequence_triplet_binding[i] = mrna.sequence_triplet_binding[i-1]
        else:
            mrna.sequence_triplet_binding[i]='R'
        if i < 10:
            for k in range(i):                          # besetzt alle Stellen ab der ersten
                mrna.sequence_triplet_binding[k] = 1         
        else:
            mrna.sequence_triplet_binding[i-10] = 0        # Stelle die frei wird
            for k in range(9):
                mrna.sequence_triplet_binding[i-(k+1)] = 1     # R gefolgt von 9 Einsen(ribosom) 
        