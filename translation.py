import processes
import molecules
import numpy as np
#import random
#import copy


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

    def __init__(self, id, name):
        # call the constructor of the base class (processes.Process in this case)
        super().__init__(id, name)
        self.__initiate_ribosomes()     #initiating one Ribosome

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

        for i in range (10):

            for mrna_id in np.random.shuffle(self.substrate_ids):           # substrate should be a list for ids of all mrnas
                prot = None                                                 # initialise prot variable
                mrna = model.states[mrna_id]                                # object of mRNA

                # ribosoms work: bind, move, initialise, elongate
                if mrna.sequence_triplet_binding[0] == 0:  # check if 1st codon is empty
                    self.bind(mrna)  # bind to first position on mRNA
                self.move(mrna)  # move ribosome that has not yet initiated translation
                    
                # put proteins in lists (state)
                if isinstance(prot, molecules.Protein):     # storing the protein in the states-dictionary
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
                self.initiate(mrna)

            self.ribosomes.count -= 1  # remove bound ribosome from list of free ribosomes

    def move(self, mrna):
        """
        move not initiated ribosomes and also other ribosomes
        

        """
        ribo_pos = [j for j, val in enumerate(mrna.sequence_triplet_binding) if val=='R'] #liste die R-positionen speichert also den index in triplet_binding 
        
        for i in np.random.shuffle(ribo_pos):  # iterate through all codons randomly (1 - occupied, 'R' - not initiated ribosome, 0 - free)  
           
            if ribo_pos[i] == len(mrna.sequence_triplet_binding)-1: # wenn das Ende der mRNA erreicht ist (letztes codon)
                self.entkoppeln(mrna, ribo_pos[i])
                self.ribosomes.count += 1   # und erhöht die Menge freier Ribosomen um 1                                                           
            
            
            if mrna.sequence_triplet_binding[ribo_pos[i]+1] == 1:  # besetzt also muss es warten
                break 
            
            elif mrna.sequence_triplet_binding[ribo_pos[i]+1]== 0:   # ...und die nächste Stelle frei ist
                if mrna.sequence[ribo_pos[i]+1] == 'AUG':     # ...und die nächste Stelle ein Startcodon ist
                    self.initiate(mrna, ribo_pos[i]+1) # gehe zu initiate !Jens!: Siehe Issue #10
                    self.elongate(mrna)  # danach kann keines mehr nur bewegt werden                    
                else:   # wenn die nächste Stelle kein Startcodon ist
                    self.occupy(mrna, ribo_pos[i]+1)  

            elif isinstance(ribosome, molecules.Protein):   # falls ein Protein synthetisiert wird
                self.elongate(mrna)  # gehe zu elongate
    

    def initiate(self, mrna, i):
        """
        Initiate translation, create protein object.

        @type mrna: MRNA
        """
        # !Jens! Da muss jemand kraeftig umbauen
        mrna.sequence_triplet_binding[i] = molecules.Protein("Protein_{}".format(mrna.mid),
                                                             "Protein_{0}".format(mrna.name.split("_")[-1]),
                                                             "",)

    def elongate(self, mrna):
        """
        Elongate the new protein by the correct amino acid. Check if an
        MRNA is bound and if ribosome can move to next codon.
        Terminate if the ribosome reaches a STOP codon.

        @type return: Protein or False
        """
        #rd_pool = copy.deepcopy(mrna.sequence_triplet_binding) #kopiere Liste aus der wir unser ribosome waehlen
        #rd_ribo = random.choice(rd_pool) #waehle Stelle die elongiert werden soll
        #rd_pool = rd_pool  #hier soll das element aus der waehlmenge raus genommen werden um nicht doppelt zu waehlen

        for i, ribosome in enumerate(mrna.sequence_triplet_binding):


                if isinstance(ribosome, molecules.Protein):

                    codon = mrna[i * 3:i * 3 + 3]
                    aa = self.code[codon]

                    if aa == "*":  # terminate at stop codon
                        return self.terminate(mrna, i)

                    if i + 1 >= len(mrna.sequence_triplet_binding):
                        return self.terminate(mrna, i)  # terminate if mrna ends

                    if not mrna.sequence_triplet_binding[i + 1]:  # if the next rna position is free
                        mrna.sequence_triplet_binding[i] += aa
                        self.occupy(mrna, i+1)                      # nächste Stelle besetzten durch ocupyfunction

        return 0
    
    def terminate(self, mrna, i):
        """
        Splits the ribosome/MRNA complex and returns a protein.
        """
        protein = mrna.sequence_triplet_binding[i]  # bound mRNA
        self.entkoppeln(mrna,i)
        self.ribosomes.count += 1
        return protein
    

    def entkoppeln(self, mrna, i):      
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


    def occupy(self, mrna, i):        #Funktion die Besetzung einrichtet

        if i < 10:
            k=0
            while k <= i:                          # besetzt alle Stellen ab der ersten
                mrna.sequence_triplet_binding[k] = 1
                k += 1  

        elif i >= 10:
            mrna.sequence_triplet_binding[i-10] = 0        # Stelle die frei wird
            k=0
            mrna.sequence_triplet_binding[i] = 'R' 
            while k <= 9:
                mrna.sequence_triplet_binding[i-k] = 1     # R gefolgt von 9 Einsen(ribosom) 
                k += 1

            

