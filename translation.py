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
        if self.debug:
            print('update')
        self.ribosomes = model.states[list(self.enzyme_ids)[0]] # call in the dictionary states for the Ribsosome-object

        sub_kopy = list(self.substrate_ids)        #in Liste umschreiben damit Zugriff auf Index moeglich
        print(sub_kopy)
        for i in range (10):
            
            if self.debug:
                print('die 10 jaegerlein')
                
            #np.random.shuffle(sub_kopy)         # List wird durch shuffle gemischt
            for mrna_id in sub_kopy:           # substrate should be a list for ids of all mrnas
                                                            # initialise prot variable
                                                # object of mRNA
                mrna = model.states[mrna_id]
                # ribosoms work: bind, move, initialise, elongate
                if mrna.sequence_triplet_binding[0] == 0:  # check if 1st codon is empty
                    print('pos 0')
                    print(mrna.sequence_triplet_binding[0])
                    self.bind(mrna)  # bind to first position on mRNA
                
                self.move(model, mrna)
                
    def bind(self, mrna):
        """
        Bind to 5'-end of mRNA --> initiate / move without protein synthesis
        
        """
        
        if np.random.poisson(self.ribosomes.count) > 1:  # check if ribosome can bind
            mrna.sequence_triplet_binding[0] = 'R'          # !Jens! codon wird 'R', Teil des Ribosomes binden an 0. Stelle binden!
            if self.debug:
                print("binding")
                print(mrna.sequence_triplet_binding)

            if mrna.sequence[0:3]=='AUG':                   # if first codon is START codon -> initiate
                self.initiate(mrna, 1)
                
            self.ribosomes.count -= 1  # remove bound ribosome from list of free ribosomes

    def move(self, model, mrna):
        """
        move not initiated ribosomes and also other ribosomes
        

        """
        ribo_pos = [j for j, val in enumerate(mrna.sequence_triplet_binding) if val=='R'] #liste die R-positionen speichert also den index in triplet_binding 
        prot_pos = [j for j, val in enumerate(mrna.sequence_triplet_binding) if isinstance(val, molecules.Protein)] #liste der Protein-Positionen
        
        pos = ribo_pos+prot_pos # Liste aller R- und Proteinpositionen
    
        np.random.shuffle(pos) # shuffles the list so random order is created
        
        for i in pos:  # iterate through all ribosome/protein positions  
            if self.debug:
                print(i)
            if i == len(mrna.sequence_triplet_binding)-1: # wenn das Ende der mRNA erreicht ist (letztes codon)
                self.entkoppeln(mrna, i)
                self.ribosomes.count += 1   # und erhöht die Menge freier Ribosomen um 1
                continue                                                           
            
            if mrna.sequence_triplet_binding[i+1] == 1:  # besetzt also muss es warten
                if self.debug:
                    print(str(i)+' wir warten')
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

            self.elongate(model, mrna, i)



    def initiate(self, mrna, i):
        """
        Initiate translation, create protein object.

        @type mrna: MRNA
        """
        # !Jens! Da muss jemand kraeftig umbauen
        
        mrna.sequence_triplet_binding[i] = molecules.Protein("Protein_{}".format(mrna.mid),
                                                             "Protein_{0}".format(mrna.name.split("_")[-1]),
                                                             "",)
        if self.debug:
            print('initiate')
            print(mrna.sequence_triplet_binding)

    def elongate(self, model, mrna, i):
        """
        Elongate the new protein by the correct amino acid. Check if an
        MRNA is bound and if ribosome can move to next codon.
        Terminate if the ribosome reaches a STOP codon.

        @type return: Protein or False
        """
        
        print('starting to elongate')
        if isinstance(mrna.sequence_triplet_binding[i], molecules.Protein):
            if self.debug:
                print(str(i)+' elongate')

            codon = mrna[i * 3:i * 3 + 3]
            aa = self.code[codon]

            if aa == "*":  # terminate at stop codon
                self.terminate(mrna, i)

            if i + 1 >= len(mrna.sequence_triplet_binding):
                self.terminate(model, mrna, i)  # terminate if mrna ends

            if mrna.sequence_triplet_binding[i + 1] == 0:  # if the next rna position is free
                if self.debug:
                    print(str(i)+' proteinbau')

                mrna.sequence_triplet_binding[i] += aa       
                self.occupy(mrna, i+1)                      # nächste Stelle besetzten durch ocupyfunction

    def terminate(self, model, mrna, i):
        """
        Splits the ribosome/MRNA complex and returns a protein.
        """
        print('terminator')
        prot = mrna.sequence_triplet_binding[i]  # bound mRNA
        self.entkoppeln(mrna,i)
        self.ribosomes.count += 1
        # put proteins in lists (state)

        if isinstance(prot, molecules.Protein):     # storing the protein in the states-dictionary
            if prot.name in model.states:
                model.states[prot.name].append(prot)
            else:
                model.states[prot.name] = [prot]
        
    

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


    def occupy(self, mrna, i): 
        """
        Movement of Ribosome OR Protein to next codon. 
        !!! Moves from position i-1 to position i TODO:Compatibility to elongate (i-1 --> i)
        Occupies the last 9 codons and frees up the 10th (backwards)
        """

           #Funktion die Besetzung einrichtet
        if self.debug:
             print(str(i)+' occupation neu')
             #print(mrna.sequence_triplet_binding)
       
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
        print(mrna.sequence_triplet_binding)
        if isinstance(mrna.sequence_triplet_binding[i], molecules.Protein):
            print(mrna.sequence_triplet_binding[i].sequence)
        