import processes
import molecules
import numpy


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
        for mrna_id in self.substrate_ids:          # substrate should be a list for ids of all mrnas
            prot = None                             # initialise prot variable
            mrna = model.states[mrna_id]            # object of mRNA

            # ribosoms work: bind, move, initialise, elongate
            if mrna.sequence_triplet_binding[0] == 0:
                self.bind(mrna)  # bind to first position on mRNA
                self.move(mrna)  # move ribosome that has not yet initiated translation
                prot = self.elongate(mrna)  # TODO:FIXME aufruf self.move() anstatt elongate
                
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
        if numpy.random.poisson(self.ribosomes.count) > 1:  # check if ribosome wants to bind
            mrna.sequence_triplet_binding[0] = 1            # check if first codon on mRNA is free
            if mrna.sequence[0:3]=='AUG':                   # if first codon is START codon -> initiate
                self.initiate(mrna)

            self.ribosomes.count -= 1  # remove bound ribosome from list of free ribosomes

    def move(self, mrna):
        """
        move not initiated ribosomes
        """
        # for i, codon in enumerate(mrna.sequence_triplet_binding):
            # if codon == 'R':

        # Ribosom läuft auf mRNA ohne protein synthese
        # # 1. Fall: nächstes codon belegt (mrna.sequence_triplet_binding[0] == 1)
        # # 2. Fall: nächstes codon frei und kein Start codon -> eins weiter bewegen
        # # 3. Fall: nächstes codon frei und Start codon -> initiate() und elongate() aufrufen und beenden
        # # 4. Fall: kein nächstes codon -> abfallen (mrna.sequence_triplet_binding[0] = 0 und self.ribosome.count += 1)
        # sobald ein Ribosom Protein ist, elongate() aufrufen und beenden

            # elif isinstance(ribosome, molecules.Protein):

        pass

    def initiate(self, mrna):
        """
        Initiate translation, create protein object.

        @type mrna: MRNA
        """

        mrna.sequence_triplet_binding[0] = molecules.Protein("Protein_{}".format(mrna.mid),
                                                             "Protein_{0}".format(mrna.name.split("_")[-1]),
                                                             "",)
 

    def elongate(self, mrna):
        """
        Elongate the new protein by the correct amino acid. Check if an
        MRNA is bound and if ribosome can move to next codon.
        Terminate if the ribosome reaches a STOP codon.

        @type return: Protein or False
        """
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
                    mrna.sequence_triplet_binding[i + 1] = mrna.sequence_triplet_binding[i]
                    mrna.sequence_triplet_binding[i] = 0
        return 0
    
    def terminate(self, mrna, i):
        """
        Splits the ribosome/MRNA complex and returns a protein.
        """
        protein = mrna.sequence_triplet_binding[i]  # bound mRNA
        mrna.sequence_triplet_binding[i] = 0
        self.ribosomes.count += 1
        return protein
    

    def entkoppeln(self, mrna, i):      #Funktion die Entkoppelt
        if i < 10:
            k=0
            while k <= i:
                mrna.sequence_triplet_binding[k] = 0
                k += 1  # k <= i ist sonst unendlich lange richtig!!

        elif i >=10:
            mrna.sequence_triplet_binding[i-10]=0  # das ist nur eien Stelle!


    def occupy(self, mrna, i):        #Funktion die Besetzung einrichtet

        if i < 10:
            k=0
            while k <= i:
                mrna.sequence_triplet_binding[k] = 1
                k += 1  # k <= i ist sonst unendlich lange richtig!!

        elif i >= 10:
            mrna.sequence_triplet_binding[i-10] = 'R'  # das ist nur eien Stelle!

            

