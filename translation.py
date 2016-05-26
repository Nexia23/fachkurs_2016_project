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
        # todo: each process class should have something like this
        return ','.join([self.name, str(type(self))])

    def __str__(self):
        # todo: each process class should have something like this
        return ','.join([self.name, str(type(self))])


    def __initiate_ribosomes(self):
        self.ribosomes = molecules.Ribosome('Ribosomes', 'Ribosome', 1)

    def update(self, model):
        """
        Update all mrnas and translate proteins.
        """
<<<<<<< HEAD
        self.ribosomes = model.states[list(self.enzyme_ids)[0]]

        for mrna_id in self.substrate_ids:
            prot = None
            mrna = model.states[mrna_id]

<<<<<<< HEAD
            # ribosomen "arbeiten"
            if not mrna.sequence_triplet_binding[0]:
                self.initiate(mrna)  # TODO:FIXME aufruf self.bind() anstatt initate
            else:  # TODO:FIXME kein else, kann initialisieren UND elongieren in einem Zeitschritt
                prot = self.elongate(mrna)  # TODO:FIXME aufruf self.move() anstatt elongate

            # wenn protein enstanden, in Listen aufnehmen    
=======
=======
        #enzymes_id -> the only one initialized ribosome
        self.ribosomes = model.states[list(self.enzyme_ids)[0]] # call in the dictionary states for the Ribsosome-object
        for mrna_id in self.substrate_ids:          #substrate should be a list for ids of all mrnas
            prot = None
            mrna = model.states[mrna_id]            #object of mRNA
>>>>>>> b175dd75d0aa1382179b9d56668269c2ebc09b61
            if not mrna.sequence_triplet_binding[0]:
                self.initiate(mrna)
            else:
                prot = self.elongate(mrna)
<<<<<<< HEAD
                
>>>>>>> fbeedac7e502267d4fdca1417c91e1fc4793051c
            if isinstance(prot, molecules.Protein):
=======
            if isinstance(prot, molecules.Protein):     #storing the protein in the states-dictionary
>>>>>>> b175dd75d0aa1382179b9d56668269c2ebc09b61
                if prot.name in model.states:
                    model.states[prot.name].append(prot)
                else:
                    model.states[prot.name] = [prot]

    def bind(self):
        # kann self.initiate() oder self.move() aufrufen
        # ribosomenanzahl um 1 reduzieren
        pass

    def move(self):
        # Ribosom läuft auf mRNA ohne protein synthese
        # ruft initiate auf, wenn start codon erreicht
        # ribosom fällt ab, wenn mRNA Ende erreicht
        # ruft self.elongate() auf
        pass

    def initiate(self, mrna):
        """
        Try to bind to a given MRNA. Binding probability corresponds to the ribosome count.

        @type mrna: MRNA
        """
        if not mrna.sequence_triplet_binding[0]:  # no ribosome bound yet and target mrna still free at pos 0
            # bind a nascent protein to the 0 codon
            if numpy.random.poisson(self.ribosomes.count) > 1:  # at least one binding event happens in time step
                mrna.sequence_triplet_binding[0] = molecules.Protein("Protein_{}".format(mrna.mid),
                                                                     "Protein_{0}".format(mrna.name.split("_")[-1]),
                                                                     "",)
                self.ribosomes.count -= 1

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
