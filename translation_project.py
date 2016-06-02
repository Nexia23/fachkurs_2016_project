import random
import sys

class BioMolecule(object):
    """
    A generic molecule that has basic attributes like id, name and
    mass.

    @type id: int
    @type name: str
    @type mass: float
    """
    def __init__(self, id, name, mass=0.):
        self.id = id
        self.name = name
        self.mass = mass
    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @property
    def mass(self):
        return self._mass

    @id.setter
    def id(self, value) :
        if not isinstance(value, int):
            raise TypeError("Mass must be Integer.")
        self._id = value

    @name.setter
    def name(self, value) :
        if not isinstance(value, str):
            raise TypeError("Mass must be String.")
        self._name = value

    @mass.setter
    def mass(self, value):
        if not isinstance(value, float):
            print(value)
            raise TypeError("Mass must be Float.")
        self._mass = value


    
    
    # 1. Write setter and getter methods for all attributes.
    #    Use @property decorators as dicussed in the lecture
    # 2. In the setter methods check for the type of each attribute.

class Polymer(BioMolecule):
    """
    A polymer molecule that has a sequence attribute which is
    accessible via indexing the object. 

    @type id: int
    @type name: str
    @type sequence: str
    @type mass: float
    """
    def __init__(self, id, name, sequence, mass=0.):
        # 3. Initialize the parent class correctly
        super().__init__( id, name, mass)
        self._sequence = sequence

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, value) :
        if not isinstance(value, str):
            raise TypeError("Mass must be String.")
        self._sequence = value

    
    # 4. Write getter and setter for sequence, again check for type
    # 5. run in ipython, instantiate this class, and test it
    def __getitem__(self, value):
        """
        Makes the sequence accessible via the indexing operators:
<        p[10] returns the tenth character in the sequence.
        """
        return self.sequence[value]

    def __setitem__(self, key, value):
        """
         Enables changing of sequence characters via the indexing operators.       
        """
        self.sequence[key] = value


class MRNA(Polymer):
    def __init__(self, id, name, sequence, mass=0.):
        # 6. Initialize the parent class correctly
        super().__init__( id, name, sequence, mass)
     
        # 7. Create a list that stores if a ribosome is bound for each
        # codon (triplet).
        bond_rib= len(sequence)//3
        self.binding = [0]*bond_rib # use this attribute for 7.
        self.mass = self.calculate_mass()

    def calculate_mass(self):
        NA_mass = {'A': 1.0, 'U': 2.2, 'G':2.1, 'C':1.3}
        self.mass = 0.

        for base in self.sequence:
            self.mass += NA_mass[base] 

        return self.mass    


        # 8. calculate the mass for the whole sequence

class Protein(Polymer):
    """Protein with Polymer features and mass calculation. A global class
    attribute counts the number of proteins that have been instantiated.
    
    A protein can be elongated using the "+" operator:
    
    >> protein = Protein(1, "prot", "MVFT")
    >> protein + "A"
    >> protein.sequence
    MVFTA

    
    
    """
    number_of_proteins = 0  # init instance counter

    def __init__(self, id, name, sequence, mass=0.):
        super().__init__(id, name, sequence, mass)
        self.__class__.number_of_proteins += 1 #  increase instance counter
        self.calculate_mass()
        

    def __add__(self, other):
        self.sequence = self.sequence + other
        return self
    # 9. implement the elongation feature described in the docstring. (__add__)

    def calculate_mass(self):
        AA_mass = {"A": 89.0929, "R": 175.208, "N": 132.118, "D": 132.094, "C": 121.158, "Q": 146.144,
                    "E": 146.121, "G": 75.0664, "H":155.154, "I":131.172, "L": 131.172, "K": 147.195,
                    "M": 149.211, "F": 165.189, "P": 115.13, "S": 105.092, "T": 119.119, "W": 204.225,
                    "Y":181.188, "V":117.146}
        self.mass = 0.
        for aa in self.sequence:
            self.mass += AA_mass[aa]
        
   

class Ribosome(BioMolecule):
    """A ribosome can bind MRNA and translate it. After translation is
    finished it produces a protein.

    During initiation the ribosome checks if a given MRNA is bound
    by another ribosome and binds only if position 0 is empty.

    Elongation checks if the next codon is unbound and only elongates
    if the ribosome can move on. If the ribosome encounters a stop
    codon ("*") translation terminates. The MRNA is detached from the
    ribosome and the finished protein is returned.
    """
    code = dict([('UCA','S'), ('UCG','S'), ('UCC','S'), ('UCU','S'),
                 ('UUU','F'), ('UUC','F'), ('UUA','L'), ('UUG','L'),
                 ('UAU','Y'), ('UAC','Y'), ('UAA','*'), ('UAG','*'),
                 ('UGU','C'), ('UGC','C'), ('UGA','*'), ('UGG','W'),
                 ('CUA','L'), ('CUG','L'), ('CUC','L'), ('CUU','L'),
                 ('CCA','P'), ('CCG','P'), ('CCC','P'), ('CCU','P'),
                 ('CAU','H'), ('CAC','H'), ('CAA','Q'), ('CAG','Q'),
                 ('CGA','R'), ('CGG','R'), ('CGC','R'), ('CGU','R'),
                 ('AUU','I'), ('AUC','I'), ('AUA','I'), ('AUG','M'),
                 ('ACA','T'), ('ACG','T'), ('ACC','T'), ('ACU','T'),
                 ('AAU','N'), ('AAC','N'), ('AAA','K'), ('AAG','K'),
                 ('AGU','S'), ('AGC','S'), ('AGA','R'), ('AGG','R'),
                 ('GUA','V'), ('GUG','V'), ('GUC','V'), ('GUU','V'),
                 ('GCA','A'), ('GCG','A'), ('GCC','A'), ('GCU','A'),
                 ('GAU','D'), ('GAC','D'), ('GAA','E'), ('GAG','E'),
                 ('GGA','G'), ('GGG','G'), ('GGC','G'), ('GGU','G')])

    def __init__(self, id, name):
        super().__init__(id, name)
        self.bound_mrna = False
        self.position = None  # position on a bound MRNA

    def initiate(self, mrna):
        """
        Tries to bind to a given MRNA.

        @type mrna: MRNA
        """
        if not self.bound_mrna and mrna.binding[0]==0:  # no mrna bound
                                                    # yet and target
                                                    # mrna still free
                                                    # at pos 0
            self.bound_mrna = mrna
            self.nascent_prot = Protein(mrna.id, 'Product','')  # 10. Instantiate a new Protein
            self.position = 0
            self.bound_mrna.binding[0]=1  # 11. Mark position 0 of MRNA to be bound by ribosome
            
    def elongate(self):
        """Elongate the new protein by the correct amino acid. Check if an
        MRNA is bound and if ribosome can move to next codon.
        Terminate if the ribosome reaches a STOP codon.

        @type return: Protein or False
        """
        if not self.bound_mrna: # can't elongate because there is no MRNA
            return False

        if self.bound_mrna.binding[1]==1:
            return False

        if self.bound_mrna.binding[0]==1:
            
            self.nascent_prot + self.code[self.bound_mrna.sequence[0:3]]
            i=1

            while self.bound_mrna.binding[i]!= 1:

                self.bound_mrna.binding[i-1] = 0
                self.bound_mrna.binding[i] = 1
                seq_i = i*3
                seq_o = seq_i+3
                self.nascent_prot + self.code[self.bound_mrna.sequence[seq_i:seq_o]]
                
                if self.code[self.bound_mrna.sequence[seq_i:seq_o]]== '*':
                    return self.terminate()
                

                i = i+1

        return 0
        # 12. Implement the described features.

    def terminate(self):
        """
        Splits the ribosome/MRNA complex and returns a protein.

        """
        if self.nascent_prot.sequence[-1]=='*':
            self.bound_mrna.binding[-1] = 0
        # 13. Dissociate the complex.
        return self.nascent_prot
        

class Cell(object):
    def __init__(self):
        self.ribosomes = [Ribosome(i, 'Ribo_{0}'.format(i)) for i in range(200)]
        self.mrnas = [MRNA(i, 'MRNA_{0}'.format(i), "UUUUUUUUUUAA") for i in range(20)]
        self.proteins = [[] for x in range(20)]

    def step(self):
        for r in self.ribosomes:
            if not r.bound_mrna:
                r.initiate(self.mrnas[random.randint(0,len(self.mrnas)-1)])
            else:
                prot = r.elongate()
                if prot:
                    self.proteins[prot.id].append(prot)

    def simulate(self, steps, p=True):
        for s in range(steps):
            self.step()
            if p:
                print([len(x) for x in self.proteins])

            
if __name__ == "__main__":  # the following is called if the module is executed
    # 14. Instantiate the Cell class and call the simulation method.
    test = Cell()
    test.simulate(10)
    
    
# 15. Generate a set of mRNA sequences to initiate the cell.
# 16. Implement protein degradation.
