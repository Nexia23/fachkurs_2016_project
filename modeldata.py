import pandas
import sqlite3
import os.path
import typecheck as tc
import random as rnd
import molecules as mol

class ModelData:
    """
    class to process data sources for usage in the model
    """
    def __init__(self):
        pass

    def get_mRNA_states(self):
        alphabet = "AUGC"
        # we could improve this by
        # 1) db information about the length
        # 2) random sequence length out of a range 20 - 500
        sequence = ''.join([rnd.choice(alphabet) for i in range(50)])
        mrnas = []
        # we could say that we want to have a random number of a specific mRNA
        for i in range(50):
             mrnas.append(('MRNA_{0}'.format(i), sequence)) # the names are unique we don't need ids
        return mrnas


    def get_states(self, molecule_class):
        """
        retrieves the information required to construct the different model molecules
        @param molecule_class: BioMolecule class
        @return: list
        """
        pass
        #if molecule_class == mol.MRNA:

