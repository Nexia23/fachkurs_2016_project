import random as rnd
import string

import translation
import molecules as mol


class ModelData:
    """
    class to process data sources for usage in the model
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

    def __init__(self):
        pass

    def get_states(self, molecule_class):
        """
        retrieves the information required to construct the different model molecules
        @param molecule_class: BioMolecule class
        @return: list
        """

        if molecule_class == mol.MRNA:
            alphabet = list(self.code.keys())
            mrnas = []
            genes = {}
            for i in range(10):
                sequence = ''.join([rnd.choice(alphabet) for i in range(rnd.randint(50, 500))])
                genes[''.join([rnd.choice(string.ascii_uppercase) for i in range(3)])] = sequence


            for gene in genes:
                for i in range(rnd.randint(1, 10)):
                    mrnas.append(("MRNA_{}_{}".format(gene, i), gene, genes[gene]))
            return mrnas
