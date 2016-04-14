import unittest
from unittest.mock import patch

import modeldata
import molecules
import translation
import model

class TestModel(unittest.TestCase):
    def setUp(self):
        self.m = model.Model()

    def test_timestep(self):
        '''
        calling step should increase the time
        @return:
        '''
        start = self.m.timestep
        self.m.step()
        self.assertGreater(self.m.timestep, start)

    def test_simulate_time(self):
        '''
        simulate should increase timestep by x
        @return:
        '''
        x = 100
        start = self.m.timestep
        self.m.simulate(x, log=False)
        self.assertEqual(start + x, self.m.timestep)


class TestData(unittest.TestCase):
    def test_data_mrna(self):
        '''
        calling get states with the mol.MRNA class should return ids with mRNA in them.
        @return:
        '''
        db = modeldata.ModelData()
        #mrnas = db.get_states(molecules.MRNA)
        mrnas = db.get_mRNA_states()
        for mrna in mrnas:
            self.assertRegex(mrna[0], "MRNA_\d+")

class TestTranslation(unittest.TestCase):
    def setUp(self):
        self.m = model.Model()
        self.t = translation.Translation("trsl", "test_translation")

    @patch('translation.numpy.random.poisson')
    def test_initiation(self, npr_mock):
        npr_mock.return_value = 2
        mrna = molecules.MRNA("mrna", "AUAUAUAUAAUG")
        self.t.initiate(mrna)
        self.assertIsInstance(mrna.sequence_triplet_binding[0], molecules.Protein)
        npr_mock.assert_called_with(1)



class TestSomething(unittest.TestCase):

    def test_model(self):
        m = model.Model()
        m.simulate(100)

if __name__ == '__main__':
    unittest.main()