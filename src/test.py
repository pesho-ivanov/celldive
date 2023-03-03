import unittest
import pandas as pd

from pandas.util.testing import assert_frame_equal

from structs.tcrs import TCRs
import utils

class TestStringMethods(unittest.TestCase):

    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

    def test_isupper(self):
        self.assertTrue('FOO'.isupper())
        self.assertFalse('Foo'.isupper())

    def test_split(self):
        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)
            
class TestUtils(unittest.TestCase):
    def test_refresh_index(self):
        tested_df = pd.DataFrame({'a': [1,2,3], 'b': [5,8,10]}, index=['08','2','14'])
        expected_df = pd.DataFrame({'a': [2,1,3], 'b': [8,5,10]}, index=['2','8','14'])
        pesho_utils.refresh_index(tested_df)
        assert_frame_equal(tested_df, expected_df)
        
class TestTCRs(unittest.TestCase):
    def test_cells(self):
        arrays = [ [10,11,12,13], ['A','B','A','B'], ['tcr1', 'tcr2', 'tcr3', 'tcr4'] ]
        tuples = list(zip(*arrays))
        index = pd.MultiIndex.from_tuples(tuples, names=['cell', 'chain', 'seq'])
        tcrs = TCRs()
        cells = tcrs.cells()
        print(cells)
        #tcrs = pd.DataFrame({}: )
        
    #def test_getEdges(self):
        
testedClasses = [ TestStringMethods,
                  #TestUtils,
                  #TestTCRs,
                ]

#unittest.main()
#unittest.main(argv=['first-arg-is-ignored'], exit=False)
def test_all():
    #suite = unittest.TestLoader().loadTestsFromTestCase(TestStringMethods)
    suite = unittest.TestSuite()
    for cls in testedClasses:
        suite.addTest(unittest.makeSuite(cls))
    #suite = unittest.TestLoader().loadTestsFromModule(TestStringMethods)
    unittest.TextTestRunner(verbosity=2).run(suite)