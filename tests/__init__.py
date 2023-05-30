import unittest
import numpy as np

from pathlib import Path
from MDAnalysis import Universe
from poke import Poke, Pocket

data_path = Path(__file__).parent / 'data'

class Tests(unittest.TestCase):

    def test_1gkc(self):
        u = Universe(data_path / '1gkc.gro')

        scouter = Poke()
        scouter.load_atoms_properties_from_mdanalysis(u.atoms)
        pockets = scouter.get_pockets(u.atoms.positions)  

        for p in pockets:
            self.assertIsInstance(p, Pocket)
            self.assertFinite(p)

        self.assertGreaterEqual(len(pockets), 1)
        self.assertLess(len(pockets), 3)
        

    @staticmethod
    def assertFinite(obj):
        for key, value in obj.__dict__.items():
            if isinstance(value, (float, np.float32, np.float64, np.ndarray)):
                if np.any(np.isnan(value)):
                    raise AssertionError(f'The member "{key}" of {obj.__class__.__name__} has some NaN value.')
                if np.any(np.isinf(value)):
                    raise AssertionError(f'The member "{key}" of {obj.__class__.__name__} has some inf value.')



if __name__ == '__main__':
    unittest.main()
