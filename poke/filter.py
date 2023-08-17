from .pocket import Pocket
from MDAnalysis import AtomGroup

class PocketFilter:
    def __init__(self, pocket: Pocket, atoms: AtomGroup):
        self.pocket = pocket
        self.atoms = atoms

    def __hasitem__(self, key: str):
        return False

    def __getitem__(self, key: str):
        try:
            try:
                return getattr(self, key)
            except AttributeError:
                return getattr(self.pocket, key)
        except:
            raise KeyError(key)
    
    @property
    def resids(self):
        return set(self.atoms[self.pocket.atoms].residues.resids)
    
    def jaccard(self, *residues: list[int]):
        intersection = len(self.resids & set(residues))
        union = len(self.resids) + len(residues) - intersection
        return intersection / union

    def tversky(self, *residues: list[int], a, b):
        intersection = len(self.resids & set(residues))
        diff1 = len(self.resids) - intersection
        diff2 = len(residues) - intersection
        return intersection / (intersection + a*diff1 + b*diff2)

    def evaluate(expr, pocket, atoms):
        return eval(expr, globals(), PocketFilter(pocket, atoms))
    

