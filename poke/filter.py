from functools import cached_property

from .pocket import Pocket
from MDAnalysis import AtomGroup

class PocketEvaluator:
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

    @cached_property
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


    def residues(self, *resids: list[int], min_common=None, max_extra=None):
        list_resids = set(resids)
        pock_resids = self.resids

        n_common = len(pock_resids & list_resids)
        n_extra = len(pock_resids) - n_common

        return (
            (min_common is None or n_common >= min_common) and
            (max_extra is None or n_extra <= max_extra)
        )

    def evaluate(expr, pocket, atoms):
        return eval(expr, globals(), PocketEvaluator(pocket, atoms))
    

