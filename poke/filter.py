from functools import cached_property
import numpy as np
from typing import Literal

from .pocket import Pocket
from .utils import Aligner
from MDAnalysis import AtomGroup

class PocketEvaluator:
    def __init__(self, pocket: Pocket, atoms: AtomGroup, reference: AtomGroup):
        self.pocket = pocket
        self.atoms = atoms
        self.reference = reference

    def __hasitem__(self, key: str):
        return False

    def __getitem__(self, key: str):
        try:
            try:
                return getattr(self, key)
            except AttributeError:
                return getattr(self.pocket, key)
        except AttributeError:
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

    # @cached_property
    def rmsd(self, selection: Literal["all", "CA"]="all"):
        if selection == "all":
            atom_ix = self.pocket.atoms

        elif selection == "CA":
            atom_to_res = self.reference.resids 
            atom_to_res -= atom_to_res.min()

            p_resix = np.unique(atom_to_res[self.pocket.atoms])
            atom_ix = self.atoms.residues[p_resix].atoms.select_atoms("name CA").ix

        else:
            raise ValueError(f'Invalid selection for rmsd calculation: {selection}')
            
        aligner = Aligner(self.reference.positions[atom_ix], atom_ix)
        rms = aligner.compute(self.atoms.positions).rmsd
        
        return rms


    def residues(self, *resids: list[int], min_common=None, max_extra=None):
        list_resids = set(resids)
        pock_resids = self.resids

        n_common = len(pock_resids & list_resids)
        n_extra = len(pock_resids) - n_common

        return (
            (min_common is None or n_common >= min_common) and
            (max_extra is None or n_extra <= max_extra)
        )

    def evaluate(expr, pocket, atoms, reference):
        return eval(expr, globals(), PocketEvaluator(pocket, atoms, reference))
    

