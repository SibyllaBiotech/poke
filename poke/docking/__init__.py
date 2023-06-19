import os
import warnings
import numpy as np
from pathlib import Path
from MDAnalysis import Universe

from .gnina import Gnina
from .molecule import Molecule


def dock_ligands(receptor: bytes, box_center: list[float], box_size: list[float], smiles: list[str], exhaustiveness=24, n_poses=5, n_cpus=1):

    docker = Gnina(protein=receptor, box_center=box_center, box_size=box_size, 
                   n_cpus=n_cpus, exhaustiveness=exhaustiveness, 
                   n_poses=n_poses, verbose=False)

    molecules = [Molecule.from_smiles(s) for s in smiles]
    return molecules

def write_box(path: str | Path, center: np.ndarray, size: np.ndarray):
    box = Universe.empty(n_atoms=9, n_residues=1, n_segments=1, atom_resindex=[0]*9, 
                         residue_segindex=[0], trajectory=True)
    box.add_TopologyAttr('name', ['O']+['C']*8)
    box.add_TopologyAttr('resname', ['BXA'])
    box.add_TopologyAttr('resid', [1])
    box.add_TopologyAttr('segid', ['B'])
    box.add_TopologyAttr('bonds', [(1,2), (3,4), (5,6), (7,8),
                                   (1,3), (2,4), (5,7), (6,8),
                                   (1,5), (2,6), (3,7), (4,8)])

    box.atoms[0].position = center
    box.atoms[1:].positions = center - size/2 + size * np.unpackbits(np.arange(8, dtype=np.uint8)).reshape(8,-1)[:,5:]
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=Warning)
        box.atoms.write(path)

def write_box_vectors(path: str | Path, origin: np.ndarray, vectors: np.ndarray):
    box = Universe.empty(n_atoms=8, n_residues=1, n_segments=1, atom_resindex=[0]*8, 
                         residue_segindex=[0], trajectory=True)
    box.add_TopologyAttr('name', ['C']*8)
    box.add_TopologyAttr('resname', ['BXA'])
    box.add_TopologyAttr('resid', [1])
    box.add_TopologyAttr('segid', ['B'])
    box.add_TopologyAttr('bonds', [(0,1), (2,3), (4,5), (6,7),
                                   (0,2), (1,3), (4,6), (5,7),
                                   (0,4), (1,5), (2,6), (3,7)])
    
    box.atoms.positions = origin + np.unpackbits(np.arange(8, dtype=np.uint8)).reshape(8,-1)[:,5:] @ vectors
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=Warning)
        box.atoms.write(path)







    
    


    
    



