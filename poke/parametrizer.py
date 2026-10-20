import json
import numpy as np

from functools import cache
from pathlib import Path

charmm = Path(__file__).absolute().parent / 'charmm36m'

class AtomsParametrizer:
    names_replacer = {
        "OC1": "OT1",
        "OC2": "OT2"
    }

    def __init__(self, ff: Path) -> None:
        self.residues, self.resnames, self.atoms_types = cached_loader(ff)

    def __call__(self, resname, atom_names):
        '''
        Compute parameters for a single residue. 
        Uses data loaded in the constructor.

        Parameters
        ------------
            resname: str
                the name of the residue
            atom_names: list[str]
                the list of atoms names of the residue.
                The order of the output parameters will be 
                in the same order of this list

        Return values
        ------------
            charges:
                charges of the atoms in electron unit
            radii:
                Van der Waals radii of the atoms
            epsilons:
                depth of the Van der Waals interaction    
        '''
        residues = self.resnames.get(resname)

        if residues is None:
            raise RuntimeError(f'Resname {resname} unknown. Add it to poke/charmm36m/resnames.json')
        
        input_atom_names = [n[1:]+n[0] if n[0].isdigit() else n for n in atom_names]
        input_atom_names = [AtomsParametrizer.names_replacer.get(n,n) for n in input_atom_names]
        input_atoms = set(input_atom_names)
        
       
        # Find the parameters for this residue using the names of the atoms
        best_types = None
        best_len = 0

        for residue_key in residues:
            residue = self.residues[residue_key]
            params = residue['params']

            for ff_atoms_list in residue["names"]:
                ff_atoms = set(ff_atoms_list)
            
                if ff_atoms >= input_atoms and (best_types is None or len(ff_atoms) < best_len):
                    best_types = [params[ff_atoms_list.index(name)] for name in input_atom_names]  
                    best_len = len(ff_atoms)
                    
                    if best_len == len(input_atoms):
                        break

        if best_types is not None:
            return self.atoms_types[best_types]

        else:
            error_msg = f'Some atoms names of residue {resname} with atoms {", ".join(sorted(input_atoms))} are not matching the forcefield:'
     
            for residue_key in residues:
                residue = self.residues[residue_key]

                for i, ff_atoms_list in enumerate(residue["names"]):
                    ff_atoms = set(ff_atoms_list)

                    error_msg += (f'\n\nForcefield option {residue_key} {i}: {", ".join(sorted(ff_atoms))}' +
                                '\nUnmatched input atoms:      ' + ", ".join(sorted(input_atoms - ff_atoms)) +
                                '\nUnmatched forcefield atoms: ' + ", ".join(sorted(ff_atoms - input_atoms)))
                                    
            raise ValueError(error_msg)
        
@cache
def cached_loader(ff: Path):
    residues = json.loads((ff / 'residues.json').read_bytes())
    resnames = json.loads((ff / 'resnames.json').read_bytes())
    atoms_types = np.loadtxt(ff / 'atoms.csv', dtype=np.float32, delimiter=',')
    return residues, resnames, atoms_types