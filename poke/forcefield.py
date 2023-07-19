import json
import numpy as np

from pathlib import Path

charmm = Path(__file__).absolute().parent / 'charmm36m'

class AtomsParametrizer:

    def __init__(self, ff: Path) -> None:
        self.residues = json.loads((ff / 'residues.json').read_bytes())
        self.rules = json.loads((ff / 'rules.json').read_bytes())
        self.residues_sets = {resname: [set(i) for i in l] for resname, l in self.residues.items()}
        self.residues_replacers = self.rules['residues_replacers']
        self.atoms_replacers = self.rules['atoms_replacers']
    
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
        resname = self.residues_replacers.get(resname, resname)
        residue_params = self.residues.get(resname)
        residue_sets = self.residues_sets.get(resname)

        if residue_params is None:
            raise RuntimeError(f'Residue {resname} unknown')

        # Fetch the replacer for this residue, use the default if there is no specific replacer
        replacer = self.atoms_replacers.get(resname, self.atoms_replacers['default'])

        # Replace the atom names of the residue
        replaced_names = [n[1:]+n[0] if n[0].isdigit() else n for n in atom_names]
        replaced_names = [replacer.get(n, n) for n in replaced_names]
        replaced_names_set = set(replaced_names)
       
        # Find the parameters for this residue using the names of the atoms
        try:
            best_params = residue_params[residue_sets.index(replaced_names_set)]
        except:
            error_msg = f'Some atoms names of residue {resname} are not matching the forcefield:'
     
            for i, trial in enumerate(residue_sets):
                error_msg += (f'\n\nForcefield option {i+1}: {", ".join(sorted(trial))}' +
                             '\nUnmatching input atoms:      ' + ", ".join(sorted(replaced_names_set - trial)) +
                             '\nUnmatching forcefield atoms: ' + ", ".join(sorted(trial - replaced_names_set)))
                                
            raise ValueError(error_msg)
  
        return np.array([best_params[atom] for atom in replaced_names], dtype=np.float32)
        


    def match_score(names, forcefield):
        '''
        Jaccard similarity between the two sets of atoms names.
        '''
        intersection = len(names & forcefield) 
        return intersection / (len(names) + len(forcefield) - intersection)

    def find_best_match(names, forcefield_sets, forcefield_params):
        best_params = None
        best_score = -1

        for trial, params in zip(forcefield_sets, forcefield_params):
            in_names = names - trial
            ff_names = trial - names

            


            # elif score > best_score:
            #     best_params = params
            #     best_score = score
        return best_params