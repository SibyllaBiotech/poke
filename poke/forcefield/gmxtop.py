import io
from pathlib import Path

class Tokenizer:
    def __init__(self, stream: io.TextIOBase) -> None:
        self.stream = stream
        self.fetch_tokens()

    def fetch_tokens(self):
        while True:
            line = self.stream.readline()
            self.tokens = line.split(';')[0].split('#')[0].split()

            if self.tokens:
                return True

            elif not line:
                self.tokens = None
                return False
        
def parse_ffnonbonded(ffpath: Path):
    atom_types = {}
    with open(ffpath / 'ffnonbonded.itp') as file:
        s = Tokenizer(file)

        if s.tokens == ['[', 'atomtypes', ']']:
            s.fetch_tokens()

            while len(s.tokens) == 7:
                atom_types[s.tokens[0]] = (round(10 * float(s.tokens[5]), 3), round(float(s.tokens[6]), 3))
                s.fetch_tokens()
    return atom_types

def parse_rtp(rtp_file: Path, atoms_types: dict[str, tuple], resnames: dict[str,list[str]]):
    residue_types = {}
    residues_map = {}

    sections = [
        'bonds',
        'impropers',
        'cmap'
    ]

    with open(rtp_file) as file:
        s = Tokenizer(file)

        if s.tokens == ['[', 'bondedtypes', ']']:
            s.fetch_tokens()
            s.fetch_tokens()

        while True:        
            while not (len(s.tokens) == 3 and s.tokens[0] == '[' and s.tokens[1] not in sections and s.tokens[2] == ']'):
                s.fetch_tokens()

                if s.tokens is None:
                    return residue_types

            resname = s.tokens[1]
            atoms = {}
            bonds = []
            
            s.fetch_tokens()
            if s.tokens == ['[', 'atoms', ']']:
                s.fetch_tokens()
               
                while len(s.tokens) == 4:
                    name, ffname, charge, _ = s.tokens
                    atoms[name] = (round(float(charge), 3), *atoms_types[ffname])
                    s.fetch_tokens()

            if s.tokens == ['[', 'bonds', ']']:
                s.fetch_tokens()
               
                while len(s.tokens) == 2:
                    atom1, atom2 = s.tokens
                    if not ('+' in atom1 or '-' in atom1 or '+' in atom2 or '-' in atom2):
                        bonds.append((atom1.strip('+-'), atom2.strip('+-')))
                    s.fetch_tokens()

            key = frozenset((tuple(atoms.items()), tuple(bonds)))

            if key in residues_map:
                other_resname = residues_map[key]

                if len(other_resname) > len(resname):
                    del residue_types[other_resname]
                    residues_map[key] = resname
                    residue_types[resname] = {"atoms": atoms, "bonds": bonds}
                    resnames[resname] = [resname]

            else:
                residues_map[key] = resname
                residue_types[resname] = {"atoms": atoms, "bonds": bonds}
                resnames[resname] = [resname]

def parse_tdb(tdb_file: Path, residue_types: dict[str,dict[str,tuple]], atoms_types: dict[str, tuple], resnames: dict[str, list[str]]):
    sections = [
        'delete',
        'replace',
        'add',
        'impropers',
        'propers'
    ]

    new_residue_types = {}
    
    with open(tdb_file) as file:
        s = Tokenizer(file)

        while True:        
            if len(s.tokens) == 3 and s.tokens[0] == '[' and s.tokens[2] == ']':

                if s.tokens[1] == 'None':
                    s.fetch_tokens()
                    continue
                
                mod_name = s.tokens[1]
                mod_delete = set()
                mod_replace = set()
                mod_add = []

                s.fetch_tokens()

                while len(s.tokens) == 3 and s.tokens[0] == '[' and s.tokens[1] in sections and s.tokens[2] == ']':
                    if s.tokens[1] == 'delete':
                        s.fetch_tokens()

                        while len(s.tokens) == 1:
                            mod_delete.add(s.tokens[0])
                            s.fetch_tokens()

                    elif s.tokens[1] == 'replace':
                        s.fetch_tokens()

                        while len(s.tokens) == 4:
                            mod_replace.add((s.tokens[0], round(float(s.tokens[3]), 3), *atoms_types[s.tokens[1]]))
                            s.fetch_tokens()

                    elif s.tokens[1] == 'add':
                        s.fetch_tokens()

                        while s.tokens is not None and s.tokens[0] != '[':
                            atom_name = s.tokens[2]
                            n_atoms = int(s.tokens[0])
                            relative_to = s.tokens[3:]
                            s.fetch_tokens()
                            mod_add.append((n_atoms, atom_name, tuple(relative_to), round(float(s.tokens[2]),3), *atoms_types[s.tokens[0]]))
                            s.fetch_tokens()

                    else:
                        while s.tokens is not None and not (len(s.tokens) == 3 and s.tokens[0] == '[' and not s.tokens[1] in sections and s.tokens[2] == ']'):
                            s.fetch_tokens()

                    if s.tokens is None:
                        break

                mod_name_toks = mod_name.split('-')

                if any(len(t) == 0 for t in mod_name_toks):
                    mod_name_toks = [mod_name]

                for resname, data in residue_types.items():
                    if len(mod_name_toks) == 2:
                        if mod_name_toks[0] not in resname:
                            continue
                        mod_name = mod_name_toks[1]

                    mod_ok = True
                    atoms = data['atoms'].copy()
                    bonds = data['bonds'].copy()

                    for atom in mod_delete:
                        if atom in atoms:
                            del atoms[atom]
                            bonds = [bond for bond in bonds if atom not in bond]
              
                    for atom, charge, sigma, epsilon in mod_replace:
                        if atom in atoms:
                            atoms[atom] = (charge, sigma, epsilon)
                        else:
                            mod_ok = False

                    for n_atoms, atom, relative_to, charge, sigma, epsilon in mod_add:
                        for i in range(n_atoms):
                            atom_name = f'{atom}{i+1}' if n_atoms > 1 else atom

                            if atom_name in atoms:
                                mod_ok = False

                            atoms[atom_name] = (charge, sigma, epsilon)

                            if any(n not in atoms for n in relative_to):
                                mod_ok = False

                            bonds.append((atom_name, relative_to[0]))

                    
                    if mod_ok and (f'{resname}-{mod_name}' not in new_residue_types or len(mod_name_toks) > 1):
                        new_residue_types[f'{resname}-{mod_name}'] = {'atoms': atoms, 'bonds': bonds}    
                        resnames[resname].append(f'{resname}-{mod_name}')        

            if s.tokens is None:
                break
    
    return new_residue_types