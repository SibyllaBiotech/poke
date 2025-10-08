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

def parse_rtp(rtp_file: Path, atoms_types: dict[str, tuple], residues_map: dict):

    sections = [
        'bonds',
        'impropers',
        'cmap',
        'angles'
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
                    return

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
                    if s.tokens is None:
                        break

            if s.tokens == ['[', 'bonds', ']']:
                s.fetch_tokens()
               
                while len(s.tokens) == 2:
                    atom1, atom2 = s.tokens
                    if not ('+' in atom1 or '-' in atom1 or '+' in atom2 or '-' in atom2):
                        bonds.append((atom1.strip('+-'), atom2.strip('+-')))
                    s.fetch_tokens()

            key = frozenset((tuple(atoms.items()), tuple(bonds)))
            residues_map.setdefault(key, [{"atoms": atoms, "bonds": bonds, "resnames": []}])[0]['resnames'].append(resname)

            if s.tokens is None:
                break


def parse_tdb(tdb_file: Path, residues_map: dict, atoms_types: dict[str, tuple], mods_spec: dict[str, dict[str, dict]]):
    
    sections = [
        'delete',
        'replace',
        'add',
        'impropers',
        'propers'
    ]

    mods_map = {}

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
                            mod_function = int(s.tokens[1])
                            n_atoms = int(s.tokens[0])
                            relative_to = s.tokens[3:]
                            s.fetch_tokens()
                            mod_add.append((mod_function, n_atoms, atom_name, tuple(relative_to), round(float(s.tokens[2]),3), *atoms_types[s.tokens[0]]))
                            s.fetch_tokens()

                    else:
                        while s.tokens is not None and not (len(s.tokens) == 3 and s.tokens[0] == '[' and not s.tokens[1] in sections and s.tokens[2] == ']'):
                            s.fetch_tokens()

                    if s.tokens is None:
                        break
            
                mod_add_condensed = []
                i = 0
                while i < len(mod_add):
                    mod_function, n_atoms, atom_name, relative_to, _, _, _ = mod_add[i]
                    new_names = [(f'{atom_name}{i+1}' if n_atoms > 1 else atom_name) for i in range(n_atoms)]
                    if mod_function == 9 and i + 1 < len(mod_add):
                        # Check if we can merge with the next if it is type 2
                        next_function, _, next_atom_name, next_relative_to, _, _, _ = mod_add[i + 1]
                        if next_function == 2 and next_relative_to[0] == new_names[1] and next_relative_to[1] == relative_to[0] and next_relative_to[2] == relative_to[1]:
                            new_names.append(next_atom_name)
                            i += 1
                    mod_add_condensed.append({
                        'type': mod_function,
                        'names': new_names,
                        'anchors': relative_to,
                        'resindex': [0] * len(relative_to)  # Assuming all are relative to the same residue
                    })
                    i += 1
                
                mod_name_toks = mod_name.split('-')
                spec_res = 'default'
                if any(len(t) == 0 for t in mod_name_toks):
                    mod_name_toks = [mod_name]
                
                if len(mod_name_toks) == 2:
                    mod_name = mod_name_toks[1]
                    spec_res = mod_name_toks[0]

                mods_spec.setdefault(mod_name, {})[spec_res] = {
                    'delete': list(mod_delete),
                    'add': mod_add_condensed,
                }

                for res_key, data_list in residues_map.items():
                    new_data_list = []
                    for data in data_list:
                    
                        if len(mod_name_toks) == 2:
                            if all(mod_name_toks[0] not in resname for resname in data['resnames']):
                                continue
                            
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

                        for _, n_atoms, atom, relative_to, charge, sigma, epsilon in mod_add:
                            for i in range(n_atoms):
                                atom_name = f'{atom}{i+1}' if n_atoms > 1 else atom

                                if atom_name in atoms:
                                    mod_ok = False

                                atoms[atom_name] = (charge, sigma, epsilon)

                                if any(n not in atoms for n in relative_to):
                                    mod_ok = False

                                bonds.append((atom_name, relative_to[0]))

                        
                        if mod_ok:
                            res_mod_name = min(data['resnames'], key=lambda x: len(x)) + '-' + mod_name
                            new_data_list.append({
                                'atoms': atoms,
                                'bonds': bonds,
                                'resnames': [res_mod_name]
                            })

                    mods_map.setdefault(res_key, []).extend(new_data_list) 

            if s.tokens is None:
                for res_key, data_list in mods_map.items():
                    residues_map[res_key].extend(data_list)
                break


def parse_hdb(hdb_file: Path, hdb_map: dict):
    
    with open(hdb_file) as file:
        s = Tokenizer(file)
        while True:   
            if len(s.tokens) == 2:
                resname, n_lines = s.tokens
                n_lines = int(n_lines)
                mods = []
                for _ in range(n_lines):
                    s.fetch_tokens()
                    if s.tokens is None:
                        break
                    n_h = int(s.tokens[0])
                    tp = int(s.tokens[1])
                    hname = s.tokens[2]
                    anchors_raw = s.tokens[3:]
                    anchors = []
                    relixs  = []
                    for ar in anchors_raw:
                        if ar.startswith('-'):
                            relixs.append(-1)
                            anchors.append(ar[1:])
                        else:
                            relixs.append(0)
                            anchors.append(ar)
                    
                    new_names = [(f'{hname}{i+1}' if n_h > 1 else hname) for i in range(n_h)]
                    
                    mods.append((tp,new_names,anchors,relixs))
                else:
                    s.fetch_tokens()

                mods_condensed = []
                i = 0
                while i < len(mods):
                    tp, new_names, anchors, relixs = mods[i]
                    if tp == 9 and i + 1 < len(mods):
                        # Check if we can merge with the next if it is type 2
                        next_tp, next_new_names, next_anchors, next_relixs = mods[i+1]
                        if next_tp == 2 and next_anchors[0] == new_names[1] and next_anchors[1] == anchors[0] and next_anchors[2] == anchors[1]:
                            new_names.extend(next_new_names)
                            i += 1
                    mods_condensed.append({
                        'type': tp,
                        'names': new_names,
                        'anchors': anchors,
                        'resindex': relixs
                    })
                    i += 1
                
                hdb_map[resname] = mods_condensed
            
            if s.tokens is None:
                break

