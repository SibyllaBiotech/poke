
amber_rules = {
    "aminoacids": {
        "default": {
            "HN": "H",
            "OT1": "OC1",
            "OT2": "OC2"
        },
        "SER": {
            "HN": "H",
            "OT1": "OC1",
            "OT2": "OC2",
            "HG1": "HG"
        },
        "CYS": {
            "HN": "H",
            "OT1": "OC1",
            "OT2": "OC2",
            "HG1": "HG"
        }
    }
}

def add_alternative_names(residues, rules: dict[str, dict[str,str]]):

    for restype, residue in residues.items():
        charmm_names = residue['names'][0]

        if set(['C', 'CA', 'N']) < set(charmm_names):
            amino_rules = rules['aminoacids']
            replacer = amino_rules.get(restype) or \
                       amino_rules.get(restype.split('-')[0]) or \
                       amino_rules.get('default')

            replaced_names = [replacer.get(name,name) for name in charmm_names]

            if replaced_names not in residue['names']:
                residue['names'].append(replaced_names)



def add_chimera_desmond(resname, names):
    amino_map = {
        'HN': 'H',
        'OT1': 'OC1',
        'OT2': 'OC2',
        "HA1": "HA2", "HA2": "HA3",
        "HB1": "HB2", "HB2": "HB3",
        "HE1": "HE2", "HE2": "HE3",
        "HG1": "HG2", "HG2": "HG3",
    }
    
    charmm_names = names[0]

    if 'N' in charmm_names and 'C' in charmm_names and 'CA' in charmm_names:
        new_names = [amino_map.get(n,n) for n in charmm_names]
     
        if 'SER' in resname or 'CYS' in resname:
            new_names = ['HG' if n == 'HG2' else n for n in new_names]

        elif 'ILE' in resname:
            ile_map = {
                'CD': 'CD1',
                'HD1': 'HD11',
                'HD2': 'HD12',
                'HD3': 'HD13',
                "HG11": "HG12", "HG12": "HG13",
            }
            new_names = [ile_map.get(n,n) for n in new_names]

        elif 'LYS' in resname:
            lys_map = {
                'HD1': 'HD2',
                'HD2': 'HD3',
            }
            new_names = [lys_map.get(n,n) for n in new_names]

        if new_names not in names:
            names.append(new_names)  


def define_atomtypes(residue_types):
    atoms = {}

    for resname, resdata in residue_types.items():
        atoms_data = resdata['atoms']

        for key in atoms_data.values():
            atoms[key] = atoms.get(key, 0) + 1

    atoms = sorted(atoms, key=lambda k: atoms[k], reverse=True)
    residues = {}
    bonds = {}

    for restype, resdata in residue_types.items():
        atoms_data = resdata['atoms']
        bonds_data = resdata['bonds']

        names = [sorted(atoms_data.keys())]
        params = [atoms.index(atoms_data[n]) for n in names[0]]

        residues[restype] = {
            "params": params,
            "names": names
        }

        bonds[restype] = [[names[0].index(n1), names[0].index(n2)] for n1, n2 in bonds_data]
    return atoms, residues, bonds