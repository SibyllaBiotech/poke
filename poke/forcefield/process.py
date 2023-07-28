
amber_rules = {
    "aminoacids": {
        "default": {
            "HN": "H",
        },
        ("SER", "CYS"): {
            "HN": "H",
            "HG1": "HG"
        },
    }
}

desmond_rules = {
    "aminoacids": {
        "default": {
            "HN": "H",
            "HA1": "HA2",
            "HA2": "HA3",
            "HA3": "HA1",
            "HB1": "HB2",
            "HB2": "HB3",
            "HB3": "HB1",
            "HD1": "HD2",
            "HD2": "HD3",
            "HD3": "HD1",
            "HE1": "HE2",
            "HE2": "HE3",
            "HE3": "HE1",
            "HG1": "HG2",
            "HG2": "HG3",
            "HG3": "HG1",
        },
        "ARGN": {
            "HN": "H",
            "HB1": "HB2",
            "HB2": "HB3",
            "HB3": "HB1",
            "HD1": "HD2",
            "HD2": "HD3",
            "HD3": "HD1",
            "HG1": "HG2",
            "HG2": "HG3",
            "HG3": "HG1",
            "NH1": "NH2",
            "NH2": "NH1",
            "HH11": "HH2",
            "HH21": "HH11",
            "HH22": "HH12",
        },
        "ASPP": {
            "HN": "H",
            "HB1": "HB2",
            "HB2": "HB3",
            "HB3": "HB1"
        },
        ("TRP", "HSE", "HSE", "HSD", "HSP", "TYR", "PHE"): {
            "HN": "H",
            "HB1": "HB2",
            "HB2": "HB3",
            "HB3": "HB1"
        },
        "SER": {
            "HN": "H",
            "HG1": "HG",
            "HB1": "HB2",
            "HB2": "HB3",
            "HB3": "HB1"
        },
        "ILE": [{
            "HN": "H",
            "HG11": "HG12",
            "HG12": "HG13",
            "HG13": "HG11",
        },{
            "HN": "H",
            "CD": "CD1",
            "HD1": "HD11",
            "HD2": "HD12",
            "HD3": "HD13",
            "HG11": "HG12",
            "HG12": "HG13",
            "HG13": "HG11",
        }],
        "GLUP": {
            "HN": "H",
            "HB1": "HB2",
            "HB2": "HB3",
            "HB3": "HB1",
            "HG1": "HG2",
            "HG2": "HG3",
            "HG3": "HG1",
        },
        "CYS": {
            "HN": "H",
            "HB1": "HB2",
            "HB2": "HB3",
            "HB3": "HB1",
            "HG1": "HG",
        },
    }
}

def add_alternative_names(residues, rules: dict[str, dict[str,str]]):
    
    for key1, value1 in rules.items():
        new_value1 = {}
        for key, value in value1.items():
            if isinstance(key, tuple):
                for k in key:
                    new_value1[k] = value
            else:
                new_value1[key] = value
        rules[key1] = new_value1
 
    for restype, residue in residues.items():
        charmm_names = residue['names'][0]

        if set(['C', 'CA', 'N']) < set(charmm_names):
            amino_rules = rules['aminoacids']
            replacers = amino_rules.get(restype) or \
                        amino_rules.get(restype.split('-')[0]) or \
                        amino_rules.get('default')
            
            if not isinstance(replacers, list):
                replacers = [replacers]

            for replacer in replacers:
                replaced_names = sorted([replacer.get(name,name) for name in charmm_names])

                if len(set(replaced_names)) < len(replaced_names):
                    raise RuntimeError(f'{restype} {charmm_names} -> {replaced_names}')

                if replaced_names not in residue['names']:
                    residue['names'].append(replaced_names)


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