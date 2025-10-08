amber_rules = {
    "aminoacids": {
        "default": {
            "HN": "H",
        },
        ("SER", "CYS"): {
            "HN": "H",
            "HG1": "HG"
        },
    },

    "others": {
        "NME" : {
            "HN": "H",
        }
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
            "OT1": "O",
            "OT2": "OXT",
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
            "OT1": "O",
            "OT2": "OXT",
        },
        "ASPP": {
            "HN": "H",
            "HB1": "HB2",
            "HB2": "HB3",
            "HB3": "HB1",
            "OT1": "O",
            "OT2": "OXT",
        },
        ("TRP", "HSE", "HSE", "HSD", "HSP", "TYR", "PHE"): {
            "HN": "H",
            "HB1": "HB2",
            "HB2": "HB3",
            "HB3": "HB1",
            "OT1": "O",
            "OT2": "OXT",
        },
        "SER": {
            "HN": "H",
            "HG1": "HG",
            "HB1": "HB2",
            "HB2": "HB3",
            "HB3": "HB1",
            "OT1": "O",
            "OT2": "OXT",
        },
        "ILE": [{
            "HN": "H",
            "HG11": "HG12",
            "HG12": "HG13",
            "HG13": "HG11",
            "OT1": "O",
            "OT2": "OXT",
        },{
            "HN": "H",
            "CD": "CD1",
            "HD1": "HD11",
            "HD2": "HD12",
            "HD3": "HD13",
            "HG11": "HG12",
            "HG12": "HG13",
            "HG13": "HG11",
            "OT1": "O",
            "OT2": "OXT",
        }],
        "GLUP": {
            "HN": "H",
            "HB1": "HB2",
            "HB2": "HB3",
            "HB3": "HB1",
            "HG1": "HG2",
            "HG2": "HG3",
            "HG3": "HG1",
            "OT1": "O",
            "OT2": "OXT",
        },
        "CYS": {
            "HN": "H",
            "HB1": "HB2",
            "HB2": "HB3",
            "HB3": "HB1",
            "HG1": "HG",
            "OT1": "O",
            "OT2": "OXT",
        },
    },

    "others": {
        "ACE": {
            "HH31": "H1",
            "HH32": "H2",
            "HH33": "H3",
        },
        "NME": {
            "CH3": "CA", 
            "HH31": "HA1", 
            "HH32": "HA2", 
            "HH33": "HA3", 
            "HN": "H", 
        },
    }
}

charmm_legacy_rules = {
    "others": {
        "NME": {
            "CH3": "CAT",
            "HH31": "HT1",
            "HH32": "HT2",
            "HH33": "HT3",
            "HN": "HNT",
            "N": "NT"
        }
    }
}


charmm_metal_ions_rules = {
    "others": {
        "AG1P": {"AG1P": "AG"},
        "AG2P": {"AG2P": "AG"},
        "AL3P": {"AL3P": "AL"},
        "AU1P": {"AU1P": "AU"},
        "AU3P": {"AU3P": "AU"},
        "BE2P": {"BE2P": "BE"},
        "BI3P": {"BI3P": "BI"},
        "CE3P": {"CE3P": "CE"},
        "CO2P": {"CO2P": "CO"},
        "CO3P": {"CO3P": "CO"},
        "CR2P": {"CR2P": "CR"},
        "CR3P": {"CR3P": "CR"},
        "CU1P": {"CU1P": "CU"},
        "CU2P": {"CU2P": "CU"},
        "DY3P": {"DY3P": "DY"},
        "ER3P": {"ER3P": "ER"},
        "EU2P": {"EU2P": "EU"},
        "EU3P": {"EU3P": "EU"},
        "FE2P": {"FE2P": "FE"},
        "FE3P": {"FE3P": "FE"},
        "GA3P": {"GA3P": "GA"},
        "GD3P": {"GD3P": "GD"},
        "HG2P": {"HG2P": "HG"},
        "HO3P": {"HO3P": "HO"},
        "IN3P": {"IN3P": "IN"},
        "LA3P": {"LA3P": "LA"},
        "LU3P": {"LU3P": "LU"},
        "MN2P": {"MN2P": "MN"},
        "MN3P": {"MN3P": "MN"},
        "ND3P": {"ND3P": "ND"},
        "NI2P": {"NI2P": "NI"},
        "PB2P": {"PB2P": "PB"},
        "PD2P": {"PD2P": "PD"},
        "PM3P": {"PM3P": "PM"},
        "PR3P": {"PR3P": "PR"},
        "PT2P": {"PT2P": "PT"},
        "PU3P": {"PU3P": "PU"},
        "RA2P": {"RA2P": "RA"},
        "RH3P": {"RH3P": "RH"},
        "SC3P": {"SC3P": "SC"},
        "SM2P": {"SM2P": "SM"},
        "SM3P": {"SM3P": "SM"},
        "SN2P": {"SN2P": "SN"},
        "SR2P": {"SR2P": "SR"},
        "TB3P": {"TB3P": "TB"},
        "TI1P": {"TI1P": "TI"},
        "TI3P": {"TI3P": "TI"},
        "TL3P": {"TL3P": "TL"},
        "U3P": {"U3P": "U"},
        "V2P": {"V2P": "V"},
        "V3P": {"V3P": "V"},
        "Y3P": {"Y3P": "Y"},
        "YB2P": {"YB2P": "YB"},
        "YB3P": {"YB3P": "YB"},
    }
}


alternative_residue_names = {
    "HEM": ("HEME",),
    "NMA": ("NME",),

    # Amber
    "ASH": ("ASPP",),
    "GLH": ("GLUP",),
    "CYX": ("CYS2",),
    "HID": ("HSD",),
    "HIE": ("HSE",),
    "HIP": ("HSP",),
    "LYN": ("LSN",),

    # Metal Ions
    "AG1": ("AG1P",),
    "AG2": ("AG2P",),
    "AL3": ("AL3P",),
    "AU1": ("AU1P",),
    "AU3": ("AU3P",),
    "BE2": ("BE2P",),
    "BI3": ("BI3P",),
    "CE3": ("CE3P",),
    "CO2": ("CO2P",),
    "CO3": ("CO3P",),
    "CR2": ("CR2P",),
    "CR3": ("CR3P",),
    "CU1": ("CU1P",),
    "CU2": ("CU2P",),
    "DY3": ("DY3P",),
    "ER3": ("ER3P",),
    "EU2": ("EU2P",),
    "EU3": ("EU3P",),
    "FE2": ("FE2P",),
    "FE3": ("FE3P",),
    "GA3": ("GA3P",),
    "GD3": ("GD3P",),
    "HG2": ("HG2P",),
    "HO3": ("HO3P",),
    "IN3": ("IN3P",),
    "LA3": ("LA3P",),
    "LU3": ("LU3P",),
    "MN2": ("MN2P",),
    "MN": ("MN2P",),
    "MN3": ("MN3P",),
    "ND3": ("ND3P",),
    "NI2": ("NI2P",),
    "PB2": ("PB2P",),
    "PD2": ("PD2P",),
    "PM3": ("PM3P",),
    "PR3": ("PR3P",),
    "PT2": ("PT2P",),
    "PU3": ("PU3P",),
    "RA2": ("RA2P",),
    "RH3": ("RH3P",),
    "SC3": ("SC3P",),
    "SM2": ("SM2P",),
    "SM3": ("SM3P",),
    "SN2": ("SN2P",),
    "SR2": ("SR2P",),
    "TB3": ("TB3P",),
    "TI1": ("TI1P",),
    "TI3": ("TI3P",),
    "TL3": ("TL3P",),
    "U3": ("U3P",),
    "V2": ("V2P",),
    "V3": ("V3P",),
    "Y3": ("Y3P",),
    "YB2": ("YB2P",),
    "YB3": ("YB3P",),
}


protonation_states = {
    "ARG": ("ARG", "ARGN"),
    "ASP": ("ASP", "ASPP"),
    "GLU": ("GLU", "GLUP"),
    "CYS": ("CYS", "CYS2"),
    "HIS": ("HSD", "HSE", "HSP"),
    "LYS": ("LYS", "LSN"),
}


def add_additional_mapping(residues_map: dict, rules: dict[str, list[str]], skip_existing: bool = False):
  
    for resname, additional in rules.items():
        if resname in residues_map and skip_existing:
            print(f"Warning: {resname} already exists in residues_map, and skip_existing is True. Skipping update!")
            continue
        
        mapping = residues_map.setdefault(resname, [])
        for add_name in additional:  
            if add_name in residues_map:
                adds = residues_map[add_name]
                for add in adds:
                    if add not in mapping:
                        mapping.append(add)
            else:
                raise ValueError(f"Residue '{add_name}' not found in residues_map")

            

def add_alternative_atom_names(residues, rules: dict[str, dict[str,str]]):
    
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
            specific_rules = rules.get('aminoacids')

        else:
            specific_rules = rules.get('others')

        if specific_rules:
            replacers = specific_rules.get(restype) or \
                        specific_rules.get(restype.split('-')[0]) or \
                        specific_rules.get('default')
            
            if replacers:                
                if not isinstance(replacers, list):
                    replacers = [replacers]

                for replacer in replacers:
                    replaced_names = [replacer.get(name,name) for name in charmm_names]

                    if len(set(replaced_names)) < len(replaced_names):
                        print(f'Warning: Duplicate names found when applying {restype} {charmm_names} -> {replaced_names}')
                        continue
                        #raise RuntimeError(f'{restype} {charmm_names} -> {replaced_names}')

                    if replaced_names not in residue['names']:
                        residue['names'].append(replaced_names)


def extract_data(residues_map: dict):
    atoms = {}

    for data_list in residues_map.values():
        for data in data_list:
            atoms_data = data['atoms']

            for key in atoms_data.values():
                atoms[key] = atoms.get(key, 0) + 1

    atoms = sorted(atoms, key=lambda k: atoms[k], reverse=True)
    resnames = {}
    residues = {}
    bonds = {}

    for data_list in residues_map.values():
        
        repr_names = []
        for data in data_list:
            repr_name = min(data['resnames'], key=lambda x: len(x))
            atoms_data = data['atoms']
            bonds_data = data['bonds']

            names = [sorted(atoms_data.keys())]
            params = [atoms.index(atoms_data[n]) for n in names[0]]

            residues[repr_name] = {
                "params": params,
                "names": names
            }

            bonds[repr_name] = [[names[0].index(n1), names[0].index(n2)] for n1, n2 in bonds_data]
            repr_names.append(repr_name)
        
        # First residue data has the names to map
        for resname in data_list[0]['resnames']:
            resnames[resname] = repr_names
    
    
    return atoms, residues, bonds, resnames