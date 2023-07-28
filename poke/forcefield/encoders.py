

def encode_residues(residues: dict):
    return '{\n' + ',\n'.join(encode_restype(restype, residues[restype]) for restype in sorted(residues)) + '\n}'

def encode_restype(resname: str, resdata: set):
    return '\t"' + resname + '": {\n\t\t"params": ' + str(list(resdata['params']))  + ',\n\t\t"names": ' + encode_names_list(resdata['names']) + '\n\t}'

def encode_names_list(names_list: list[list[str]]):
    return '[\n' + ',\n'.join(f'\t\t\t{encode_names(names)}' for names in names_list) + '\n\t\t]'

def encode_names(names: list[str]):
    return "[" + ', '.join(f'"{n}"' for n in names) + "]"

def encode_bonds(bonds: dict):
    return '{\n' + ',\n'.join(encode_resname_bonds(restype, bonds[restype]) for restype in sorted(bonds)) + '\n}'

def encode_resname_bonds(resname: str, resdata: set):
    return '\t"' + resname + '": ' + str(resdata)