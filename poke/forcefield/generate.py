import sys
import argparse
import numpy as np
from pathlib import Path
from .gmxtop import parse_ffnonbonded, parse_rtp, parse_tdb
from .process import define_atomtypes, add_alternative_names, amber_rules, desmond_rules
from .encoders import encode_residues, encode_bonds

parser = argparse.ArgumentParser()
parser.add_argument('gmx_forcefield', type=Path, help="Path of the gromacs forcefield")
args = parser.parse_args()

ffpath = Path(args.gmx_forcefield)
out_dir = Path(__file__).parent.parent / 'charmm36m'

resnames = {}
atom_types = parse_ffnonbonded(ffpath)
aminoacids_residue_types = parse_rtp(ffpath / 'aminoacids.rtp', atoms_types=atom_types, resnames=resnames)
nter_aminoacids_residue_types = parse_tdb(ffpath / 'aminoacids.n.tdb', residue_types=aminoacids_residue_types, atoms_types=atom_types, resnames=resnames)
cter_aminoacids_residue_types = parse_tdb(ffpath / 'aminoacids.c.tdb', residue_types=aminoacids_residue_types, atoms_types=atom_types, resnames=resnames)

na_residue_types = parse_rtp(ffpath / 'na.rtp', atoms_types=atom_types, resnames=resnames)
nter_na_residue_types = parse_tdb(ffpath / 'na.n.tdb', residue_types=na_residue_types, atoms_types=atom_types, resnames=resnames)
cter_na_residue_types = parse_tdb(ffpath / 'na.c.tdb', residue_types=na_residue_types, atoms_types=atom_types, resnames=resnames)

residue_types = {
    **aminoacids_residue_types, 
    **nter_aminoacids_residue_types, 
    **cter_aminoacids_residue_types,
    **nter_na_residue_types,
    **cter_na_residue_types
}

resnames['ARG'].extend(resnames['ARGN'])
resnames['ASP'].extend(resnames['ASPP'])
resnames['CYS'].extend(resnames['CYS2'])
resnames['GLU'].extend(resnames['GLUP'])
resnames['HIS'] = [*resnames['HSD'], *resnames['HSE'], *resnames['HSP']]
resnames['LYS'].extend(resnames['LSN'])
resnames['HEM'] = resnames['HEME']

atoms, residues, bonds = define_atomtypes(residue_types)


add_alternative_names(residues, amber_rules)
add_alternative_names(residues, desmond_rules)

np.savetxt(out_dir / 'atoms.csv', atoms, fmt='%.3f', delimiter=',')

with open(out_dir / 'residues.json', 'w') as file:
    file.write(encode_residues(residues))

with open(out_dir / 'bonds.json', 'w') as file:
    file.write(encode_bonds(bonds))