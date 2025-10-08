import sys
import argparse
import numpy as np
from pathlib import Path
from .gmxtop import parse_ffnonbonded, parse_rtp, parse_tdb, parse_hdb
from .process import (extract_data, add_alternative_atom_names, add_additional_mapping, 
                      amber_rules, desmond_rules, charmm_legacy_rules, charmm_metal_ions_rules, alternative_residue_names, protonation_states)
from .encoders import json_dumps

parser = argparse.ArgumentParser()
parser.add_argument('gmx_forcefield', type=Path, help="Path of the gromacs forcefield")
args = parser.parse_args()

ffpath = Path(args.gmx_forcefield)
out_dir = Path(__file__).parent.parent / 'charmm36m'
out_dir.mkdir(parents=True, exist_ok=True)

termini = {}
hydrogens = {}

print(f'Processing forcefield from {ffpath} to {out_dir}')
print('Reading atom types...')
atom_types = parse_ffnonbonded(ffpath)

print('Reading residues...')
residues_map = {}
parse_rtp(ffpath / 'aminoacids.rtp',   atoms_types=atom_types, residues_map=residues_map)
parse_tdb(ffpath / 'aminoacids.n.tdb', atoms_types=atom_types, residues_map=residues_map, mods_spec=termini)
parse_tdb(ffpath / 'aminoacids.c.tdb', atoms_types=atom_types, residues_map=residues_map, mods_spec=termini)
parse_hdb(ffpath / 'aminoacids.hdb', hdb_map=hydrogens)

print('Reading nucleic acids...')
na_map = {}
parse_rtp(ffpath / 'na.rtp',   atoms_types=atom_types, residues_map=na_map)
parse_tdb(ffpath / 'na.n.tdb', atoms_types=atom_types, residues_map=na_map, mods_spec=termini)
parse_tdb(ffpath / 'na.c.tdb', atoms_types=atom_types, residues_map=na_map, mods_spec=termini)
parse_hdb(ffpath / 'na.hdb', hdb_map=hydrogens)

print('Reading carbohydrates...')
carb_map = {}
parse_rtp(ffpath / 'carb.rtp', atoms_types=atom_types, residues_map=carb_map)
parse_hdb(ffpath / 'carb.hdb', hdb_map=hydrogens)

print('Reading metals...')
metals_map = {}
parse_rtp(ffpath / 'metals.rtp', atoms_types=atom_types, residues_map=metals_map)

print('Reading solvent...')
solvent_map = {}
parse_rtp(ffpath / 'solvent.rtp', atoms_types=atom_types, residues_map=solvent_map)
parse_hdb(ffpath / 'solvent.hdb', hdb_map=hydrogens)


all_map = {
    **residues_map,
    **na_map,
    **carb_map,
    **metals_map,
    **solvent_map
}


# Also insert PRO terminus in NH3+, to have a single name to refer to charged termini
termini['NH3+']['PRO'] = termini['NH2+']['PRO']
del termini['NH3+']['GLY'] # This is duplicated for our purposes

atoms, residues, bonds, resnames = extract_data(all_map)

# Add additional residue mappings
add_additional_mapping(resnames, alternative_residue_names, skip_existing=True)
add_additional_mapping(resnames, protonation_states)

# Add alternative atom naming 
add_alternative_atom_names(residues, amber_rules)
add_alternative_atom_names(residues, desmond_rules)
add_alternative_atom_names(residues, charmm_legacy_rules)
add_alternative_atom_names(residues, charmm_metal_ions_rules)

# Save all forcefield data
np.savetxt(out_dir / 'atoms.csv', atoms, fmt='%.3f', delimiter=',')

Path(out_dir / 'residues.json').write_text(
    json_dumps(residues, sort=True, inline_limit=0)
)

Path(out_dir / 'bonds.json').write_text(
    json_dumps(bonds, sort=True)
)

Path(out_dir / 'resnames.json').write_text(
    json_dumps(resnames, sort=True)
)

Path(out_dir / 'termini.json').write_text(
    json_dumps(termini, sort=True)
)

Path(out_dir / 'hydrogens.json').write_text(
    json_dumps(hydrogens, sort=True)
)

print('Forcefield data generated successfully!')