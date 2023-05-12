import os
import csv
import gzip
import argparse
import warnings
import pickle
from pathlib import Path

from poke import Poke, Pocket
from poke.writer import PocketWriter
from MDAnalysis import Universe, AtomGroup, Writer

structure_suffixes = ('.pdb', '.gro', '.pdbqt')
trajectory_suffixes = ('.xtc', '.trr')
list_suffixes = ('.txt',)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('structures', nargs='+', type=Path, help=
                        'input strctures (.pdb, .gro, .pdbqt, .xtc, ...).\n'
                        'List of proteins to search for pockets.\n'
                        'When coordinate-only file is passed the previous stucture '
                        ' is used to open it. '
                        'A txt file can be passed containing additional structures one per line.')     
    parser.add_argument('-o', '--output', type=Path, metavar='', default=Path(''), help=
                        'output path where to create the folders with the results.')
    parser.add_argument('-d', '--descriptors', action='store_true', default=False, help=
                        'compute the descriptors and save them in a csv file.')
    parser.add_argument('-v', '--volumes', action='store_true', default=False, help=
                        'save the volumes and save them in a multi-model PDB file.')
    parser.add_argument('-e', '--extract', action='store_true', default=False, help=
                        'save the results relative to a new trajectory containing only frames with at least one pocket.')
    parser.add_argument('-p', '--pymol', action='store_true', default=False, help=
                        'save a pymol script usable to open the pockets found.')
    parser.add_argument('-k', '--pickle', action='store_true', default=False, help=
                        'save a pickle with a list of pockets for each frame.')
    parser.add_argument('-q', '--quiet', action='store_true', default=False, help=
                        'avoid to writing the progess for long trajectories.')
    parser.add_argument('-f', '--filter', type=str, metavar='', default=None, help=
                        'filter the pockets with a condition.')
    parser.add_argument('-z', '--compress', action='store_true', default=False, help=
                        'compress the results using gzip, most of the programs can read them without problem.')
    

    args = parser.parse_args()

    try:
        for structures in split_structures(args.structures):
            out_dir = args.output / structures[0].stem.split('.', maxsplit=1)[0]
            out_dir.mkdir(parents=True, exist_ok=True)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                u = Universe(*structures)
            atoms = u.atoms

            scouter = Poke().load_atoms_properties_from_mdanalysis(atoms)
            frames = []
            pockets = []
            
            for i, ts in enumerate(u.trajectory):
                frame_pockets = scouter.get_pockets(ts.positions)

                if i >= 10 and i % 10 == 0 and not args.quiet:
                    print(f'{" ".join(map(str, structures))} frame {i}', end='\r')
                
                if args.filter is not None:
                    frame_pockets = [p for p in frame_pockets if eval(args.filter, p.__dict__, {"volume": p.volume, "surface": p.surface})]

                if not args.extract or len(frame_pockets) > 0: 
                    frames.append(i)
                    pockets.append(frame_pockets)
            
            if i >= 10 and not args.quiet:
                print(f'{" ".join(map(str, structures))} {i} frames', end='\r')

            if args.extract:
                write_coordinates(out_dir, u, frames)

            if args.descriptors:
                write_desciptors(out_dir, atoms, pockets, frames)

            if args.volumes:
                if args.extract:
                    write_points(args, out_dir, [structures[0], out_dir / 'frames.xtc'], pockets)
                else:
                    write_points(args, out_dir, structures, pockets)

            if args.pickle:
                write_pickle(args, out_dir, pockets)
        
    except UserError as error:
        print(error)


def split_structures(structures: list[Path]):
    output_list = None

    for i, s in enumerate(structures, 1):
        if s.suffix in structure_suffixes:
            if output_list:
                yield output_list            
            output_list = [s]
        
        elif s.suffix in trajectory_suffixes:
            if output_list:
                output_list.append(s)
            else:
                raise UserError('Before a trajectory file, a structure file must be specified.')

        elif s.suffix in list_suffixes:
            structures[i:i] = [Path(l.strip()) for l in s.read_text().splitlines() if l.strip()] 
        
        else:
            raise UserError(f'{s} has an unknown extension.')

    if output_list:
        yield output_list

def write_coordinates(out_dir, u: Universe, frames: list[int]):
    with Writer(str(out_dir / 'frames.xtc'), multiframe=True) as stream:
        for frame in frames:
            u.trajectory[frame]
            stream.write(u.atoms)

def write_desciptors(out_dir, atoms: AtomGroup, pockets: list[Pocket], frames: list[int]):
    with open(out_dir / 'descriptors.csv', 'wt') as stream:
        writer = csv.writer(stream)

        writer.writerow(('OutputFrame', 'OriginalFrame', 'Volume', 'Surface', 'Lenght1', 'Length2', 'Length3', 'Quadrupole1', 'Quadrupole2', 'Quadrupole3',
                        'Hydrophobicity', 'Enclousure', '#Regions', '#PolarAtoms', '#ChargedAtoms', 
                        'Residues', 'Atoms'))

        for frame, pks in enumerate(pockets):
            writer.writerows((frame, frames[frame], f'{p.volume:.2f}', f'{p.surface:.2f}', 
                              *(f'{l:.2f}' for l in p.lengths), *(f'{l:.2f}' for l in p.quadrupoles), 
                              f'{p.hydrophobicity:.2f}', f'{p.enclousure:.2f}', p.n_regions, 
                              p.n_polar_atoms, p.n_charges_atoms, 
                              ' '.join(map(str, atoms[p.atoms].residues.resids)),
                              ' '.join(map(str, p.atoms))) for p in pks)


def write_points(args, out_dir, structures: list[Path],  pockets: list[list[Pocket]]):
    with PocketWriter(open_stream(out_dir / 'volumes.pdb', 'wt', args.compress)) as writer:
        for pks in pockets:
            writer.write_frame(pks)

        if args.pymol:
            writer.write_pymol_script(out_dir / 'pymol.pml', structures)#, {'ligand': Path('1gkc/ligand.pdb')})

def write_pickle(args, out_dir, pockets: list[list[Pocket]]):
    with open_stream(out_dir / 'pockets.pickle', 'wb', args.compress) as stream:
        pickle.dump(pockets, stream)


def open_stream(path: os.PathLike,  mode: str, compressed: bool):
    if compressed:
        return gzip.open(str(path) + '.gz', mode)
    else:
        return open(path, mode)


class UserError(RuntimeError):
    pass