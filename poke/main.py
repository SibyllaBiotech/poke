import os
import csv
import gzip
import argparse
import warnings
import pickle
from pathlib import Path

from poke import Poke, Pocket
from poke.filter import PocketEvaluator
from poke.writer import PocketWriter
from MDAnalysis import Universe, AtomGroup, Writer

warnings.filterwarnings('ignore', category=UserWarning, module='MDAnalysis')

structure_suffixes = ('.pdb', '.gro', '.pdbqt')
trajectory_suffixes = ('.xtc', '.trr')
list_suffixes = ('.txt',)
resid_start = 800

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
    parser.add_argument('-a', '--atoms', type=str, metavar='', default='all', help=
                        'selection of atoms used by the scouter.')
    parser.add_argument('-e', '--extract', action='store_true', default=False, help=
                        'save the results relative to a new trajectory containing only frames with at least one pocket.')
    parser.add_argument('-s', '--sort', type=str, metavar='', default='all', help=
                        'selection of atoms used by the scouter.')
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
        # args.structures contains all the inputs in sequence
        # split_structures produce a iterator of chunks to be passed to the Universe constructor
        for structures in split_structures(args.structures):

            # Suppress the warning during the opening of the structures
            u = Universe(*structures)
            atoms = u.select_atoms(args.atoms)

            # Initialize the scouter with the atoms 
            scouter = Poke().load_atoms_properties_from_mdanalysis(atoms)
            
            frames = []
            pockets = []
            
            # Run the scouter for each frame
            for i, ts in enumerate(u.trajectory):
                frame_pockets = scouter.get_pockets(atoms.positions)

                # When there are many frame write a progress if not quiet
                if i >= 10 and i % 10 == 0 and not args.quiet:
                    print(f'{" ".join(map(str, structures))} frame {i}', end='\r')
                
                # Apply a filter to the pockets evaluating the args.filter on each pocket
                if args.filter is not None:
                    frame_pockets = [p for p in frame_pockets if PocketEvaluator.evaluate(args.filter, p, atoms)]

                # If sorting is active only the trajectory will be written with only one pocket per frame
                if args.sort:
                    for p in frame_pockets:
                        frames.append(i)
                        pockets.append([p])

                # If the args.extract is active, the element is added only if at least one pocket was found
                elif not args.extract or len(frame_pockets) > 0:
                    frames.append(i)
                    pockets.append(frame_pockets)
        
            # Complete the progress writing the total number of frames analyzed
            if i >= 10 and not args.quiet:
                print(f'{" ".join(map(str, structures))} {i} frames', end='\n')

            # Sort the pockets
            if args.sort:
                frames, pockets = zip(*sorted(zip(frames, pockets), key=lambda fp: PocketEvaluator.evaluate(args.sort, fp[1][0], atoms), reverse=True))

            # Define the output folder for this chunk using the stem name of the 
            # first structure in the chunk that usually is the reference (.pdb / .gro)
            out_dir = args.output / structures[0].stem.split('.', maxsplit=1)[0]
            out_dir.mkdir(parents=True, exist_ok=True)

            # If args.extract is set, write a new trajectory with all
            # the frames collected in the list (with at least a pocket)
            if args.extract or args.sort:
                write_coordinates(out_dir, atoms, frames)

            if args.descriptors:
                write_desciptors(out_dir, atoms, pockets, frames)

            if args.volumes:                
                if args.extract:
                    # If the args.extract flag is active the pymol script needs a 
                    # trajectory with only the frames written by write_coordinates
                    write_points(args, out_dir, [out_dir / 'structure.pdb', out_dir / 'frames.xtc'], pockets)
                else:
                    write_points(args, out_dir, structures, pockets)

            if args.pickle:
                write_pickle(args, out_dir, pockets)
        
    except UserError as error:
        print(error)    

def split_structures(structures: list[Path]):
    '''
        Split the list of inputs in chunks.
        Every chunk start with a file with structure_suffixes.
        and continues until the next structure_suffixes.
        If it finds a .txt file, add its lines to the input list.
    '''
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

def write_coordinates(out_dir, atoms: AtomGroup, frames: list[int]):
    '''
        Write an xtc trajectory with a subset of 
        the frames specified in the list frames.
    '''
    atoms.universe.trajectory[0]
    atoms.write(str(out_dir / 'structure.pdb'))

    with Writer(str(out_dir / 'frames.xtc'), multiframe=True) as stream:
        for frame in frames:
            atoms.universe.trajectory[frame]
            stream.write(atoms)

def write_desciptors(out_dir, atoms: AtomGroup, pockets: list[Pocket], frames: list[int]):
    '''
        Write a csv with the descriptors for each pocket in the list.
    '''
    with open(out_dir / 'descriptors.csv', 'wt') as stream:
        writer = csv.writer(stream)

        header_written = False
        for frame, pks in enumerate(pockets):
            for id, p in enumerate(pks, resid_start):
                props = p.props()
                
                if not header_written:
                    writer.writerow((
                        'OutputFrame', 
                        'OriginalFrame', 
                        'PocketResid', 
                        *props.keys() ,
                        'Residues', 
                        'Atoms'
                    ))
                    
                    header_written = True

                writer.writerow((
                    frame, 
                    frames[frame], 
                    id,  
                    *(f'{v:.2g}' for v in props.values()),
                    ' '.join(map(str, atoms[p.atoms].residues.resids)),
                    ' '.join(map(str, p.atoms))
                ))


def write_points(args, out_dir, structures: list[Path],  pockets: list[list[Pocket]]):
    with PocketWriter(open_stream(out_dir / 'volumes.pdb', 'wt', args.compress), resid_start=resid_start) as writer:
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