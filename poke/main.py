import os
import sys
import csv
import gzip
import argparse
import warnings
import pickle
import multiprocessing
import numpy as np
from pathlib import Path
import shlex

from poke import Poke, Pocket, DEFAULT_MIN_VOLUME, DEFAULT_ACCESSIBLE_CUTOFF
from poke.filter import PocketEvaluator
from poke.writer import PocketWriter

from .algorithms import add_pocket_to_density

from .preparator import Preparator, MDAUniverse
from .parametrizer import charmm

from MDAnalysis import Universe, AtomGroup, Writer, Merge
from .utils import Aligner
from rich import progress
from rich.text import Text

warnings.filterwarnings('ignore', category=UserWarning, module='MDAnalysis')

structure_suffixes = ('.pdb', '.gro', '.pdbqt')
trajectory_suffixes = ('.xtc', '.trr', '.dcd')
list_suffixes = ('.txt',)
resid_start = 800

def build_parser() -> argparse.ArgumentParser:
    '''Build the command-line parser. Used by main() and by the documentation.'''
    parser = argparse.ArgumentParser(prog='poke', add_help=False)
    parser.add_argument('structures', nargs='+', type=Path, help=
                        'input structures (.pdb, .gro, .pdbqt, .xtc, ...). '
                        'List of proteins to search for pockets. '
                        'When a coordinate-only file is passed, the previous structure '
                        'is used to open it. '
                        'A .txt file can be passed containing additional structures, one per line.')

    # The help option is registered first so that it leads the usage line, as users expect;
    # the groups are created in the order their sections appear in the help text.
    output = parser.add_argument_group('output options')
    general = parser.add_argument_group('general options')
    general.add_argument('-h', '--help', action='help', help='show this help message and exit')
    output.add_argument('-d', '--descriptors', action='store_true', default=False, help=
                        'compute the descriptors and save them in a CSV file.')
    output.add_argument('--write-atoms', action='store_true', default=False, help=
                        'add the pocket atoms to the descriptors CSV file (requires -d).')
    output.add_argument('-v', '--volumes', action='store_true', default=False, help=
                        'save the volumes in a multi-model PDB file.')
    output.add_argument('-p', '--pymol', action='store_true', default=False, help=
                        'save a PyMOL script to open the pockets found.')
    output.add_argument('-dx', '--density', action='store_true', default=False, help=
                        'save the density map of the pockets in a DX file.')
    output.add_argument('-e', '--extract', action='store_true', default=False, help=
                        'save the results relative to a new trajectory containing only frames with at least one pocket.')
    output.add_argument('-k', '--pickle', action='store_true', default=False, help=
                        'save a pickle with a list of pockets for each frame.')
    output.add_argument('-g', '--graph', action='store_true', default=False, help=
                        'plot distributions of pocket descriptors.')
    output.add_argument('-ls', '--launch-script', action='store_true', default=False, help=
                        'save the launch command in a bash script to be able to rerun with the same parameters.')

    general.add_argument('-o', '--output', type=Path, metavar='PATH', default=Path(''), help=
                         'output path where to create the folders with the results.')
    general.add_argument('-q', '--quiet', action='store_true', default=False, help=
                         'avoid writing the progress for long trajectories.')
    general.add_argument('-nc','--ncores', type=int, metavar='N', default=None, help=
                         'number of cores used to parallelize the computation of trajectories or multi-model PDB files.')
    general.add_argument('-z', '--compress', action='store_true', default=False, help=
                         'compress the results using gzip; most programs can read them without problem.')

    selection = parser.add_argument_group('atom selection options')
    selection.add_argument('-ss', '--selection-scouter', type=str, metavar='SEL', default=None, help=
                           'selection of atoms used by the scouter.')
    selection.add_argument('-so', '--selection-output', type=str, metavar='SEL', default='all', help=
                           'selection of atoms written in the output.')

    alignment = parser.add_argument_group('alignment options')
    alignment.add_argument('-a', '--align', type=str, metavar='SEL', default=None, help=
                           'align the given atom selection to the reference.')
    alignment.add_argument('-r', '--reference', type=str, metavar='PATH', default=None, help=
                           'reference used for alignment. Default is the first structure.')

    detection = parser.add_argument_group('pocket detection options')
    detection.add_argument('-mv', '--min-volume', type=float, metavar='V', default=DEFAULT_MIN_VOLUME, help=
                           'minimum volume of the pocket core in cubic angstroms. '
                           'Pockets with a smaller volume are discarded. Default is %(default)s A^3.')
    detection.add_argument('-ac', '--accessible-cutoff', type=float, metavar='R', default=DEFAULT_ACCESSIBLE_CUTOFF, help=
                           'distance in angstroms around the pocket core within which the low-energy points '
                           'are counted as accessible volume. Default is %(default)s A.')

    scoring = parser.add_argument_group('scoring options')
    scoring.add_argument('-sf', '--scoring-function', type=str, metavar='NAME', default='model', help=
                         'scoring function used.')

    filtering = parser.add_argument_group('filtering and sorting options')
    filtering.add_argument('-f', '--filter', type=str, metavar='EXPR', default=None, help=
                           'filter the pockets with a condition.')
    filtering.add_argument('-s', '--sort', type=str, metavar='EXPR', default=None, help=
                           'sorting by the value of an expression.')
    filtering.add_argument('-n', '--max-output-frames', type=int, metavar='N', default=None, help=
                           'maximum number of pockets in the output. Only when --sort is used.')
    filtering.add_argument('-ks', '--keep-subpockets', action='store_true', default=False, help=
                           'keep all pockets even if they share regions, instead of only the highest scoring one.')

    preparation = parser.add_argument_group('structure preparation options')
    preparation.add_argument('--prepare', action='store_true', default=False, help=
                        'prepare the structures before scouting (select altLocs, add hydrogens and caps). '
                        'Default is to NOT prepare the structures.')
    preparation.add_argument('--keep-hydrogens', action='store_true', default=False, help=
                        'try to keep the existing hydrogens during preparation, adding only the missing ones.')
    preparation.add_argument('--selection-prepare', type=str, metavar='SEL', default=None, help=
                        'selection of atoms to prepare. Other atoms are kept as they are. '
                        'If not provided, only proteins, nucleic acids, '
                        'and common cofactors and structural ions are prepared.')

    capping = parser.add_argument_group('capping options')
    capping.add_argument('--caps-types', nargs='+', type=str, metavar='TYPE N_TERM C_TERM',
                        default='protein NH3+ COO- nucleic_acid 5TER 3TER'.split(), help=
                        'types of caps to add to the structures.')
    capping.add_argument('--caps-exclude', nargs='*', type=str, metavar='NAME', 
                        default='ACE NME'.split(), help=
                        'list of residue names to exclude from capping.')
    capping.add_argument('--ignore-segid', action='store_true', default=False, help=
                        'do not use the segid (chain) value to split the structure for capping.')
    capping.add_argument('--max-other-len', type=int, metavar='LEN', default=5, help=
                        'maximum number of residues of type "other" before a chain is considered interrupted.')
    capping.add_argument('--resid-threshold', type=int, metavar='LEN', default=50, help=
                        'jump in residue numbers to consider a chain interrupted (used only if segid info is absent).')
    capping.add_argument('--ter-indices', type=int, nargs='*', metavar='IX', default=[], help=
                        'explicit list of 0-based indices of residues to be considered as terminal residues ' 
                        'for splitting.')

    protonation = parser.add_argument_group('protonation state options')
    protonation.add_argument('--his', type=str, metavar='TYPE', default='auto', help=
                        'type of histidine to use. Options are "auto", "bond", "HSP", "HSE", "HSD". '
                        'Default is "auto" which uses the presence of HD1 and HE2 or hydrogen bonding detection.')
    protonation.add_argument('--lys', type=str, metavar='TYPE', default='auto', help=
                        'type of lysine to use. Options are "auto", "LYS", "LSN". Default is "auto".')
    protonation.add_argument('--arg', type=str, metavar='TYPE', default='auto', help=
                        'type of arginine to use. Options are "auto", "ARG", "ARGN". Default is "auto".')
    protonation.add_argument('--glu', type=str, metavar='TYPE', default='auto', help=
                        'type of glutamate to use. Options are "auto", "GLU", "GLUP". Default is "auto".')
    protonation.add_argument('--asp', type=str, metavar='TYPE', default='auto', help=
                        'type of aspartate to use. Options are "auto", "ASP", "ASPP". Default is "auto".')
    protonation.add_argument('--cys', type=str, metavar='TYPE', default='auto', help=
                        'type of cysteine to use. Options are "auto", "bond", "CYS", "CYS2". ' \
                        'Default is "auto" which uses the presence of HG1 or covalent bonding detection.')
    protonation.add_argument('--custom-restype', nargs='*', type=str, metavar='IX TYPE', help=
                        'custom residue types to use at specific 0-based residue indices.')
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    args.caps_types = {t: (cn, cc) for t, cn, cc in zip(args.caps_types[::3], args.caps_types[1::3], args.caps_types[2::3])}
    if args.custom_restype:
        args.custom_restype = {int(ix): t for ix, t in zip(args.custom_restype[::2], args.custom_restype[1::2])}
    
    # Precompile filter and sort expressions
    if args.filter:
        args.filter = compile(args.filter, __file__, 'eval')
    if args.sort:
        args.sort = compile(args.sort, __file__, 'eval')

    # Disable parallelization of numpy operations
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
    os.environ['NUMEXPR_NUM_THREADS'] = '1'

    try:
        # args.structures contains all the inputs in sequence
        # split_structures produce a iterator of chunks to be passed to the Universe constructor
        for structures in split_structures(args.structures):

            # Define the output folder for this chunk using the stem name of the 
            # first structure in the chunk that usually is the reference (.pdb / .gro)
            out_dir = args.output / '-'.join(fn.stem.split('.', maxsplit=1)[0] for fn in structures)

            # Run the main function in structures
            run(out_dir, structures, args)
        
    except UserError as error:
        print(error) 


def run(out_dir: Path, structures: list[Path], args: argparse.Namespace):
 
    # Print the current structures
    print(f'{" ".join(map(str, structures))} ▶ ./{out_dir}')

    if not args.prepare:
        u = Universe(*structures)
        if args.selection_scouter is None:
            args.selection_scouter = 'all'
            
    else:
        u = MDAUniverse(*structures,
                        selection=args.selection_prepare,
                        forcefield=charmm,
                        keep_existing=args.keep_hydrogens,
                        strictness=Preparator.StrictLevel.NONE,
                        caps_types=args.caps_types,
                        max_interruption_size=args.max_other_len,
                        on_segid=not args.ignore_segid,
                        resid_threshold=args.resid_threshold,
                        skip=args.caps_exclude,
                        ter= args.ter_indices,
                        his=args.his,
                        lys=args.lys,
                        arg=args.arg,
                        glu=args.glu,
                        asp=args.asp,
                        cys=args.cys,
                        specific=args.custom_restype)
        
        if args.selection_scouter is None:
            args.selection_scouter = u.selection
            print(f'\nStructure is being processed and no selection for scouting has been provided. '
                  f'Selection "{args.selection_scouter}" will be used for scouting.\n')

        # Preparation always implys extraction to write protonated structures
        args.extract = True 

    
    atoms_scouter = u.select_atoms(args.selection_scouter)
    atoms_output = u.select_atoms(args.selection_output)
    scouted_output_mask = np.isin(atoms_output.ix, atoms_scouter.ix)
    scouted_output_indexes = np.where(scouted_output_mask)[0]
    extra_output_indexes = np.where(~scouted_output_mask)[0]

    # Initialize copy of atom scouter and reference (needed for rmsd filtering)
    atoms_scouter_copy = atoms_scouter.universe.copy().atoms

    if args.reference is not None:
        reference = Universe(args.reference).atoms
    else:
        u.trajectory[0]
        reference = u.copy().atoms

    # Initialize the scouter with the atoms 
    scouter = Poke.from_mdanalysis(atoms_scouter, scoring_model=args.scoring_function, keep_subpockets=args.keep_subpockets,
                                   min_volume=args.min_volume, accessible_cutoff=args.accessible_cutoff)

    frames = []
    pockets = []

    with progress.Progress(
            progress.TextColumn('   '),
            progress.BarColumn(),
            progress.TimeElapsedColumn(),
            progress.TextColumn('eta'),
            progress.TimeRemainingColumn(),
            progress.MofNCompleteColumn(),
            FrameSpeedColumn(),
            PocketsFoundColumn(),
        ) as progress_bar:

        task = progress_bar.add_task("scouting", total=u.trajectory.n_frames, n_pockets=0, n_kept_pockets=0)
        n_kept_pockets = 0
        n_pockets = 0
    
        with Pool(scouter=scouter, n_frames=u.trajectory.n_frames, n_cores=args.ncores) as pool: 
            for i, frame_pockets in enumerate(pool.get_pockets(atoms_scouter.positions for _ in u.trajectory)):
                atoms_scouter_copy.universe.trajectory[i]
                n_pockets += len(frame_pockets)
                # Apply a filter to the pockets evaluating the args.filter on each pocket
                if args.filter is not None:
                    frame_pockets = [p for p in frame_pockets if PocketEvaluator.evaluate(args.filter, p, atoms_scouter_copy, reference)]
                n_kept_pockets += len(frame_pockets)

                # If sorting is active only the trajectory will be written with only one pocket per frame
                if args.sort:
                    for p in frame_pockets:
                        frames.append(i)
                        pockets.append([p])

                # If the args.extract is active, the element is added only if at least one pocket was found
                elif not args.extract or len(frame_pockets) > 0:
                    frames.append(i)
                    pockets.append(frame_pockets)
                
                if (i+1) % pool.chunksize == 0:
                    progress_bar.update(task, advance=pool.chunksize, n_kept_pockets=n_kept_pockets, n_pockets=n_pockets)
        
        progress_bar.update(task, completed=u.trajectory.n_frames, n_kept_pockets=n_kept_pockets, n_pockets=n_pockets)
    
    # print(f'Frames with pockets found: {len(frames)}')
    # print(f'Pockets found: {sum(len(pks) for pks in pockets)}')

    # Sort the pockets
    if args.sort and frames and pockets:
        frames, pockets = zip(*sorted(zip(frames, pockets), key=lambda fp: PocketEvaluator.evaluate(args.sort, fp[1][0], atoms_scouter, reference), reverse=True))

    if args.max_output_frames is not None:
        frames = frames[:args.max_output_frames]
        pockets = pockets[:args.max_output_frames]

    # Create the output directory
    out_dir.mkdir(parents=True, exist_ok=True)

    # Write atoms anv volumes
    if (args.extract or args.sort or args.volumes) and not args.density:
        write_atoms_and_pockets(args, structures, out_dir, atoms_output, frames, pockets, scouted_output_indexes, extra_output_indexes, reference)

    if (args.extract or args.sort or args.volumes) and args.density:
        write_atoms_and_pockets_and_density(args, structures, out_dir, atoms_output, frames, pockets, scouted_output_indexes, extra_output_indexes, reference)

    if args.descriptors:
        write_descriptors(out_dir, atoms_scouter, pockets, frames, write_atoms=args.write_atoms)

    if args.pickle:
        write_pickle(args, out_dir, pockets)

    if args.graph:
        plot_descriptors(out_dir, pockets)
    
    if args.launch_script:
        with open_stream(out_dir / 'launch.sh', 'wt', compressed=False) as stream:
            for env_variable in ['POKE_FONT', 'POKE_FIGSIZE', 'POKE_COLOR', 'POKE_QUANTILES']:
                if env_variable in os.environ:
                    stream.write(f'export {env_variable}="{os.environ[env_variable]}"\n')
            stream.write('poke ' + shlex.join(sys.argv[1:]) + '\n')



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

def write_atoms_and_pockets_and_density(args, structures: list[Path], out_dir: Path, atoms: AtomGroup, frames: list[int], pockets: list[list[Pocket]], scouted_output_indexes: list[int], extra_output_indexes: list[int], reference: AtomGroup):

    if args.density:
        try: 
            from gridData import Grid 
        except ImportError:
            raise UserError('The gridData package is required to compute the density map. Install it with "pip install gridDataFormats"')
        grid_margin = 5.
        grid_spacing = 1.
        bbox = reference.atoms.bbox()
        bbox = np.array([[x-grid_margin for x in bbox[0]], [x+grid_margin for x in bbox[1]]], dtype=np.float32)
        bbox_start = bbox[0] # - grid_margin
        shape = np.ceil((bbox[1] - bbox[0]) / grid_spacing).astype(np.int32)
        grid = np.zeros(shape, dtype=np.float32)

    if args.align is not None:
        reference_pos = reference.select_atoms(args.align).positions
        indexes = atoms.universe.select_atoms(args.align).ix
        # Ensure that the reference and the atoms have the same number of atoms
        # Better matching of the structures can be implemented in the future
        reference_pos = reference_pos[:min(len(reference_pos), len(indexes))]
        indexes = indexes[:min(len(reference_pos), len(indexes))]
        
        aligner = Aligner(reference=reference_pos, indexes=indexes)

    try:
        atoms_writer = None
        pockets_writer = None

        if args.extract or args.sort or args.align:
            atoms_writer = Writer(str(out_dir / 'frames.xtc'), multiframe=True)
            ts = atoms.universe.trajectory[0]

            if args.align:
                atoms.positions = aligner.compute(ts.positions).apply(atoms.positions)

            atoms.write(str(out_dir / 'structure.pdb'))
        
        if args.reference is not None:
            reference.atoms.write(str(out_dir / 'reference.pdb'))
        
        if args.volumes:
            pockets_writer = PocketWriter(open_stream(out_dir / 'volumes.pdb', 'wt', args.compress), resid_start=resid_start)

        for frame, pks in zip(frames, pockets):
            ts = atoms.universe.trajectory[frame]

            if args.align:
                aligner.compute(ts.positions)
                
            if atoms_writer:
                if args.align:
                    atoms.positions = aligner.apply(atoms.positions)
                atoms_writer.write(atoms)

            if pockets_writer or args.density:
                if args.align:
                    for p in pks:
                        p.translate(-aligner.positions_cog)
                        p.apply_matrix(aligner.rotation)
                        p.translate(aligner.reference_cog)
                if pockets_writer:
                    pockets_writer.write_frame(pks)
                if args.density:
                    nout = 0
                    for p in pks:
                        vp = p.volume_positions
                        add_pocket_to_density(grid, vp, bbox_start, grid_spacing, nout)
        
        if args.density:
            print(f'Maximum density value: {grid.max()}')
            # print(f'Minimum density value: {grid.min()}')
            # print(f'Minimum density value different from zero: {grid[grid != 0].min()}')
            print(f'Number of volume points out of the grid: {nout}')
            mdgrid = Grid(grid, origin=bbox_start+(grid_spacing*0.5), delta=np.repeat(grid_spacing, 3))
            mdgrid.export(str(out_dir / 'density.dx'), type='double') # type double maximizes pymol compatibility
            try:
                import mrcfile
                with mrcfile.new(str(out_dir / 'density.mrc'), overwrite=True) as mrc:
                    # Set data and header attributes
                    mrc.set_data(grid.T)
                    # Set the grid spacing along the 3 axes
                    mrc.voxel_size = grid_spacing
                    # Set the origin of the grid
                    mrc.header['origin']['x'] = bbox_start[0] + grid_spacing * 0.5
                    mrc.header['origin']['y'] = bbox_start[1] + grid_spacing * 0.5
                    mrc.header['origin']['z'] = bbox_start[2] + grid_spacing * 0.5
                    mrc.update_header_from_data()
                    mrc.update_header_stats()
            except:
                print('Could not save density map in mrc format. Only dx format was saved.')

        if args.pymol:
            if args.extract or args.sort or args.align:
                structures_script = [out_dir / 'structure.pdb', out_dir / 'frames.xtc']
            else:
                structures_script = structures
            additional_strucutres = {} if args.reference is None else {'reference': (out_dir / 'reference.pdb', 'brightorange', False)}
            pockets_writer.write_pymol_script(out_dir / 'pymol.pml', structures_script, additional_strucutres=additional_strucutres, scouted_output_indexes=scouted_output_indexes, extra_output_indexes=extra_output_indexes)

    finally:
        if atoms_writer:
            atoms_writer.close()
        if pockets_writer:
            pockets_writer.output.close()


def write_atoms_and_pockets(args, structures: list[Path], out_dir: Path, atoms: AtomGroup, frames: list[int], pockets: list[list[Pocket]], scouted_output_indexes: list[int], extra_output_indexes: list[int], reference: AtomGroup):
    
    if args.align is not None:
        reference_pos = reference.select_atoms(args.align).positions
        indexes = atoms.universe.select_atoms(args.align).ix
        # Ensure that the reference and the atoms have the same number of atoms
        # Better matching of the structures can be implemented in the future
        reference_pos = reference_pos[:min(len(reference_pos), len(indexes))]
        indexes = indexes[:min(len(reference_pos), len(indexes))]
        
        aligner = Aligner(reference=reference_pos, indexes=indexes)

    try:
        atoms_writer = None
        pockets_writer = None

        if args.extract or args.sort or args.align:
            atoms_writer = Writer(str(out_dir / 'frames.xtc'), multiframe=True)
            ts = atoms.universe.trajectory[0]

            if args.align:
                atoms.positions = aligner.compute(ts.positions).apply(atoms.positions)

            atoms.write(str(out_dir / 'structure.pdb'))

        if args.reference is not None:
            reference.atoms.write(str(out_dir / 'reference.pdb'))
        
        if args.volumes:
            pockets_writer = PocketWriter(open_stream(out_dir / 'volumes.pdb', 'wt', args.compress), resid_start=resid_start)
            
        for frame, pks in zip(frames, pockets):
            ts = atoms.universe.trajectory[frame]

            if args.align:
                aligner.compute(ts.positions)
                
            if atoms_writer:
                if args.align:
                    atoms.positions = aligner.apply(atoms.positions)
                atoms_writer.write(atoms)

            if pockets_writer:
                if args.align:
                    for p in pks:
                        p.translate(-aligner.positions_cog)
                        p.apply_matrix(aligner.rotation)
                        p.translate(aligner.reference_cog)

                pockets_writer.write_frame(pks)

        if args.pymol:
            if args.extract or args.sort or args.align:
                structures_script = [out_dir / 'structure.pdb', out_dir / 'frames.xtc']
            else:
                structures_script = structures

            additional_strucutres = {} if args.reference is None else {'reference': (out_dir / 'reference.pdb', 'brightorange', False)}
            pockets_writer.write_pymol_script(out_dir / 'pymol.pml', structures_script, additional_strucutres=additional_strucutres, scouted_output_indexes=scouted_output_indexes, extra_output_indexes=extra_output_indexes)

    finally:
        if atoms_writer:
            atoms_writer.close()
        if pockets_writer:
            pockets_writer.output.close()

def write_descriptors(out_dir, atoms: AtomGroup, pockets: list[Pocket], frames: list[int], write_atoms=False):
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
                    header = [ 
                        'OutputFrame', 
                        'OriginalFrame', 
                        'PocketResid', 
                        *props.keys() ,
                        'Residues'
                    ]
                    
                    if write_atoms:
                        header.append('Atoms')

                    writer.writerow(header)
                    header_written = True

                if write_atoms:
                    writer.writerow((
                        frame, 
                        frames[frame], 
                        id,  
                        *(f'{v:.2g}' for v in props.values()),
                        ' '.join(map(str, atoms[p.atoms].residues.resids)),
                        ' '.join(map(str, atoms[p.atoms].ids))
                    ))

                else:
                    writer.writerow((
                        frame, 
                        frames[frame], 
                        id,  
                        *(f'{v:.2g}' for v in props.values()),
                        ' '.join(map(str, atoms[p.atoms].residues.resids)),
                    ))

def plot_descriptors(out_dir, pockets: list[Pocket]):
    '''
        Plot the descriptors for each pocket in the list.
    '''
    import matplotlib.pyplot as plt

    # Check for POKE_FONT environment variable
    if 'POKE_FONT' in os.environ:
        plt.rcParams['font.family'] = os.environ['POKE_FONT']
    
    # Check for POKE_FIGSIZE environment variable. Deafult are 6.4 x 4.8 
    if 'POKE_FIGSIZE' in os.environ:
        plt.rcParams['figure.figsize'] = tuple(float(s.strip()) for s in os.environ['POKE_FIGSIZE'].split(','))

    color = os.environ.get('POKE_COLOR', None)
    do_quantiles = os.environ.get('POKE_QUANTILES', 'yes').lower() not in ('no', 'false', '0')    
        
    (out_dir / 'histograms').mkdir(parents=True, exist_ok=True)

    descriptors = np.array([list(p.props().values()) for pks in pockets for p in pks])
    names = Pocket.properties.keys()

    for name, values in zip(names, descriptors.T):
        fig, ax = plt.subplots()
        kwargs = Pocket.properties[name][1]
        # If an arg is a callable, call it passing values
        kwargs = {k: v(values) if callable(v) else v for k, v in kwargs.items()}
        # Update the color with the one from the environment
        if color: kwargs['color'] = color
        # Set default arguments
        kwargs.setdefault('bins', 'auto')
        
        ax.hist(values, **kwargs)
        if do_quantiles:
            q25, q50, q75 = np.quantile(values, [0.25, 0.5, 0.75])
            ax.axvline(q25, color='black', ls='--', label=f'Q25: {q25:.2g}')
            ax.axvline(q50, color='black', ls='-', label=f'Q50: {q50:.2g}')
            ax.axvline(q75, color='black', ls='--', label=f'Q75: {q75:.2g}')
            ax.legend()
        ax.set_xlabel(name)
        ax.set_ylabel('Counts')
        fig.tight_layout()
        fig.savefig(out_dir / 'histograms' / f'{name}.png', dpi=300)
        plt.close(fig)

def write_pickle(args, out_dir, pockets: list[list[Pocket]]):
    with open_stream(out_dir / 'pockets.pickle', 'wb', args.compress) as stream:
        pickle.dump(pockets, stream)


def open_stream(path: os.PathLike,  mode: str, compressed: bool):
    if compressed:
        return gzip.open(str(path) + '.gz', mode)
    else:
        return open(path, mode)

class FrameSpeedColumn(progress.ProgressColumn):
    def render(self, task: progress.Task) -> Text:
        speed = task.finished_speed or task.speed
        if speed is None:
            return Text("?", style="progress.data.speed")
        return Text(f"{speed:.0f} frames/s", style="progress.data.speed")
    
class PocketsFoundColumn(progress.ProgressColumn):
    def render(self, task: progress.Task) -> Text:
        n_filtered = task.fields['n_pockets'] - task.fields['n_kept_pockets']
 
        text = Text(f" {task.fields['n_kept_pockets']} pockets", style="magenta")

        if n_filtered > 0:
            text += Text(f" {n_filtered} discarted", style="cyan")
        return text


class Pool:
    def __init__(self, scouter: Poke, n_frames: int, n_cores: int=None) -> None:
        self.scouter = scouter

        # Set number of cores, if not provided by user
        # Criteria: 
        # - if trajectory is short, ensure at least 8 frames for each core
        # - for longer trajectories, ensure the number of cores does not exceed available cores,
        #   up to a maximum of 32
        if n_cores is None:
            n_cores = max(1, min(32, os.cpu_count(), n_frames // 8))
        self.n_cores = n_cores

        if n_cores > 1:
            self.pool = multiprocessing.Pool(n_cores)

            # Set chunksize accordingly 
            # for low number of cores, chunksize is smaller to update progress bar more frequently
            self.chunksize = max(8, min(n_frames // n_cores, n_cores*2))
        else:
            self.pool = None
            self.chunksize = 1

    def __enter__(self):
        if self.pool:
            self.pool.__enter__()
        return self
        
    def __exit__(self, exc_type, exc_value, traceback):
        if self.pool:
            self.pool.__exit__(exc_type, exc_value, traceback)

    def get_pockets(self, frames):
        if self.pool:
            # Run the first frame in the main process to avoid 
            # the numba compilation in the forked processes 
            yield self.scouter.get_pockets(next(frames))
            yield from self.pool.imap(self.scouter.get_pockets, frames, chunksize=self.chunksize)
        else:
            yield from map(self.scouter.get_pockets, frames)
             
class UserError(RuntimeError):
    pass
