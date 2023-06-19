import io
import os
import warnings
import tempfile
import subprocess
import MDAnalysis
import numpy as np

from pathlib import Path
from MDAnalysis import AtomGroup
from openbabel import pybel

from .molecule import Molecule, SDFWriter, SDFReader
from .docker import Docker


class Gnina(Docker):

    def __init__(self, protein: AtomGroup, box_origin: list[np.ndarray], box_vectors: list[np.ndarray], padding=3,
                        executable='gnina', n_cpus: int=None, verbose=False, exhaustiveness: int=24, n_poses: int=4, env=os.environ.copy(), gpu=False):
        self.executable = executable
        self.protein = protein
        self.n_cpus = n_cpus
        self.box_origin = box_origin
        self.box_vectors = box_vectors
        self.exhaustiveness = exhaustiveness
        self.n_poses = n_poses
        self.verbose = verbose
        self.environment = env

        self.box_size = np.linalg.norm(box_vectors, axis=1) 
        self.transform = box_vectors.T / self.box_size
        self.box_center = (box_origin @ self.transform) + 0.5 * self.box_size

        self.box_size += padding

        if not gpu:
            env.setdefault("CUDA_VISIBLE_DEVICES", "-1")

    def get_receptor(self):

        # Save the original positions 
        original_positions = self.protein.positions
        rotated_positions = original_positions @ self.transform

        # Initialize a fake text file
        string_io = io.StringIO()

        # Write the atoms group to the string_io    
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)    
            with MDAnalysis.Writer(MDAnalysis.lib.util.NamedStream(string_io, filename=f'name.gro')) as writer:
                self.protein.positions = rotated_positions
                writer.write(self.protein)
                self.protein.positions = original_positions
            
        return pybel.readstring('gro', string_io.getvalue()).write(format="pdbqt", opt={"r": "r"})

    
    def dock_list(self, molecules: list[Molecule]):
        if len(molecules) == 0:
            return []

        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            receptor = directory / 'receptor.pdbqt'
            lig_input = directory / 'lig_input.sdf'
            lig_output = directory / 'lig_output.sdf'

            args = [
                self.executable,
                '--receptor', receptor.as_posix(),
                '--ligand', lig_input.as_posix(), 
                '--exhaustiveness', str(self.exhaustiveness),
                '--addH', 'off',
                '--stripH', 'off',
                '--num_modes', str(self.n_poses),
                '--cnn', 'general_default2018',

                '--center_x', str(self.box_center[0]),
                '--center_y', str(self.box_center[1]),
                '--center_z', str(self.box_center[2]),

                '--size_x', str(self.box_size[0]),
                '--size_y', str(self.box_size[1]),
                '--size_z', str(self.box_size[2]),

                '--out', lig_output.as_posix()
            ]

            if self.n_cpus is not None:
                args.extend(['--cpu', str(self.n_cpus)])

            stdout = subprocess.PIPE if not self.verbose else None
            stderr = subprocess.PIPE if not self.verbose else None

            receptor.write_text(self.get_receptor())

            output_molecules = []

            with SDFWriter(lig_input) as writer:
                for i, m in enumerate(molecules):
                    try:
                        m = m.add_hydrogens().generate_3d_coordinates().apply_matrix(self.transform.T)
                        m['gnina_batch_index'] = str(i)
                        writer.write(m)
                    finally:
                        output_molecules.append([])
                    
            result = subprocess.run(args, stdout=stdout, stderr=stderr, env=self.environment)

            if result.returncode:
                raise RuntimeError(result.stderr.decode("utf8") if stderr is not None else '')

            with SDFReader(lig_output, sanitize=False) as reader:
                docked = [m for m in reader] 

            for m in docked:
                m['Energy'] = '{:.3f}'.format(-1.3635606572452885 * float(m['CNNaffinity']))
                m.apply_matrix(self.transform)
                output_molecules[int(m['gnina_batch_index'])].append(m)
                del m['gnina_batch_index']

            return output_molecules