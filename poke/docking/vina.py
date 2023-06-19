import io
import warnings
import MDAnalysis

from tempfile import NamedTemporaryFile

from MDAnalysis import AtomGroup
from vina import Vina
from openbabel import pybel

from .molecule import Molecule
from .docker import Docker


class AutodockVina(Docker):

    def __init__(self, protein: AtomGroup, box_center: list[float], box_size: list[float], 
                        sf_name='vinardo', cpus: int=None, verbosity=False, exhaustiveness: int=24, n_poses: int=4):
        self.exhaustiveness = exhaustiveness
        self.n_poses = n_poses

        self.vina = Vina(sf_name=sf_name, cpu=cpus or 0, verbosity=int(verbosity))
        self.load_receptor(protein, box_center=box_center, box_size=box_size)

        
    def load_receptor(self, protein: AtomGroup, box_center: list[float], box_size: list[float]):

        # Initialize a fake text file
        string_io = io.StringIO()

        # Write the atoms group to the string_io    
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)    
            with MDAnalysis.Writer(MDAnalysis.lib.util.NamedStream(string_io, filename=f'name.gro')) as writer:
                writer.write(protein)
            
        pdbqt = pybel.readstring('gro', string_io.getvalue()).write(format="pdbqt", opt={"r": "r"})

        # Pass the receptor pdb to vina 
        with NamedTemporaryFile(mode='w', suffix='.pdbqt') as receptor_file:

            # Write the pdb to the temporary file
            receptor_file.write(pdbqt)
            receptor_file.flush()

            # Load the pdb into vina
            self.vina.set_receptor(receptor_file.name)
            self.vina.compute_vina_maps(center=list(map(float, box_center)), box_size=list(map(float, box_size)))

    
    def dock(self, molecule: Molecule) -> Molecule:
        molecule.add_hydrogens().generate_3d_coordinates()

        self.vina.set_ligand_from_string(molecule.pdbqt)
        self.vina.dock(exhaustiveness=self.exhaustiveness, n_poses=self.n_poses)

        return self.vina.poses()