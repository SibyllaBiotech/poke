"""Poke scouter for protein druggable pockets"""

__version__ = "0.0.1"

import pickle
import numpy as np

from pathlib import Path
from .pocket import Pocket
from .parametrizer import AtomsParametrizer, charmm
from .algorithms import get_pockets

class Poke:
    scoring_models = {}

    def __init__(self, charges: np.ndarray, 
                 sigmas: np.ndarray, epsilons: np.ndarray,
                 spacing: float=0.5, cutoff: float=7.5, 
                 electro_force: float=1120, probe_sigma: float=3.8,
                 energy_threshold: float=-14.5, min_volume: float=15, 
                 min_spot_volume: float=2.7, max_bulk_distance: float=1.8, 
                 max_spot_distance: float=3.6, scoring_model: str='model'):
        '''
        Poke object found drugguable pockets.
        Atom properties must be provided using one of the methods available.
        Use the method .get_pockets() to search for pockets.
        The same instance can be reused multiple times with different positions.

        Parameters
        ------------
            charges: np.ndarray
                Array with the charges of the atoms.
            sigmas: np.ndarray
                Array with the Van der Waals radius of the atoms.
            epsilons: np.ndarray
                Array with the strength of the Van der Waals interactions.
            spacing: float
                Distance between points in the grid.
                A shorter distance corresponds to higher precision, 
                but it also results in longer computational time.
            cutoff: float
                Cutoff radius for atoms interactions.
                A bigger value corresponds to higher precision, 
                but it also results in longer computational time.
            electro_force: float
                Importance of the electrostatic potential in defining 
                the shape of the pocket. 
                With higher values, regions with high electrostatic potential are favored 
                and hydrophilic regions on the surface can appear as pockets.
            probe_sigma: float
                Size of the probe used to define the shape of the pocket.
            energy_threshold: float
                After computing the energy map, regions below this threshold
                are kept as pockets. 
                By selecting lower values, the resultant pockets 
                will be of reduced size yet possess stronger affinity.
            min_volume: float
                Minimum volume of a pocket in cubic angstroms.
                Pockets with a smaller volume are discarded.
            max_bulk_distance: float
                Regions within this distance are merged together to form the bulk in a
                recursive manner. Increasing this value less pockets are found, but they are bigger.
            max_spot_distance: float
                Regions within this distance relative to the bulk are merged in the pocket.
                The assignment is not recoursive.
        '''

        # Store the parameters
        self.spacing = spacing
        self.inverse_spacing = 1 / spacing
        self.cutoff = cutoff
        self.grid_cutoff = cutoff * self.inverse_spacing
        self.electro_force = electro_force
        self.probe_sigma = probe_sigma
        self.energy_threshold = energy_threshold
        self.min_points = min_volume * self.inverse_spacing**3
        self.min_spot_points = min_spot_volume * self.inverse_spacing**3
        self.max_bulk_distance = max_bulk_distance * self.inverse_spacing
        self.max_spot_distance = max_spot_distance * self.inverse_spacing

        # Store the atoms' properties
        self.charges = charges
        self.sigmas = sigmas
        self.epsilons = epsilons

        # Compute the combined sigmas and epsilons
        # using the probe sigma and the sigma and epsilon of the atoms
        combined_sigmas = 0.5 * (sigmas + self.probe_sigma)
        self.probe_sigmas2 = (combined_sigmas * self.inverse_spacing)**2 
        self.probe_epsilons = 6.26 * np.sqrt(epsilons)
        self.sigma_charges = self.electro_force * self.charges / combined_sigmas**6

        # Load the scoring model if it is not already loaded
        # The scoring model is a machine learning model that predicts the druggability of a pocket
        if scoring_model not in Poke.scoring_models:
            Poke.scoring_models[scoring_model] = pickle.loads((Path(__file__).parent / f'{scoring_model}.pickle').read_bytes())
        self.scoring_model = Poke.scoring_models[scoring_model]

    @classmethod
    def from_mdanalysis(cls, atoms, min_volume: float=15, scoring_model: str='model'):
        '''
        Load atoms' properties from an AtomGroup of MDAnalysis.
        Retrieve the names of the residues and atoms, and internally 
        it uses charmm to obtain charges, sigmas and epsilons. 
    
        Parameters
        ------------
            atoms: MDAnalysis.AtomGroup
                The atom group used to set the parameters
        '''
        resnames = []
        atom_names = []
        for r in atoms.residues:
            resnames.append(r.resname)
            atom_names.append(r.atoms.names)
        
        return cls.from_names(resnames, atom_names, min_volume=min_volume, scoring_model=scoring_model)
       
    @classmethod
    def from_names(cls, resnames: list[str], atom_names: list[list[str]],
                   min_volume: float=15, scoring_model: str='model') -> 'Poke':
        '''
        Retrieves atoms' properties based on the names of residues and atoms.
        Internally it uses charmm forcefield to obtain charges, sigmas and epsilons. 

        Parameters
        ------------
            resnames: list[str]
                List with the residues' names
            atom_names: list[list[str]]
                List containing the names of atoms for each residue.
                atom_names and resnames must have the same length. 
        '''
        parametrizer = AtomsParametrizer(charmm)

        charges, sigmas, epsilons = np.concatenate([parametrizer(resname, names) 
                        for resname, names in zip(resnames, atom_names)], axis=0).T

        return cls(charges=charges, sigmas=sigmas, 
                   epsilons=epsilons, min_volume=min_volume, scoring_model=scoring_model)

        
    def get_pockets(self, positions: np.ndarray) -> list[Pocket]:
        '''
        Compute the pockets for the conformation specified.
        Before calling this method, atoms' properties must be loaded
        and the number of atoms must correspond to the number
        of positions provided to this method.    

        Parameters
        ------------
            positions: np.ndarray

        Return values
        ------------
            pockets: list[Pocket]
                the list of pockets
        '''

        # Call the numba function to compute the pockets
        pockets_raw = get_pockets(positions=positions, inverse_spacing=self.inverse_spacing,
                        grid_cutoff=self.grid_cutoff, energy_threshold=self.energy_threshold,
                        probe_sigmas2=self.probe_sigmas2, probe_epsilons=self.probe_epsilons,
                        sigma_charges=self.sigma_charges, min_spot_points=self.min_spot_points,
                        max_spot_distance=self.max_spot_distance, max_bulk_distance=self.max_bulk_distance, 
                        min_points=self.min_points)

        # Create the pockets objects
        pockets = [
            Pocket(
                origin=origin, 
                spacing=np.diag(np.full(3, self.spacing, dtype=np.float32)), 
                volume_points=volume_points, 
                surface_points=surface_points,
                regions=regions,
                energy=energy,
                electro=electro,
                atoms=atoms,
                charges=self.charges[atoms],
                lengths=lengths,
                enclosure=enclosure
            )
            for origin, volume_points, surface_points, regions, atoms, enclosure, lengths, energy, electro in pockets_raw
        ]

        # Compute the score for each pocket
        for p in pockets:
            p.score = self.scoring_model(np.array(p.features()))[0]

        # Remove pockets that are contained in other pockets if they have a lower score
        pockets = [p for p in pockets if all(p.score >= q.score for q in pockets if p.regions < q.regions)]
        return pockets