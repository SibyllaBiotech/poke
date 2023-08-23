"""Poke scouter for protein druggable pockets"""

__version__ = "0.0.1"

import pickle
import numpy as np

from pathlib import Path
from .pocket import Pocket
from .parametrizer import AtomsParametrizer, charmm
from .algorithms import compute_low_resolution_grid, enhance_volumes, \
                        intersect_box, is_reachable, get_region_atoms

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
            min_body_volume: float
                A region under the threshold can be used as a pocket only
                if it has a volume greater of this value.
                Using lower values will result in the generation of smaller pockets.
            min_spot_volume: float
                Small low-energy regions (spots) around a pocket, which are greater than min_spot_volume 
                are included in the output pocket. 
            max_spot_distance: float
                Only regions within a distance less than this value are included in
                the output. 
        '''
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

        self.charges = charges
        self.sigmas = sigmas
        self.epsilons = epsilons

        combined_sigmas = 0.5 * (sigmas + self.probe_sigma)

        self.probe_sigmas2 = (combined_sigmas * self.inverse_spacing)**2 
        self.probe_epsilons = 6.26 * np.sqrt(epsilons)
        self.sigma_charges = self.electro_force * self.charges / combined_sigmas**6

        if scoring_model not in Poke.scoring_models:
            Poke.scoring_models[scoring_model] = pickle.loads((Path(__file__).parent / f'{scoring_model}.pickle').read_bytes())
        self.scoring_model = Poke.scoring_models[scoring_model]

    @classmethod
    def from_mdanalysis(cls, atoms, min_volume: float=15):
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
        
        return cls.from_names(resnames, atom_names, min_volume=min_volume)
       
    @classmethod
    def from_names(cls, resnames: list[str], atom_names: list[list[str]],
                   min_volume: float=15) -> 'Poke':
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
                   epsilons=epsilons, min_volume=min_volume)

        
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
        # Compute the origin and the shape of the map
        pos_min = np.min(positions, axis=0)
        pos_max = np.max(positions, axis=0)

        # Define the origin and the shape of the grid
        origin = pos_min
        shape = np.round((pos_max - pos_min) * self.inverse_spacing).astype(np.int32) + 1
        shape_low = shape // 4 + 2

        # Check the size of the box
        if any(l > 1000000 for l in shape):
            raise RuntimeError('The box surrounding the system is too big.')

        # Convert the positions of the atoms in grid indexes (floating point)
        grid_positions = (positions - origin) * self.inverse_spacing

        # Allocate arrays on the low-resolution grid
        neighbours = np.zeros((*shape_low, 512), dtype=np.int32)

        # Allocate the map for the bulk energy
        energy_map = np.full((*shape, 2), 100, dtype=np.float32)

        # Compute the bulk energy on the low-resolution grid
        compute_low_resolution_grid(energy_map, self.grid_cutoff, neighbours,
                     grid_positions, self.probe_sigmas2, 
                     self.probe_epsilons, self.sigma_charges)

        # Increase the resolution of the low-energy regions 
        _, pockets_points, pockets_slicer, surfaces_points, surfaces_slicer = enhance_volumes(
            energy_map, self.energy_threshold, 
            neighbours, grid_positions,
            self.probe_sigmas2, self.probe_epsilons, self.sigma_charges)

        # Convert the regions into objects
        regions = [Region(pockets_points[pockets_slicer[i]:pockets_slicer[i+1]],
                          surfaces_points[surfaces_slicer[i]:surfaces_slicer[i+1]]) 
                    for i in range(pockets_slicer.shape[0]-1)
                        if pockets_slicer[i+1]-pockets_slicer[i] >= self.min_spot_points]

        # Cluster the spots with the big regions
        return [
            self.compute_descriptors(volume_points=volume_points, surface_points=surface_points,
                                     origin=origin, energy_map=energy_map, 
                                     grid_positions=grid_positions, neighbours=neighbours)
            for volume_points, surface_points in self.cluster_regions(regions, grid_positions, neighbours)
        ]

    def cluster_regions(self, regions: list['Region'], 
                        grid_positions: np.ndarray,
                        neighbours: np.ndarray) -> list['Pocket']:
        '''
            Divide the regions into big sections (bodies) and small sections (spots). 
            The spots that are in close proximity to a body will be included within it.
            For each (body + near spots) a Pocket object will be returned.
        '''
        # big_regions: set[Region] = set()
        # small_regions: set[Region] = set()
        regions1 = set(regions)
        regions2 = regions1.copy()
        grouped_regions: list[list[Region]] = []

        # Divide the regions in big regions (bodies) and small regions (spots)
        # for r in regions:
        #     (big_regions if r.volume >= self.min_body_points else small_regions).add(r)
        
        # Merge big regions in groups if they are close enough
        while regions1:
            group = [regions1.pop()]
            stack = {group[0]}
            while stack:
                r1 = stack.pop()
                for r2 in list(regions1):
                    if intersect_box(r1.box, r2.box, self.max_bulk_distance):
                        if is_reachable(r1.surface_points, r2.surface_points, self.max_bulk_distance, 
                                    grid_positions, self.probe_sigmas2, neighbours):
                            regions1.remove(r2)
                            group.append(r2)
                            stack.add(r2)

            grouped_regions.append(group)

        for br_group in grouped_regions:
            volume_points = []
            surface_points = []
            regions2_tmp = regions2.copy()
            for br in br_group:
                regions2_tmp.remove(br)
            
            for br in br_group:
                volume_points.append(br.volume_points)
                surface_points.append(br.surface_points)
                for sr in list(regions2_tmp):
                    if intersect_box(br.box, sr.box, self.max_spot_distance):
                         if is_reachable(br.surface_points, sr.surface_points, self.max_spot_distance, 
                                         grid_positions, self.probe_sigmas2, neighbours):
                            volume_points.append(sr.volume_points)
                            surface_points.append(sr.surface_points)
                            regions2_tmp.remove(sr)
            
            volume_points = np.concatenate(volume_points)
            if volume_points.shape[0] >= self.min_points:
                surface_points = np.concatenate(surface_points)
                yield volume_points, surface_points

    def compute_descriptors(self, volume_points: np.ndarray, 
                                  surface_points: np.ndarray,
                                  origin: tuple[float, float, float],
                                  energy_map: np.ndarray,
                                  grid_positions: np.ndarray,
                                  neighbours: np.ndarray) -> Pocket:
        '''
            Compute descriptors of a region to build a Pocket.
            Most properties are computed using the surface_points
            due to their greater efficiency.

            Parameters
            ------------
                volume_points: np.ndarray
                    Points within the regions for region computation
                surface_points: np.ndarray
                    Points on the surface of the same region
                    
        '''
        # n_regions = len(volume_points)

        # Concatenate the list of points
        # volume_points = np.concatenate(volume_points)
        # surface_points = np.concatenate(surface_points)

        # Compute the list of the atoms surrounding the pocket
        atoms, net_vector = get_region_atoms(surface_points, neighbours, grid_positions, k_nearest=int(20 * self.spacing**2 + 0.5))

        # The enclosure is computed using the ratio between 
        # the number of atoms and the number of surface points
        enclosure = np.exp(-0.3 * net_vector**2)

        # Compute three orthogonal lengths of the pocket diagonalizing the autocorrelation matrix
        residuals = surface_points - np.mean(surface_points, axis=0)
        variances = np.linalg.eigvalsh(np.dot(residuals.T, residuals)) / residuals.shape[0]
        lengths = np.sqrt(variances)[::-1] * self.spacing

        # Extract the electrostatic potential
        energy = energy_map[volume_points[:,0], volume_points[:,1], volume_points[:,2], 0] 
        electro = energy_map[volume_points[:,0], volume_points[:,1], volume_points[:,2], 1]

        # Count polar and apolar atoms
        charges = self.charges[atoms]

        # The pocket requires a spacing as a full-matrix 
        spacing_matrix = np.diag(np.full(3, self.spacing, dtype=np.float32))

        # Create a pocket with the computed properties
        pocket = Pocket(origin, spacing_matrix, volume_points=volume_points, 
                      surface_points=surface_points, energy=energy, electro=electro, 
                      atoms=atoms, charges=charges, lengths=lengths, enclosure=enclosure)

        # Compute the score
        pocket.score = self.scoring_model(np.array(pocket.features()))[0]
        return pocket

class Region:
    def __init__(self, volume_points: np.ndarray, surface_points: np.ndarray):
        self.volume_points = volume_points
        self.surface_points = surface_points

        # Total volume of the region
        self.surface = self.surface_points.shape[0]
        self.volume = self.volume_points.shape[0]
        
        # Box containing all the points
        self.box = np.array([np.min(surface_points, axis=0), np.max(surface_points, axis=0)])

    def ix(self):
        ''' Useful to extract values from a map using map[region.ix()] '''
        return self.volume_points[:,0], self.volume_points[:,1], self.volume_points[:,2]

    def six(self):
        ''' Useful to extract values from a map using map[region.six()] '''
        return self.surface_points[:,0], self.surface_points[:,1], self.surface_points[:,2]