"""Poke scouter for protein druggable pockets"""

__version__ = "0.0.1"

import os
import numpy as np

from .pocket import Pocket
from .algorithms import compute_low_resolution_grid, enhance_volumes, \
                        intersect_box, is_reachable, get_region_atoms

class Poke:
    def __init__(self, spacing: float=0.5, cutoff: float=7.5, 
                 electro_force: float=2000, probe_sigma: float=3.83,
                 energy_threshold: float=-15, min_body_volume: float=25, 
                 min_spot_volume: float=2, max_spot_distance: float=8):
        '''
        Poke object found drugguable pockets.
        Atom properties must be provided using one of the methods available.
        Use the method .get_pockets() to search for pockets.
        The same instance can be reused multiple times with different positions.

        Parameters
        ------------
            spacing: float
                distance between points in the grid.
                A smaller distance means more precision but 
                longer computational times.
            cutoff: float
                The cutoff radius for atoms interactions.
                A bigger value means more precision but longer 
                computational times.
            electro_force: float
                Importance of the electrostatic potential in defining 
                the shape of the pocket. With higher values 
                regions with high electrostatic potential are favored 
                and hydrophilic regions on the surface can appear as pockets.
            probe_sigma: float
                Size of the probe used to define the shape of the pocket.
            energy_threshold: float
                After computing the energy map regions below this threshold
                are kept as pockets. Putting lower values the resulting
                pockets will be smaller but with higher affinity.
            min_body_volume: float
                A region under the threshold can be used as a pocket only
                if it has a volume greater of this value.
                Using lower values more small pockets will apper.
            min_spot_volume: float
                Small low-energy regions (spots) around a pocket are included 
                in the output pocket. Only regions greater min_spot_volume
                are included in the pocket.
            max_spot_distance: float
                Only small regions closer than this value are included in
                the output. 
        '''
        self.spacing = spacing
        self.inverse_spacing = 1 / spacing
        self.cutoff = cutoff
        self.grid_cutoff = cutoff * self.inverse_spacing
        self.electro_force = electro_force
        self.probe_sigma = probe_sigma
        self.energy_threshold = energy_threshold
        self.min_body_points = min_body_volume / spacing**3
        self.min_spot_points = min_spot_volume / spacing**3
        self.max_spot_distance = max_spot_distance / spacing

        self.charges: np.ndarray = None
        self.sigmas: np.ndarray = None
        self.epsilons: np.ndarray = None

        self.probe_sigmas2: np.ndarray = None
        self.probe_epsilons: np.ndarray = None
        self.sigma_charges: np.ndarray = None

        self.origin: np.ndarray = None
        self.shape = (0,0,0)
        self.shape_low = (0,0,0)

        self.grid_positions: np.ndarray = None
        self.neighbours: np.ndarray = None
        self.energy_map: np.ndarray = None

        self.parametrizer = None
        
    def get_pockets(self, positions: np.ndarray) -> list[Pocket]:
        '''
        Compute the pockets for the conformation specified.
        Before calling this method atom properties must be loaded
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
        # Compute origin and shape of the map
        pos_min = np.min(positions, axis=0)
        pos_max = np.max(positions, axis=0)

        # Define the origin and shape of the grid
        self.origin = pos_min
        self.shape = np.round((pos_max - pos_min) * self.inverse_spacing).astype(np.int32) + 1
        self.shape_low = self.shape // 4 + 2

        # Check the size of the box
        if any(l > 1000000 for l in self.shape):
            raise RuntimeError('The box surrounding the system is too big.')

        # Convert the positions of the atoms in grid indexes (floating point)
        self.grid_positions = (positions - self.origin) * self.inverse_spacing

        # Allocate arrays on the low-resolution grid
        self.neighbours = np.zeros((*self.shape_low, 512), dtype=np.int32)

        # Allocate the map for the bulk energy
        self.energy_map = np.full((*self.shape, 2), 100, dtype=np.float32)

        # Compute the bulk energy on the low-resolution grid
        compute_low_resolution_grid(self.energy_map, self.grid_cutoff, self.neighbours,
                     self.grid_positions, self.probe_sigmas2, 
                     self.probe_epsilons, self.sigma_charges)

        # Increase the resolution of the low-energy regions 
        _, pockets_points, pockets_slicer, surfaces_points, surfaces_slicer = enhance_volumes(
            self.energy_map, self.energy_threshold, 
            self.neighbours, self.grid_positions,
            self.probe_sigmas2, self.probe_epsilons, self.sigma_charges)

        # Converts the regions into objects
        regions = [Region(pockets_points[pockets_slicer[i]:pockets_slicer[i+1]],
                          surfaces_points[surfaces_slicer[i]:surfaces_slicer[i+1]]) 
                    for i in range(pockets_slicer.shape[0]-1)
                        if pockets_slicer[i+1]-pockets_slicer[i] >= self.min_spot_points]

        # Cluster the spots with the big regions
        pockets = self.cluster_regions(regions)

        # Create a grid with atoms types
        return pockets

    def cluster_regions(self, regions: list['Region']) -> list['Pocket']:
        '''
            Divide the regions in big regions (bodies) and small regions (spots).
            The spots that are near a body will be added to it.
            For each body + near spots a Pocket object will be returned.
        '''
        big_regions: list[Region] = []
        small_regions: list[Region] = []

        # Divide the regions in big regions (bodies) and small regions (spots)
        for r in regions:
            (big_regions if r.volume >= self.min_body_points else small_regions).append(r)

        pockets = []      
        for br in big_regions:
            volume_points = [br.volume_points]
            surface_points = [br.surface_points]

            # Search for small regions near the current big region
            for sr in small_regions:
                if intersect_box(br.box, sr.box, self.max_spot_distance):
                    if is_reachable(br.surface_points, sr.surface_points, self.max_spot_distance, 
                                    self.grid_positions, self.probe_sigmas2,
                                    self.neighbours):
                        volume_points.append(sr.volume_points)
                        surface_points.append(sr.surface_points)

            # Compute descriptors and add it to the current list of pockets
            pockets.append(self.compute_descriptors(volume_points, surface_points))
        return pockets

    def compute_descriptors(self, volume_points: np.ndarray, surface_points: np.ndarray) -> Pocket:
        '''
            Compute descriptors of a region to build a Pocket.
            The majority of the properties are computed using the surface_points
            because it is more efficent.

            Parameters
            ------------
                volume_points: np.ndarray
                    Points of the regions to compute regions
                surface_points: np.ndarray
                    Points on the surface of the same region
                    
        '''
        n_regions = len(volume_points)

        # Concatenate the list of points
        volume_points = np.concatenate(volume_points)
        surface_points = np.concatenate(surface_points)
        
        # Compute the list of the atoms surrounding the pocket
        atoms = get_region_atoms(surface_points, self.neighbours, self.grid_positions)

        # The enclousure is computed using the ratio between 
        # the number of atoms and the number of surface points
        enclousure = 1 - np.exp(-3 * atoms.shape[0] / surface_points.shape[0])

        # Compute three orthogonal lengths of the pocket diagonalizing the autocorrelation matrix
        residuals = surface_points - np.mean(surface_points, axis=0)
        variances = np.linalg.eigvalsh(np.dot(residuals.T, residuals)) / residuals.shape[0]
        lengths = np.sqrt(variances)[::-1] * self.spacing

        # Extract the electrostatic potential
        electro = self.energy_map[volume_points[:,0], volume_points[:,1], volume_points[:,2], 1]
        electro_surface = self.energy_map[surface_points[:,0], surface_points[:,1], surface_points[:,2], 1]

        # Compute the multipole expansion for the electrostatic potential in the points of the surface
        monopole = np.sum(electro_surface)
        dipole = np.linalg.norm(np.sum(electro_surface * residuals.T, axis=1))
        # quadrupoles = np.linalg.eigvalsh(np.dot(electro_surface * residuals.T, residuals)) / residuals.shape[0]
        # quadrupoles = np.sqrt(np.abs(quadrupoles) + 1e-8)[::-1] * self.spacing        

        # Counts polar and apolar atoms
        abs_charges = np.abs(self.charges[atoms])
        n_apolar_atoms, n_polar_atoms, n_charges_atoms = np.histogram(abs_charges, bins=[0, 0.3, 0.7, 10])[0]
        hydrophobicity = n_apolar_atoms / atoms.shape[0]

        # The pocket requires a spacing as a full-matrix 
        spacing_matrix = np.diag(np.full(3, self.spacing, dtype=np.float32))

        # Create a pocket with properties computed
        return Pocket(self.origin, spacing_matrix, volume_points=volume_points, 
                      surface_points=surface_points, electro=electro, 
                      monopole=monopole, dipole=dipole, n_regions=int(n_regions), 
                      atoms=atoms, lengths=lengths, enclousure=enclousure,
                      hydrophobicity=hydrophobicity, n_polar_atoms=n_polar_atoms, 
                      n_charges_atoms=n_charges_atoms)


    def write_bulk_map(self, out: os.PathLike):
        from gridData import Grid
        Grid(grid=np.clip(self.bulk_map.astype(np.float16), -25, 3), 
             origin=self.origin, 
             delta=[self.spacing]*3).export(out, file_format='dx')


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

