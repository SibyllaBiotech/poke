import numpy as np
import numba as nb

@nb.njit(cache=True, fastmath=True)
def compute_low_resolution_grid(energy_map: np.ndarray, grid_cutoff: float, 
                 neighbours: np.ndarray, grid_positions: np.ndarray, 
                 probe_sigmas2: np.ndarray, probe_epsilons: np.ndarray, sigma_charges: np.ndarray):
    """
        Compute the low-resolution energy_map and build the  neighbours lists.
        In the low-resolution energy_map only the values energy_map[::4,::4,::4] are computed.
        The remaining values are kept unchanged.

        For each point in the low-resolution grid a list of neighbour atoms is build in the
        neighbours ndarray. The neighbours array must have a quarter of the shapes of the grid.
        The first value in the list is the number of neighbours, the following are the indexes of the atoms.

        The sigma_charges are the charge values transformed to be used with distances in unit of measure of sigma.
    """

    # Set to zero the values of the map that will be computed
    energy_map[::4,::4,::4] = 0

    # Precompute the cutoff squared. It will be used in the hot-cycle 
    grid_cutoff2 = grid_cutoff * grid_cutoff

    # For each atom, its potential will be added in the neighbour cells.
    # At the same time the index of the atom will be added in the neighbours lists.
    for i in range(grid_positions.shape[0]):

        # Van de Walls composition rules between the atom and the probe
        charge = sigma_charges[i]
        sigma2 = probe_sigmas2[i]
        epsilon = probe_epsilons[i]

        # Floating point position on the grid
        nx = grid_positions[i,0]
        ny = grid_positions[i,1]
        nz = grid_positions[i,2]

        # Compute the influence boundaries of this atom
        min_kx = np.int32(max(nx-grid_cutoff, 0) * 0.25 + 0.5) << 2
        min_ky = np.int32(max(ny-grid_cutoff, 0) * 0.25 + 0.5) << 2
        min_kz = np.int32(max(nz-grid_cutoff, 0) * 0.25 + 0.5) << 2
        max_kx = np.int32(min(nx+grid_cutoff, energy_map.shape[0]) * 0.25 + 0.5) << 2
        max_ky = np.int32(min(ny+grid_cutoff, energy_map.shape[1]) * 0.25 + 0.5) << 2
        max_kz = np.int32(min(nz+grid_cutoff, energy_map.shape[2]) * 0.25 + 0.5) << 2

        # The step of 4 needs to match the low-resolution grid
        for kx in range(min_kx, max_kx, 4):
            kx4 = kx >> 2 
            dx2 = (kx-nx)*(kx-nx)

            for ky in range(min_ky, max_ky, 4):
                ky4 = ky >> 2
                dy2 = (ky-ny)*(ky-ny)
                dx2y2 = dx2 + dy2

                for kz in range(min_kz, max_kz, 4):
                    dz2 = (kz-nz)*(kz-nz)

                    # Distance squared between the atoms 
                    # and the point on the grid
                    r2 = dx2y2 + dz2 + 1e-2

                    # Avoid the computation if the atom
                    # is too far away
                    if r2 >= grid_cutoff2:
                        continue
                
                    kz4 = kz >> 2
                    ir2 = sigma2 / r2
                    ir6 = ir2*ir2*ir2
                    ir12 = ir6*ir6

                    # Add the atom in the neighbours list
                    # of this point in the grid                       
                    n_neighbours = neighbours[kx4,ky4,kz4,0]
                    i_neighbour = n_neighbours + 1
                    if i_neighbour < neighbours.shape[3]:
                        neighbours[kx4,ky4,kz4,0] = i_neighbour
                        neighbours[kx4,ky4,kz4,i_neighbour] = i

                    # Add the energy contribution
                    energy_map[kx,ky,kz,0] += epsilon * (ir12 - ir6)
                    energy_map[kx,ky,kz,1] += charge * ir6

    # Add the strength (absolute value) of the eletrostatic potential to the Van der Waals energy 
    energy_map[::4,::4,::4, 0] -= np.abs(energy_map[::4,::4,::4, 1])
