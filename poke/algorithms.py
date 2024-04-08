import numpy as np
import numba as nb

from .utils import compute_box, average_3d, fetch_3d, intersect_box, is_reachable

@nb.njit(cache=True, fastmath=True, nogil=True)
def compute_low_resolution_grid(energy_map: np.ndarray, grid_cutoff: float, 
                 neighbours: np.ndarray, grid_positions: np.ndarray, 
                 probe_sigmas2: np.ndarray, probe_epsilons: np.ndarray, sigma_charges: np.ndarray):
    """
        Compute the low-resolution energy_map and build the neighbours lists.
        In the low-resolution energy_map, only energy_map[::4,::4,::4] values 
        are calculated, while the rest remain unchanged.
        
        For each point in the low-resolution grid, a list of neighboring atoms is built in the
        neighbours ndarray. The neighbours array must have a quarter of the dimensions of the grid.
        The list begins with the number of neighbors, followed by the atoms' indexes. 

        The sigma_charges refer to the charge values which have been transformed to be used with 
        the distances measured in units of sigma.
    """

    # Set to zero the values of the map that will be computed
    energy_map[::4,::4,::4] = 0

    # Precompute the cutoff squared. It will be used in the hot-cycle 
    grid_cutoff2 = grid_cutoff * grid_cutoff

    # For each atom, its potential will be added to the neighbour cells.
    # At the same time the index of the atom will be added to the neighbours lists.
    for i in range(grid_positions.shape[0]):

        # Van der Waals composition rules between the atom and the probe
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

        # Use a value of 4 to match the step of the low-resolution grid,
        # as it's one-fourth the resolution of the original grid.
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
                    # for this grid point                       
                    n_neighbours = neighbours[kx4,ky4,kz4,0]
                    i_neighbour = n_neighbours + 1
                    if i_neighbour < neighbours.shape[3]:
                        neighbours[kx4,ky4,kz4,0] = i_neighbour
                        neighbours[kx4,ky4,kz4,i_neighbour] = i

                    # Add the energy contribution
                    energy_map[kx,ky,kz,0] += epsilon * (ir12 - ir6)
                    energy_map[kx,ky,kz,1] += charge * ir6

@nb.njit(cache=True, fastmath=True, nogil=True)
def enhance_volumes(energy_map: np.ndarray, isovalue: float, 
                    neighbours: np.ndarray, grid_positions: np.ndarray,
                    probe_sigmas2: np.ndarray, probe_epsilons: np.ndarray, sigma_charges: np.ndarray,
                    min_spot_points: int):
    '''
        Expand the low-energy regions found in the low-resolution map. This expansion is 
        carried out on the high-resolution map by computing the potential in the missing positions.
        This function returns the set of points that fall below a certain threshold of total energy.
    '''
    # To store the points that will be included in pockets_points
    inserted = np.zeros((energy_map.shape[0], energy_map.shape[1], energy_map.shape[2]), dtype=np.uint16)

    # Create lists of points, each corresponding to a specific pocket 
    pockets_points = np.empty((energy_map.size, 3), dtype=np.int16) # 65536 maximum dimension of the grid along each axis
    surfaces_points = np.empty((energy_map.size // 2, 3), dtype=np.int16) 

    # For each pocket, mark the starting index of the volume/surface points array (excluding the first one, which starts at zero).
    # These are used as slicers for the pocket_points and surface_points arrays, eliminating the need to computate the sum of the volumes array.
    # Example: the indices [17,69,666] can be used as pockets_points[17:69] to extract the points belonging to the second pocket.
    pockets_slicer  = np.zeros(65536, dtype=np.uint32)
    surfaces_slicer = np.zeros(65536, dtype=np.uint32)

    # Pocket counter, needed for ids
    # This is also used to keep track of the boundaries of volumes and surfaces arrays
    n_pockets = 0

    n_visited = 0
    n_inserted = 0  
    n_surfaces = 0

    for i in range(0, energy_map.shape[0], 4):
        for j in range(0, energy_map.shape[1], 4):
            for k in range(0, energy_map.shape[2], 4):

                # Find first point belonging to new pocket (root of graph)
                if energy_map[i,j,k,0] < isovalue and not inserted[i,j,k]:

                    if n_pockets+1 >= pockets_slicer.shape[0]:
                        break

                    # New pocket found
                    n_pockets += 1

                    inserted[i,j,k] = n_pockets
                    pockets_points[n_inserted] = (i,j,k)
                    n_inserted += 1                   
                    
                    while n_visited < n_inserted:

                        # Visit a new point
                        px = pockets_points[n_visited,0]
                        py = pockets_points[n_visited,1]
                        pz = pockets_points[n_visited,2]
                        n_visited += 1
                        
                        n_internal_faces = 0

                        # Compute the energies of nearby points 
                        # that will be used in the following if clauses
                        energy_xp, energy_xm, energy_yp, energy_ym, energy_zp, energy_zm = \
                            compute_star_points(px,py,pz, energy_map, 
                                            neighbours, grid_positions, 
                                            probe_sigmas2, probe_epsilons, sigma_charges)
                        
                        # Decide which points are below the threshold and
                        # must be added to the returning lists
                        if energy_xp < isovalue:
                            n_internal_faces += 1
                            if not inserted[px+1,py,pz]:
                                inserted[px+1,py,pz] = n_pockets
                                pockets_points[n_inserted] = (np.int16(px+1),py,pz)
                                n_inserted += 1

                        if energy_xm < isovalue:
                            n_internal_faces += 1
                            if not inserted[px-1,py,pz]:
                                inserted[px-1,py,pz] = n_pockets
                                pockets_points[n_inserted] = (np.int16(px-1),py,pz)
                                n_inserted += 1

                        if energy_yp < isovalue:
                            n_internal_faces += 1
                            if not inserted[px,py+1,pz]:
                                inserted[px,py+1,pz] = n_pockets
                                pockets_points[n_inserted] = (px,np.int16(py+1),pz)
                                n_inserted += 1

                        if energy_ym < isovalue:
                            n_internal_faces += 1
                            if not inserted[px,py-1,pz]:
                                inserted[px,py-1,pz] = n_pockets
                                pockets_points[n_inserted] = (px,np.int16(py-1),pz)
                                n_inserted += 1

                        if energy_zp < isovalue:
                            n_internal_faces += 1
                            if not inserted[px,py,pz+1]:
                                inserted[px,py,pz+1] = n_pockets
                                pockets_points[n_inserted] = (px,py,np.int16(pz+1))
                                n_inserted += 1

                        if energy_zm < isovalue:
                            n_internal_faces += 1
                            if not inserted[px,py,pz-1]:
                                inserted[px,py,pz-1] = n_pockets
                                pockets_points[n_inserted] = (px,py,np.int16(pz-1))
                                n_inserted += 1

                        if n_internal_faces <= 4:
                            surfaces_points[n_surfaces] = (px,py,pz)
                            n_surfaces += 1

                    pockets_slicer[n_pockets] = n_inserted
                    surfaces_slicer[n_pockets] = n_surfaces
        
    # Convert the regions into objects
    i_region: int = 0

    # List of tuples of the form (region_id, volume_points, surface_points, box, (atoms, net_vector, net_count))
    regions: list[tuple[int, np.ndarray, np.ndarray, np.ndarray, tuple[set[int], np.ndarray, int]], np.ndarray] = nb.typed.List()
    # The regions are separated by the pockets_slicer and surfaces_slicer arrays
    # which contain the indices of the first point of each region
    for i in range(n_pockets):
        if pockets_slicer[i+1]-pockets_slicer[i] >= min_spot_points:
            volume_points = pockets_points[pockets_slicer[i]:pockets_slicer[i+1]]
            surface_points = surfaces_points[surfaces_slicer[i]:surfaces_slicer[i+1]]

            # Build tuples of the form (region_id, volume_points, surface_points, box, (atoms, net_vector, net_count))
            regions.append((
                i_region,
                volume_points,
                surface_points,
                compute_box(surface_points),
                get_region_atoms(surface_points, neighbours=neighbours, grid_positions=grid_positions, k_nearest=5),
                compute_accessible_volume(surface_points, energy_map=energy_map, cutoff=2)
            ))

            i_region += 1

    return regions


@nb.njit(cache=True, fastmath=True, nogil=True)
def compute_star_points(px: int, py: int, pz: int, energy_map: np.ndarray,
                        neighbours: np.ndarray, grid_positions: np.ndarray, 
                        probe_sigmas2: np.ndarray, probe_epsilons: np.ndarray, sigma_charges: np.ndarray):
    '''
        Compute the potential of the 6 points around the central point (px, py, pz).
        The computation only takes into consideration the atoms present in the neighbours list.
    '''

    # Fetch energy values from the map
    energy_xp = energy_map[px+1,py,pz,0] if px+1 < energy_map.shape[0] else 1000
    energy_xm = energy_map[px-1,py,pz,0] if px > 0 else 1000

    energy_yp = energy_map[px,py+1,pz,0] if py+1 < energy_map.shape[1] else 1000
    energy_ym = energy_map[px,py-1,pz,0] if py > 0 else 1000

    energy_zp = energy_map[px,py,pz+1,0] if pz+1 < energy_map.shape[2] else 1000
    energy_zm = energy_map[px,py,pz-1,0] if pz > 0 else 1000

    # Check whether those values have not been previously calculated    
    compute_xp, energy_xp, elec_xp = (True, 0, 0) if energy_xp == 100 else (False, energy_xp, 0) 
    compute_xm, energy_xm, elec_xm = (True, 0, 0) if energy_xm == 100 else (False, energy_xm, 0) 
    compute_yp, energy_yp, elec_yp = (True, 0, 0) if energy_yp == 100 else (False, energy_yp, 0) 
    compute_ym, energy_ym, elec_ym = (True, 0, 0) if energy_ym == 100 else (False, energy_ym, 0) 
    compute_zp, energy_zp, elec_zp = (True, 0, 0) if energy_zp == 100 else (False, energy_zp, 0) 
    compute_zm, energy_zm, elec_zm = (True, 0, 0) if energy_zm == 100 else (False, energy_zm, 0) 

    kx4 = (px + 2) >> 2
    ky4 = (py + 2) >> 2
    kz4 = (pz + 2) >> 2
    n_neighbours = neighbours[kx4,ky4,kz4,0]
    i_neighbours = neighbours[kx4,ky4,kz4,1:]

    for i_neighbour in range(n_neighbours):
        i = i_neighbours[i_neighbour]

        # Van der Waals composition rules between the atom and the probe
        charge = sigma_charges[i]
        sigma2 = probe_sigmas2[i]
        epsilon = probe_epsilons[i]

        # Floating point position on the grid
        nx = grid_positions[i,0]
        ny = grid_positions[i,1]
        nz = grid_positions[i,2]

        px_nx = px-nx
        py_ny = py-ny
        pz_nz = pz-nz

        r2 = px_nx*px_nx + py_ny*py_ny + pz_nz*pz_nz + 1.002

        # Add the energy contribution of this atom 
        # to the 6 points around
        if compute_xp:
            ir2 = sigma2 / (r2 + 2 * px_nx)
            ir6 = ir2*ir2*ir2
            ir12 = ir6*ir6
            energy_xp += epsilon * (ir12 - ir6)
            elec_xp += charge * ir6

        if compute_xm:
            ir2 = sigma2 / (r2 - 2 * px_nx)
            ir6 = ir2*ir2*ir2
            ir12 = ir6*ir6
            energy_xm += epsilon * (ir12 - ir6)
            elec_xm += charge * ir6

        if compute_yp:
            ir2 = sigma2 / (r2 + 2 * py_ny)
            ir6 = ir2*ir2*ir2
            ir12 = ir6*ir6
            energy_yp += epsilon * (ir12 - ir6)
            elec_yp += charge * ir6

        if compute_ym:
            ir2 = sigma2 / (r2 - 2 * py_ny)
            ir6 = ir2*ir2*ir2
            ir12 = ir6*ir6
            energy_ym += epsilon * (ir12 - ir6)
            elec_ym += charge * ir6

        if compute_zp:
            ir2 = sigma2 / (r2 + 2 * pz_nz)
            ir6 = ir2*ir2*ir2
            ir12 = ir6*ir6
            energy_zp += epsilon * (ir12 - ir6)
            elec_zp += charge * ir6

        if compute_zm:
            ir2 = sigma2 / (r2 - 2 * pz_nz)
            ir6 = ir2*ir2*ir2
            ir12 = ir6*ir6
            energy_zm += epsilon * (ir12 - ir6)
            elec_zm += charge * ir6


    # Save the computed values in the map
    if compute_xp:
        energy_map[px+1,py,pz,0] = energy_xp
        energy_map[px+1,py,pz,1] = elec_xp
    if compute_xm:
        energy_map[px-1,py,pz,0] = energy_xm
        energy_map[px-1,py,pz,1] = elec_xm
    if compute_yp:
        energy_map[px,py+1,pz,0] = energy_yp
        energy_map[px,py+1,pz,1] = elec_yp
    if compute_ym:
        energy_map[px,py-1,pz,0] = energy_ym
        energy_map[px,py-1,pz,1] = elec_ym
    if compute_zp:
        energy_map[px,py,pz+1,0] = energy_zp
        energy_map[px,py,pz+1,1] = elec_zp
    if compute_zm:
        energy_map[px,py,pz-1,0] = energy_zm
        energy_map[px,py,pz-1,1] = elec_zm

    return energy_xp, energy_xm, energy_yp, energy_ym, energy_zp, energy_zm
    
@nb.njit(cache=True, fastmath=True, nogil=True)   
def get_region_atoms(points, neighbours, grid_positions, k_nearest):
    '''
        This Function identifies the k-nearest atoms for every point within the region.
        After collecting all the atoms for each point, it returns the unique sorted list.
    '''
    atoms = set()

    # Vector used to accumulate data to compute the enclosure
    net_vector = np.zeros(shape=3, dtype=np.float32)
    net_count = 0

    # These arrays are used multiple times to store the k-nearest atoms
    # for each atom. Atoms collection is achieved through an insertion sort 
    # using these arrays
    k_r2 = np.empty(k_nearest, dtype=np.float32)
    k_i  = np.empty(k_nearest, dtype=np.int32)

    for p in range(points.shape[0]):
        px = points[p,0]
        py = points[p,1]
        pz = points[p,2]

        kx4 = (px + 2) >> 2
        ky4 = (py + 2) >> 2
        kz4 = (pz + 2) >> 2

        n_neighbours = neighbours[kx4,ky4,kz4,0]
        i_neighbours = neighbours[kx4,ky4,kz4,1:]

        # Reset the k-nearest arrays
        k_r2[:] = 1000
        k_i[:] = -1

        for i_neighbour in range(n_neighbours):
            i = i_neighbours[i_neighbour]

            # Floating point position on the grid
            nx = grid_positions[i,0]
            ny = grid_positions[i,1]
            nz = grid_positions[i,2]

            r2 = (px-nx)*(px-nx) + (py-ny)*(py-ny) + (pz-nz)*(pz-nz)

            # Insertion sort
            k = 0
            while r2 < k_r2[k] and k < k_nearest:
                k_r2[k] = k_r2[k+1]
                k_i[k] = k_i[k+1]
                k += 1

            if k > 0:
                k_i[k-1] = i
                k_r2[k-1] = r2

        # Transfer the k-nearest arrays into the set of atoms
        for i in range(k_i.shape[0]):
            if k_i[i] != -1:
                atoms.add(k_i[i])

                # Update the accumulation vector with the 
                # vector from the pocket point to the atom              
                net_vector += grid_positions[k_i[i]] - points[p]
                net_count += 1
                

    return atoms, net_vector, net_count

@nb.njit(cache=True, fastmath=True, nogil=True)
def compute_accessible_volume(points, energy_map: np.ndarray, cutoff=5):
    accessible_points = set()

    for p in range(points.shape[0]):
        px = points[p,0]
        py = points[p,1]
        pz = points[p,2]

        px4 = ((px + 2) >> 2) << 2
        py4 = ((py + 2) >> 2) << 2
        pz4 = ((pz + 2) >> 2) << 2

        # Compute the influence boundaries of this atom
        min_kx = max(px4-4, 0)
        min_ky = max(py4-4, 0)
        min_kz = max(pz4-4, 0)
        max_kx = min(px4+5, energy_map.shape[0])
        max_ky = min(py4+5, energy_map.shape[1])
        max_kz = min(pz4+5, energy_map.shape[2])

        # Use a value of 4 to match the step of the low-resolution grid,
        # as it's one-fourth the resolution of the original grid.
        for kx in range(min_kx, max_kx, 4): 
            dx2 = (kx-px)*(kx-px)

            for ky in range(min_ky, max_ky, 4): 
                dy2 = (ky-py)*(ky-py)
                dx2y2 = dx2 + dy2

                for kz in range(min_kz, max_kz, 4):
                    dz2 = (kz-pz)*(kz-pz)

                    # Distance squared between the atoms 
                    # and the point on the grid
                    r2 = dx2y2 + dz2 + 1e-2

                    # Avoid the computation if the atom
                    # is too far away
                    if r2 > cutoff**2:
                        continue

                    if energy_map[kx,ky,kz,0] < 0:
                        accessible_points.add((kx,ky,kz))

    out = np.empty((len(accessible_points), 3), dtype=np.int16)

    for i, e in enumerate(accessible_points):
        out[i,0] = e[0]
        out[i,1] = e[1]
        out[i,2] = e[2]

    return out


@nb.njit(cache=True, fastmath=True, nogil=True) 
def get_pockets(positions: np.ndarray, inverse_spacing: float, grid_cutoff: float, energy_threshold: float,
                probe_sigmas2: np.ndarray, probe_epsilons: np.ndarray, sigma_charges: np.ndarray,
                min_spot_points: int, max_spot_distance: float, max_bulk_distance: float, min_points: int):
    
    box = compute_box(positions)

    # Define the origin and the shape of the grid
    origin = box[0]
    shape = ((box[1] - box[0]) * inverse_spacing + 0.5).astype(np.int32) + 1
    shape_low = shape // 4 + 2

    # Convert the positions of the atoms in grid indexes (floating point)
    grid_positions = (positions - origin) * inverse_spacing

    # Allocate arrays on the low-resolution grid
    neighbours = np.zeros((shape_low[0], shape_low[1], shape_low[2], 512), dtype=np.int32)

    # Allocate the map for the bulk energy
    energy_map = np.full((shape[0], shape[1], shape[2], 2), 100, dtype=np.float32)

    # Compute the bulk energy on the low-resolution grid
    compute_low_resolution_grid(energy_map, grid_cutoff, neighbours,
                    grid_positions, probe_sigmas2, 
                    probe_epsilons, sigma_charges)
    
    # Compute the energy on the high-resolution grid and separate the regions
    regions = enhance_volumes(
            energy_map, energy_threshold, 
            neighbours, grid_positions,
            probe_sigmas2, probe_epsilons, sigma_charges, min_spot_points=min_spot_points)

    groups = cluster_regions(regions=regions, grid_positions=grid_positions, neighbours=neighbours,
                    probe_sigmas2=probe_sigmas2, max_bulk_distance=max_bulk_distance,
                    max_spot_distance=max_spot_distance, min_points=min_points)

    return [compute_descriptors(origin, group, energy_map=energy_map) for group in groups]


@nb.njit(cache=True, fastmath=True, nogil=True)
def cluster_regions(regions: list[tuple], grid_positions: np.ndarray, neighbours: np.ndarray,
                    probe_sigmas2: np.ndarray, max_bulk_distance: float, max_spot_distance: float,
                    min_points: int) -> tuple[np.ndarray, np.ndarray]:
    '''
        This function clusters the regions into groups. The regions are considered connected if
        they are close enough with respect to the bulk distance and if they are reachable from each other.
        Additionally spots regions are added to the groups if they are close enough with respect to the spot distance.
    '''
    # Initialize the connectivity matrix between the regions
    connectivity_bluk = np.zeros((len(regions), len(regions)), dtype=np.int8)
    connectivity_spot = np.zeros((len(regions), len(regions)), dtype=np.int8)

    # Compute the connectivity between the regions
    for i in range(len(regions)):
        
        for j in range(i):

            # Check if the two regions are close enough with respect to the spot distance
            if intersect_box(regions[i][3], regions[j][3], max_spot_distance):
                if is_reachable(regions[i][2], regions[j][2], max_spot_distance, grid_positions, probe_sigmas2, neighbours):
                    connectivity_spot[i,j] = 1
                    connectivity_spot[j,i] = 1

                    # Check if the two regions are close enough with respect to the bulk distance
                    if intersect_box(regions[i][3], regions[j][3], max_bulk_distance):
                        if is_reachable(regions[i][2], regions[j][2], max_bulk_distance, grid_positions, probe_sigmas2, neighbours):
                            connectivity_bluk[i,j] = 1
                            connectivity_bluk[j,i] = 1
        connectivity_bluk[i,i] = 1
        connectivity_spot[i,i] = 1
                    

    # Find all possible groups of connected regions
    groups = []
    for region in regions:

        # Groups that are reachable from the current region
        reachable_groups = []

        for group in groups:

            # Check if group is reachable comparing all the regions of the group
            for group_region in group:
                if connectivity_bluk[region[0], group_region[0]]:
                    reachable_groups.append(group)
                    break
            
        # Security check to avoid too many groups
        reachable_groups = reachable_groups[:4]

        # Produce all the partitions of reachable_groups
        for i in range(2 ** len(reachable_groups)):
            group = [region]

            p = 0
            while i > 0:
                if i % 2:
                    group.extend(reachable_groups[p])
                i = i // 2
                p += 1

            groups.append(group)

    # Add spots that are reachable from the groups
    for region in regions:
        for group in groups:
            
            # Stupid numba way to check if region in group
            region_in_group = False
            for group_region in group:
                if region[0] == group_region[0]:
                    region_in_group = True
                    break

            if not region_in_group:
                for group_region in group:
                    if connectivity_spot[region[0], group_region[0]]:
                        groups.append(group + [region])
                        break

    # Stupid numba way to do big_groups = [group.n_points >= min_points for group in groups]                     
    big_groups = []

    for group in groups:
        n_points = 0

        for region in group:
            n_points += region[1].shape[0]

        if n_points >= min_points:
            big_groups.append(group)

    # Remove duplicate groups
    final_groups = []
    unique_groups = []

    for group in big_groups:
        group_set = set([region[0] for region in group])

        if group_set not in unique_groups:
            unique_groups.append(group_set)
            final_groups.append(group)

    return final_groups

@nb.njit(cache=True, fastmath=True, nogil=True)
def compute_descriptors(origin: np.ndarray, group: list[tuple], energy_map: np.ndarray):
    '''
        This function computes the descriptors of a group of regions.

        - volume_points: the points inside the regions
        - surface_points: the points on the surface of the regions
        - regions: the indexes of the regions
        - atoms: the indexes of the atoms around the surface_points
        - enclosure: the enclosure of the surface_points
        - lengths: the lengths of the surface_points
        - energy: the energy of the volume_points
        - electro: the electrostatic potential of the volume_points
    '''
    # Concatenate the points of the regions
    n_volume_points = 0
    n_surface_points = 0
    n_accessible_points = 0

    for region in group:
        n_volume_points += region[1].shape[0]
        n_surface_points += region[2].shape[0]
        n_accessible_points += region[5].shape[0]

    volume_points = np.empty((n_volume_points, 3), dtype=np.int16)
    surface_points = np.empty((n_surface_points, 3), dtype=np.int16)
    accessible_points = np.empty((n_accessible_points, 3), dtype=np.int16)
    regions = set([region[0] for region in group])

    i_volume_points = 0
    i_surface_points = 0
    i_accessible_points = 0
    for region in group:
        volume_points[i_volume_points:i_volume_points+region[1].shape[0]] = region[1]
        surface_points[i_surface_points:i_surface_points+region[2].shape[0]] = region[2]
        accessible_points[i_accessible_points:i_accessible_points+region[5].shape[0]] = region[5]
        i_volume_points += region[1].shape[0]
        i_surface_points += region[2].shape[0]
        i_accessible_points += region[5].shape[0]

    # Concatenate the atoms of the regions
    atoms = set()
    net_vector = np.zeros(shape=3, dtype=np.float32)
    net_count = 0

    for region in group:
        atoms.update(region[4][0])
        net_vector += region[4][1]
        net_count += region[4][2]

    net_norm = np.linalg.norm(net_vector / net_count)

    # The enclosure is computed using the length of the net vector
    # when the net vector is short the enclosure is near 1
    # when the net vector is long the enclosure is near 0
    enclosure = np.exp(-0.3 * net_norm**2)

    # Compute three orthogonal lengths of the pocket diagonalizing the autocorrelation matrix
    residuals = surface_points - average_3d(surface_points)
    variances = np.linalg.eigvalsh(np.dot(residuals.T, residuals)) / residuals.shape[0]
    lengths = np.sqrt(variances)[::-1]

    # Extract the electrostatic potential
    energy = fetch_3d(energy_map, volume_points, 0) 
    electro = fetch_3d(energy_map, volume_points, 1)

    # Return the descriptors
    return origin.copy(), volume_points, surface_points, regions, np.array(sorted(atoms)), enclosure, lengths, energy, electro, accessible_points

@nb.njit(cache=True, fastmath=True)
def add_pocket_to_density(grid: np.ndarray, vp: np.ndarray, grid_origin: np.ndarray, grid_spacing: float, nout: int):
    '''
        This function adds the volume points of a pocket to a density grid.
    '''
    nx, ny, nz = grid.shape
    for i in range(vp.shape[0]):
        # convert vp coordinates to grid points:
        vpidx_x = np.int32((vp[i,0] - grid_origin[0]) // grid_spacing)
        vpidx_y = np.int32((vp[i,1] - grid_origin[1]) // grid_spacing)
        vpidx_z = np.int32((vp[i,2] - grid_origin[2]) // grid_spacing)
        # check if the point is inside the grid
        if 0 <= vpidx_x< nx and 0 <= vpidx_y < ny and 0 <= vpidx_z < nz:
            grid[vpidx_x, vpidx_y, vpidx_z] += 1
        else:
            nout += 1