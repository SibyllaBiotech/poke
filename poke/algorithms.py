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
    # energy_map[::4,::4,::4, 0] -= np.abs(energy_map[::4,::4,::4, 1])


@nb.njit(cache=True, fastmath=True)
def enhance_volumes(energy_map: np.ndarray, isovalue: float, 
                    neighbours: np.ndarray, grid_positions: np.ndarray,
                    probe_sigmas2: np.ndarray, probe_epsilons: np.ndarray, sigma_charges: np.ndarray):
    '''
        Expand the low-energy regions found in the low-resolution map. The expansion is done 
        on the high-resolution map computing the potential in the missing positions.
        This function return the set of points under a threshold of total energy.
    '''
    # Points that have already been inserted in pockets_points
    inserted = np.zeros((energy_map.shape[0], energy_map.shape[1], energy_map.shape[2]), dtype=np.uint16)

    # Create lists of points, one for each pocket 
    pockets_points = np.empty((energy_map.size, 3), dtype=np.int16) # 65536 maximum dimension of the grid along each axis
    surfaces_points = np.empty((energy_map.size // 2, 3), dtype=np.int16) 

    # For each pocket, mark the beginning of the volume/surface points array (excluding the first one, which starts in zero)
    # these serve as slicers of the pocket_points and surface_points arrays, avoids computation of the sum of the volumes array
    # example [17,69,666] can be used as pockets_points[17:69] to extract the points belonging to the second pocket
    pockets_slicer  = np.zeros(65536, dtype=np.uint32)
    surfaces_slicer = np.zeros(65536, dtype=np.uint32)

    # Pocket counter, needed for ids
    # This also serves to keep track of the boundaries of the volumes and surfaces arrays
    n_pockets = 0

    n_visited = 0
    n_inserted = 0  
    n_surfaces = 0

    for i in range(0, energy_map.shape[0], 4):
        for j in range(0, energy_map.shape[1], 4):
            for k in range(0, energy_map.shape[2], 4):

                # Find first point belonging to new pocket (root of graph)
                if energy_map[i,j,k,0] < isovalue and not inserted[i,j,k]:

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

                        # Compute the energies of the near points 
                        # that will be used in the following ifs
                        energy_xp, energy_xm, energy_yp, energy_ym, energy_zp, energy_zm = \
                            compute_star_points(px,py,pz, energy_map, 
                                            neighbours, grid_positions, 
                                            probe_sigmas2, probe_epsilons, sigma_charges)
                        
                        # Decide which points are under the threshold and
                        # must be added in the returning lists
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
                    
                    # check for IndexError & return 
                    if n_pockets+1 >= pockets_slicer.shape[0]:
                        return (inserted, pockets_points[:n_inserted],  pockets_slicer, 
                                          surfaces_points[:n_surfaces], surfaces_slicer)
    
    return (inserted, pockets_points[:n_inserted],   pockets_slicer[:n_pockets+1], 
                      surfaces_points[:n_surfaces], surfaces_slicer[:n_pockets+1])


@nb.njit(cache=True, fastmath=True)
def compute_star_points(px: int, py: int, pz: int, energy_map: np.ndarray,
                        neighbours: np.ndarray, grid_positions: np.ndarray, 
                        probe_sigmas2: np.ndarray, probe_epsilons: np.ndarray, sigma_charges: np.ndarray):
    '''
        Compute the potential of the 6 points around the central point (px, py, pz).
        The computation take in consideration only the atoms in the neighbour list.
    '''

    # Fetch energy values from the map
    energy_xp = energy_map[px+1,py,pz,0] if px+1 < energy_map.shape[0] else 1000
    energy_xm = energy_map[px-1,py,pz,0] if px > 0 else 1000

    energy_yp = energy_map[px,py+1,pz,0] if py+1 < energy_map.shape[1] else 1000
    energy_ym = energy_map[px,py-1,pz,0] if py > 0 else 1000

    energy_zp = energy_map[px,py,pz+1,0] if pz+1 < energy_map.shape[2] else 1000
    energy_zm = energy_map[px,py,pz-1,0] if pz > 0 else 1000

    # Checks if those values were never been computed    
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

        # Van de Walls composition rules between the atom and the probe
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
        # on the 6 points around
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
        # energy_xp -= np.abs(elec_xp)
        energy_map[px+1,py,pz,0] = energy_xp
        energy_map[px+1,py,pz,1] = elec_xp
    if compute_xm:
        # energy_xm -= np.abs(elec_xm)
        energy_map[px-1,py,pz,0] = energy_xm
        energy_map[px-1,py,pz,1] = elec_xm
    if compute_yp:
        # energy_yp -= np.abs(elec_yp)
        energy_map[px,py+1,pz,0] = energy_yp
        energy_map[px,py+1,pz,1] = elec_yp
    if compute_ym:
        # energy_ym -= np.abs(elec_ym)
        energy_map[px,py-1,pz,0] = energy_ym
        energy_map[px,py-1,pz,1] = elec_ym
    if compute_zp:
        # energy_zp -= np.abs(elec_zp)
        energy_map[px,py,pz+1,0] = energy_zp
        energy_map[px,py,pz+1,1] = elec_zp
    if compute_zm:
        # energy_zm -= np.abs(elec_zm)
        energy_map[px,py,pz-1,0] = energy_zm
        energy_map[px,py,pz-1,1] = elec_zm

    return energy_xp, energy_xm, energy_yp, energy_ym, energy_zp, energy_zm

@nb.njit(cache=True, fastmath=True)
def compute_real_elec(points: np.ndarray, positions: np.ndarray, charges: np.ndarray):
    elec = np.empty(points.shape[0])
    for i in range(points.shape[0]):
        el = 0
        for j in range(positions.shape[0]):
            el += charges[j] / np.linalg.norm(points[i] - positions[j])**2
        elec[i] = el
    return elec

    

@nb.njit(cache=True, fastmath=True)
def intersect_box(box1: np.ndarray, box2: np.ndarray, threshold: float):
    ''' 
        Return true if box1 and box2 have an intersection.
    '''
    return (box1[0,0] - threshold < box2[1,0] and box1[1,0] + threshold > box2[0,0] and
            box1[0,1] - threshold < box2[1,1] and box1[1,1] + threshold > box2[0,1] and
            box1[0,2] - threshold < box2[1,2] and box1[1,2] + threshold > box2[0,2])

@nb.njit(cache=True, fastmath=True)   
def is_reachable(points1, points2, threshold, grid_positions, probe_sigmas2, neighbours):
    '''
        Return True if exists a segment shorter than threshold that connects points1 with points2
        without crossing a barrier in the energy_map higher than max_barrier. 
    '''
    threshold2 = threshold**2

    for p1 in range(points1.shape[0]):
        p1x = points1[p1,0]
        p1y = points1[p1,1]
        p1z = points1[p1,2]

        for p2 in range(points2.shape[0]):
            p2x = points2[p2,0]
            p2y = points2[p2,1]
            p2z = points2[p2,2]

            dx = p2x-p1x
            dy = p2y-p1y
            dz = p2z-p1z
            d2 = dx*dx + dy*dy + dz*dz

            if d2 < threshold2:
                # Check that the segment connecting p1 and p2 does not cross 
                # a barrier in the energy map. 
                max_delta = max(np.abs(dx),np.abs(dy),np.abs(dz))

                # Rescale the direction vector making equal to 1 the largest coordinate
                # It's an equivalent version of the Bresenham's line algorithm
                slope_x = dx / max_delta
                slope_y = dy / max_delta
                slope_z = dz / max_delta
                
                good = True
                for k in range(2, max_delta-1):

                    # Compute the point on the segment to test
                    px = p1x + slope_x * k
                    py = p1y + slope_y * k
                    pz = p1z + slope_z * k

                    kx4 = np.int32(px * 0.25 + 0.5)
                    ky4 = np.int32(py * 0.25 + 0.5)
                    kz4 = np.int32(pz * 0.25 + 0.5)

                    n_neighbours = neighbours[kx4,ky4,kz4,0]
                    i_neighbours = neighbours[kx4,ky4,kz4,1:]

                    for i_neighbour in range(n_neighbours):
                        i = i_neighbours[i_neighbour]

                        # Fetch the atom position on the grid coordinates
                        nx = grid_positions[i,0]
                        ny = grid_positions[i,1]
                        nz = grid_positions[i,2]

                        # Van de Walls composition rules between the atom and the probe
                        sigma2 = probe_sigmas2[i]

                        px_nx = px-nx
                        py_ny = py-ny
                        pz_nz = pz-nz

                        r2 = px_nx*px_nx + py_ny*py_ny + pz_nz*pz_nz

                        # Is this atom intesected?
                        if r2 < 0.8 * sigma2:
                            good = False
                            break

                    # We don't need to continue on this segment
                    # we can short-circuit and continue to the next
                    if not good:
                        break
                
                # An entire line has been completed without 
                # intersecting any atoms? So the regions are connected!
                if good:
                    return True
    
    # Connecting the regions was not be possible
    return False


@nb.njit(cache=True, fastmath=True)   
def get_region_atoms(points, neighbours, grid_positions, k_nearest):
    '''
        For each point of the region, this function find the k-nearest atoms.
        After collecting all the atoms for each point, it returns the unique sorted list.
    '''
    atoms = set()

    # Accumulation vector for computing the enclosure
    net_vector = np.zeros(shape=3, dtype=np.float32)
    net_count = 0

    # These arrays are used multiple times to store the k-nearest atoms
    # for each atom. The atoms are collected using a insert sort 
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

        # Dump the k-nearest arrays into the atoms set
        for i in range(k_i.shape[0]):
            if k_i[i] != -1:
                atoms.add(k_i[i])

                # Add to the accumulation vector the vector from the pocket point to the atom              
                net_vector += grid_positions[k_i[i]] - points[p]
                net_count += 1
                

    return np.array(sorted(atoms)), np.linalg.norm(net_vector / net_count)