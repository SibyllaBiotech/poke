import numpy as np
import numba as nb

@nb.njit(cache=True, fastmath=True, nogil=True)
def compute_box(array: np.ndarray):
    '''
    Compute the minimal box that contains all the points in the array.

    Parameters
    ------------
        array: np.ndarray (N,3)
            The array of points

    Return values
    ------------
        box: np.ndarray
            The box containing all the points
    '''
    box = np.empty((2,3), dtype=array.dtype)

    box[0] = array[0]
    box[1] = array[1]

    for i in range(1, array.shape[0]):
        box[0,0] = min(box[0,0], array[i,0])
        box[1,0] = max(box[1,0], array[i,0])

        box[0,1] = min(box[0,1], array[i,1])
        box[1,1] = max(box[1,1], array[i,1])

        box[0,2] = min(box[0,2], array[i,2])
        box[1,2] = max(box[1,2], array[i,2])

    return box

@nb.njit(cache=True, fastmath=True, nogil=True)
def average_3d(array: np.ndarray):
    '''
    Compute the average of the points in the array.

    Parameters
    ------------
        array: np.ndarray (N,3)
            The array of points
        
    Return values
    ------------
        average: np.ndarray (3,)
            The average of the points
    '''
    sum = np.zeros(3, dtype=array.dtype)

    for i in range(array.shape[0]):
        sum[0] += array[i,0]
        sum[1] += array[i,1]
        sum[2] += array[i,2]

    return sum / array.shape[0]

@nb.njit(cache=True, fastmath=True, nogil=True)
def fetch_3d(map: np.ndarray, points: np.ndarray, i: int):
    '''
    Fetch the values of the map at the specified points.

    Parameters
    ------------
        map: np.ndarray (N,M,P)
            The map
        points: np.ndarray (K,3)
            The points
        i: int
            The index of the map to fetch

    Return values
    ------------
        values: np.ndarray (K,)
            The values of the map at the specified points
    '''
    values = np.empty(points.shape[0], dtype=map.dtype)

    for i in range(points.shape[0]):
        values[i] = map[points[i,0], points[i,1], points[i,2], i]
    return values


@nb.njit(cache=True, fastmath=True, nogil=True)
def intersect_box(box1: np.ndarray, box2: np.ndarray, threshold: float):
    ''' 
        Return true if there is an intersection between box1 and box2
    '''
    return (box1[0,0] - threshold < box2[1,0] and box1[1,0] + threshold > box2[0,0] and
            box1[0,1] - threshold < box2[1,1] and box1[1,1] + threshold > box2[0,1] and
            box1[0,2] - threshold < box2[1,2] and box1[1,2] + threshold > box2[0,2])

@nb.njit(cache=True, fastmath=True, nogil=True)   
def is_reachable(points1, points2, threshold, grid_positions, probe_sigmas2, neighbours):
    '''
        Return True if there is a segment shorter than the threshold connecting 
        points1 with points2, and this segment does not intersect any atom. 
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

                # Rescale the direction vector making the largest coordinate equal to 1 
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

                        # Retrieve the atom's position in grid coordinates
                        nx = grid_positions[i,0]
                        ny = grid_positions[i,1]
                        nz = grid_positions[i,2]

                        # Van der Waals composition rules between the atom and the probe
                        sigma2 = probe_sigmas2[i]

                        px_nx = px-nx
                        py_ny = py-ny
                        pz_nz = pz-nz

                        r2 = px_nx*px_nx + py_ny*py_ny + pz_nz*pz_nz

                        # Does the line intersect the atom?
                        if r2 < 0.8 * sigma2:
                            good = False
                            break

                    # There is no need to proceed with this segment 
                    # if the line intersects the atom.
                    # we can skip and move on to the next segment.
                    if not good:
                        break
                
                # If the entire line has been traversed without intersecting any atoms, 
                # then the regions are connected.
                if good:
                    return True
    
    # Connecting the regions was unfeasible
    return False