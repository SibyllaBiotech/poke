import numpy as np
from numba import njit

@njit(cache=True, fastmath=True)
def sp2_1h(
    central: np.ndarray, 
    neighbor1: np.ndarray, 
    neighbor2: np.ndarray, 
    h_dist: float = 1.0
) -> np.ndarray:
    """
    Calculate the position of a hydrogen atom attached to a sp2 atom (eg. C₂CH, in rings). GMX case 1
    """
    v_cn1 = central - neighbor1
    v_cn2 = central - neighbor2
    v_sum = v_cn1/np.linalg.norm(v_cn1) + v_cn2/np.linalg.norm(v_cn2)
    
    return central + h_dist * v_sum / np.linalg.norm(v_sum)


@njit(cache=True, fastmath=True)
def sp3_lone_1h(
    central: np.ndarray, 
    neighbor: np.ndarray, 
    neighbor_neighbor: np.ndarray, 
    h_dist: float = 1.0
) -> np.ndarray:
    """
    Calculate the position of a hydrogen atom attached to a sp3 atom with two lone pairs (eg. -OH). GMX case 2
    """
    # α = 109.47°, angle between the new bond and the already existing bond
    cos_α = -1/3
    sin_α = 2 * np.sqrt(2) / 3
    
    v_nc = neighbor - central
    v_nc /= np.linalg.norm(v_nc)

    v_nnn = neighbor - neighbor_neighbor

    v_norm = np.cross(v_nnn, v_nc)
    v_norm /= np.linalg.norm(v_norm)

    v_k = np.cross(v_nc, v_norm)

    return central + h_dist * (cos_α * v_nc + sin_α * v_k)


@njit(cache=True, fastmath=True)
def sp2_2h(
    central: np.ndarray, 
    neighbor: np.ndarray, 
    neighbor_neighbor: np.ndarray, 
    h_dist: float = 1.0
) -> np.ndarray:
    """
    Calculate the position of two hydrogen atoms attached to a sp2 atom (eg. –NH₂). GMX case 3
    """
    # α = 120°, angle between the new bond and the already existing bond
    cos_α = -1/2
    sin_α = np.sqrt(3) / 2
    
    v_nc = neighbor - central
    v_nc /= np.linalg.norm(v_nc)

    v_nnn = neighbor - neighbor_neighbor

    v_norm = np.cross(v_nnn, v_nc)
    v_norm /= np.linalg.norm(v_norm)

    v_k = np.cross(v_nc, v_norm)

    h_pos = np.empty((2, 3), dtype=np.float32)
    h_pos[0] = central + h_dist * (cos_α * v_nc + sin_α * v_k)
    h_pos[1] = central + h_dist * (cos_α * v_nc - sin_α * v_k)

    return h_pos


@njit(cache=True, fastmath=True)
def sp3_3h(
    central: np.ndarray, 
    neighbor: np.ndarray, 
    neighbor_neighbor: np.ndarray, 
    h_dist: float = 1.0
) -> np.ndarray:
    """
    Calculate the position of three hydrogen atoms attached to a sp3 atom (eg. -CH₃). GMX case 4
    """
    # α = 109.47°, angle between the new bond and the already existing bond
    # β = 120°
    cos_α = -1/3
    sin_α = 2 * np.sqrt(2) / 3
    sin_β = np.sqrt(3) / 2
    
    v_nc = neighbor - central
    v_nc /= np.linalg.norm(v_nc)

    v_nnn = neighbor - neighbor_neighbor

    v_norm = np.cross(v_nnn, v_nc)
    v_norm /= np.linalg.norm(v_norm)

    v_k = np.cross(v_nc, v_norm)

    h_pos = np.empty((3, 3), dtype=np.float32)
    h_pos[0] = central + h_dist * (cos_α * v_nc + sin_α * v_k)
    h_pos[1] = central + h_dist * (cos_α * v_nc + sin_α * v_norm * sin_β - sin_α * v_k / 2)
    h_pos[2] = central + h_dist * (cos_α * v_nc - sin_α * v_norm * sin_β - sin_α * v_k / 2)

    return h_pos
    

@njit(cache=True, fastmath=True)
def sp3_1h(
    central: np.ndarray, 
    neighbor1: np.ndarray, 
    neighbor2: np.ndarray, 
    neighbor3: np.ndarray, 
    h_dist: float = 1.0
) -> np.ndarray:
    """
    Calculate the position of a hydrogen atom attached to a sp3 atom (eg. C₃CH). GMX case 5
    """
    centroid = (neighbor1 + neighbor2 + neighbor3) / 3
    v_cc = central - centroid

    return central + h_dist * v_cc / np.linalg.norm(v_cc)


@njit(cache=True, fastmath=True)
def sp3_2h(
    central: np.ndarray, 
    neighbor1: np.ndarray, 
    neighbor2: np.ndarray, 
    h_dist: float = 1.0
) -> np.ndarray:
    """
    Calculate the position of two hydrogen atoms attached to a sp3 atom (eg. -CH₂-). GMX case 6
    """
    # α = 109.47°, angle between the two new bonds. Here α / 2 = 54.735° is needeed. 
    cos_half_α = 1 / np.sqrt(3)
    sin_half_α = np.sqrt(2 / 3)
    
    v_cm = central - (neighbor1 + neighbor2) / 2
    v_cm /= np.linalg.norm(v_cm)
    
    v_norm  = np.cross(central - neighbor1, central - neighbor2)
    v_norm /= np.linalg.norm(v_norm)

    h_pos = np.empty((2, 3), dtype=np.float32)
    h_pos[0] = central + h_dist * (cos_half_α * v_cm + sin_half_α * v_norm)
    h_pos[1] = central + h_dist * (cos_half_α * v_cm - sin_half_α * v_norm)

    return h_pos


@njit(cache=True, fastmath=True)
def sp3_2o(
    central: np.ndarray, 
    neighbor: np.ndarray, 
    neighbor_neighbor: np.ndarray, 
    o_dist: float = 1.36
) -> np.ndarray:
    """
    Calculate the position of two oxygen atoms attached to a sp3 atom (eg. -COO⁻). GMX case 8
    """
    # α = 117°, angle between the new bond and the already existing bond
    cos_α = np.cos(np.deg2rad(117))
    sin_α = np.sin(np.deg2rad(117))
    
    v_nc = neighbor - central
    v_nc /= np.linalg.norm(v_nc)

    v_nnn = neighbor - neighbor_neighbor

    v_norm = np.cross(v_nnn, v_nc)
    v_norm /= np.linalg.norm(v_norm)

    v_k = np.cross(v_nc, v_norm)

    o_pos = np.empty((2, 3), dtype=np.float32)
    o_pos[0] = central + o_dist * (cos_α * v_nc + sin_α * v_k)
    o_pos[1] = central + o_dist * (cos_α * v_nc - sin_α * v_k)

    return o_pos


@njit(cache=True, fastmath=True)
def sp3_ooh(
    central: np.ndarray, 
    neighbor: np.ndarray, 
    neighbor_neighbor: np.ndarray, 
    o_dist: float = 1.23, 
    oh_dist: float = 1.25, 
    h_dist: float = 1.0
) -> np.ndarray:
    """
    Calculate the position of two oxygen and a hydrogen atoms attached to a sp3 atom (eg. -COOH). GMX case 9
    """
    # α = 121° 
    cos_α = np.cos(np.deg2rad(121)) # for O⁻ bond
    sin_α = np.sin(np.deg2rad(121)) # for O⁻ bond
    cos_β = np.cos(np.deg2rad(115)) # for OH bond
    sin_β = np.sin(np.deg2rad(115)) # for OH bond

    
    v_nc = neighbor - central
    v_nc /= np.linalg.norm(v_nc)

    v_nnn = neighbor - neighbor_neighbor

    v_norm = np.cross(v_nnn, v_nc)
    v_norm /= np.linalg.norm(v_norm)

    v_k = np.cross(v_nc, v_norm)

    oh_pos = np.empty((3, 3), dtype=np.float32)
    oh_pos[0] = central + o_dist * (cos_α * v_nc + sin_α * v_k)
    oh_pos[1] = central + oh_dist * (cos_β * v_nc - sin_β * v_k)
    oh_pos[2] = sp3_lone_1h(oh_pos[1], central, neighbor, h_dist)

    return oh_pos