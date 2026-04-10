import numpy as np
import numba as nb


@nb.jit(nopython=True, cache=True, fastmath=True)
def is_covbonded(
    first: np.ndarray, 
    candidates: np.ndarray, 
    cov_dist: float = 1.5
) -> int:
    """
    Check if a covalent bond is formed between the first atom and any of the candidates, 
    and return the index of the candidate.
    """

    dist2_th = cov_dist ** 2
    for ci, candidate in enumerate(candidates):
        if 0.81 * dist2_th < np.sum((candidate - first) ** 2) < 1.21 * dist2_th:
            return ci
    return -1


@nb.jit(nopython=True, cache=True, fastmath=True)
def is_hbonded(
    donor: np.ndarray, 
    h_pos: np.ndarray, 
    acceptors: np.ndarray, 
    hb_dist: float = 3.0, 
    angle: float = 135.0
) -> bool:
    """
    Check if a hydrogen bond is formed between a donor and any acceptor.
    """

    dist2_th = hb_dist ** 2
    cos_angle_th = np.cos(np.deg2rad(angle))

    for acceptor in acceptors:
        if np.sum((acceptor - donor) ** 2) <= dist2_th:
            v1 = donor - h_pos
            v2 = acceptor - h_pos

            cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            if cos_angle <= cos_angle_th:
                return True

    return False
    