import numpy as np

def box_coordinates(box_origin: np.ndarray, box_vectors: np.ndarray, points: np.ndarray):
    inverse = np.linalg.inv(box_vectors.T)
    return (points - box_origin) @ inverse.T 

def is_inside(box_origin: np.ndarray, box_vectors: np.ndarray, points: np.ndarray):
    box_coord = box_coordinates(box_origin, box_vectors, points)
    return np.all((box_coord >= 0) & (box_coord <= 1), axis=1)

    
