import numpy as np

class Pocket:
    '''
        A pocket with its properties.
        The points of the pockets are described as points in a grid with
        an origin and a spacing. Electrostatic potential for each point is 
        present.              
    '''
    def __init__(self, origin: np.ndarray, spacing: np.ndarray, 
                 volume_points: np.ndarray, surface_points: np.ndarray, electro: np.ndarray,
                 monopole: np.ndarray, dipole: np.ndarray,  
                 n_regions: int, atoms: np.ndarray, lengths: np.ndarray,
                 hydrophobicity: float, enclousure: float,
                 n_polar_atoms: int, n_charges_atoms: int) -> None:
        '''
            Initialize the pocket with the specified values (points and descriptors)   
        '''
        self.origin = origin
        self.spacing = spacing
        self.electro = electro
        self.monopole = monopole
        self.dipole = dipole
        self.volume_points = volume_points
        self.surface_points = surface_points
        self.n_regions = n_regions
        self.atoms = atoms
        self.lengths = lengths
        self.hydrophobicity = hydrophobicity
        self.enclousure = enclousure
        self.n_polar_atoms = n_polar_atoms
        self.n_charges_atoms = n_charges_atoms

    @property
    def volume_unit(self) -> float:
        return np.linalg.det(self.spacing)

    @property
    def volume(self) -> float:
        return self.volume_points.shape[0] * self.volume_unit
    
    @property
    def surface(self):
        return self.surface_points.shape[0] * self.volume_unit**(2/3)
    
    @property
    def volume_positions(self):
        '''
        Points of the pocket as floating point coordinate 
        '''
        return self.origin + np.dot(self.spacing, self.volume_points.T).T

    @property
    def electro_histogram(self):
        '''
        Compute a large-bins histogram for potential energy of the volume points
        '''
        return np.histogram(self.electro, bins=[-3.5, -1.5, -0.5, +0.5, +1.5, +3.5])
    
    def translate(self, vector: np.ndarray):
        '''
        Translate the pocket along the specified vector
        '''
        self.origin += vector

    def apply_matrix(self, transformation_matrix: np.ndarray):
        '''
        Transform the points of the pocket applying a linear transformation.
        '''
        self.spacing = np.dot(transformation_matrix, self.spacing) 

    def minimal_bounding_box_rotation(self):
        p = self.volume_positions

        θ = 0.1 
        cosθ = np.cos(θ)
        sinθ = np.sin(θ)

        Rx = np.array([
            [    1,     0,     0],
            [    0,  cosθ,  sinθ],
            [    0, -sinθ,  cosθ]
        ])
        iRx = np.linalg.inv(Rx)

        Ry = np.array([
            [ cosθ,     0, -sinθ],
            [    0,     1,     0],
            [ sinθ,     0,  cosθ]
        ])
        iRy = np.linalg.inv(Ry)

        Rz = np.array([
            [ cosθ, -sinθ,     0],
            [ sinθ,  cosθ,     0],
            [    0,    0,      1]
        ])
        iRz = np.linalg.inv(Rz)
        rot_steps = [Rx, iRx, Ry, iRy, Rz, iRz]

        best_rotation = np.eye(3)
        best_volume   = np.prod(np.max(p, axis=0) - np.min(p, axis=0))

        for _ in range(20):
            for rot_step in rot_steps:
                trial_rot = rot_step @ best_rotation
                tp = np.einsum("ij,kj->ki", trial_rot, p)
                trial_vol = np.prod(np.max(tp, axis=0) - np.min(tp, axis=0))
                
                if trial_vol < best_volume:
                    best_rotation = trial_rot
                    best_volume = trial_vol

        return best_rotation

        
        
    def features(self):
        return np.array([1,2,3])

    def __repr__(self) -> str:
        return f'<Pocket volume={self.volume:.1f} #atoms={len(self.atoms)}>'