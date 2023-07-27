import pickle
import numpy as np

from pathlib import Path

class Pocket:
    '''
        A pocket with its properties.
        The points of the pockets are described as points in a grid with
        an origin and a spacing. Electrostatic potential for each point is 
        present.              
    '''

    # Machine learning model to score 
    # the pocket starting from descriptors
    scoring_model = None


    def __init__(self, origin: np.ndarray, spacing: np.ndarray, 
                 volume_points: np.ndarray, surface_points: np.ndarray, 
                 energy: np.ndarray, electro: np.ndarray,
                 atoms: np.ndarray, charges: np.ndarray, lengths: np.ndarray, enclousure: float) -> None:
        '''
            Initialize the pocket with the specified values (points and descriptors)   
        '''
        self.origin = origin
        self.spacing = spacing
        self.energy = energy
        self.electro = electro
        self.volume_points = volume_points
        self.surface_points = surface_points
        self.atoms = atoms
        self.charges = charges
        self.lengths = lengths
        self.enclousure = enclousure
        self.n_apolar_atoms, self.n_polar_atoms, self.n_charges_atoms = np.histogram(np.abs(charges), bins=[0, 0.3, 0.7, 10])[0]
        self.hydrophobicity = self.n_apolar_atoms / self.atoms.shape[0]
        self.score = self.compute_score()

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
        return self.origin + self.volume_points @ self.spacing.T
    
    @property
    def surface_positions(self):
        '''
        Points of the pocket as floating point coordinate 
        '''
        return self.origin + self.surface_points @ self.spacing.T

    @property
    def electro_histogram(self):
        '''
        Compute a large-bins histogram for potential energy of the volume points
        '''
        return np.histogram(self.electro, bins=[-1.5, -1.0, -0.5, +0.5, +1.0, +1.5])[0]
    
    @property
    def charges_histogram(self):
        '''
        Compute a large-bins histogram for atoms charges
        '''
        return np.histogram(self.charges, bins=[-10, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 10])[0]
    
    @property
    def energy_histogram(self):
        '''
        Compute a large-bins histogram for potential energy of the volume points
        '''
        return np.histogram(self.energy, bins=[-100, -25, -24, -23, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13, 0])[0]
    
    @property
    def energy_quants(self):
        '''
        Compute a large-bins histogram for potential energy of the volume points
        '''
        return np.quantile(self.energy, q=[0.05, 0.1, 0.3, 0.5, 0.6, 0.9])

    def translate(self, vector: np.ndarray):
        '''
        Translate the pocket along the specified vector
        '''
        self.origin += vector

    def apply_matrix(self, transformation_matrix: np.ndarray):
        '''
        Transform the points of the pocket applying a linear transformation.
        '''
        self.origin = np.dot(transformation_matrix, self.origin)
        self.spacing = np.dot(transformation_matrix, self.spacing) 

    def minimal_bounding_box(self):
        p = self.surface_positions

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
        best_min      = np.min(p, axis=0)
        best_max      = np.max(p, axis=0)
        best_volume   = np.prod(best_max - best_min)

        for _ in range(20):
            for rot_step in rot_steps:
                trial_rot = rot_step @ best_rotation
                tp = np.einsum('ij,kj->ki', trial_rot, p)
                trial_min = np.min(tp, axis=0)
                trial_max = np.max(tp, axis=0)
                trial_vol = np.prod(trial_max - trial_min)
                
                if trial_vol < best_volume:
                    best_rotation = trial_rot
                    best_volume = trial_vol 
                    best_min = trial_min
                    best_max = trial_max

        origin = best_min @ best_rotation
        box_vectors = np.diag(best_max - best_min) @ best_rotation
        return origin, box_vectors

    def features(self):
        return [self.volume, self.enclousure, *self.energy_quants, *self.charges_histogram]

    def props(self):
        return {
            'Score': self.score, 
            'Volume': self.volume, 
            'Surface': self.surface, 
            'Lenght1': self.lengths[0],
            'Length2': self.lengths[1], 
            'Length3': self.lengths[2], 
            'Hydrophobicity': self.hydrophobicity,
            'Enclousure': self.enclousure,
            '#PolarAtoms': self.n_polar_atoms, 
            '#ChargedAtoms': self.n_charges_atoms,
        }
    
    def compute_score(self):
        if self.scoring_model is None:
            model_path = Path(__file__).parent / 'model.pickle'
            self.scoring_model = pickle.loads(model_path.read_bytes())
        return self.scoring_model(np.array(self.features()))[0]

    def __repr__(self) -> str:
        return f'<Pocket volume={self.volume:.1f} #atoms={len(self.atoms)}>'