import numpy as np

class Pocket:
    '''
        A pocket with its properties.
        The points of the pockets are described as points in a grid with
        an origin and a spacing. Electrostatic potential for each point is 
        present.              
    '''
    def __init__(self, origin: np.ndarray, spacing: np.ndarray, 
                 volume_points: np.ndarray, surface_points: np.ndarray, accessible_points: float, 
                 regions: set[int], energy: np.ndarray, electro: np.ndarray,
                 atoms: np.ndarray, charges: np.ndarray, lengths: np.ndarray, enclosure: float) -> None:
        '''
            Initialize the pocket with the specified values (points and descriptors)   
        '''
        self.origin = origin
        self.spacing = spacing
        self.energy = energy
        self.electro = electro
        self.volume_points = volume_points
        self.surface_points = surface_points
        self.accessible_points = accessible_points
        self.regions = regions
        self.atoms = atoms
        self.charges = charges
        self.lengths = lengths
        self.enclosure = enclosure
        self.n_apolar_atoms, self.n_polar_atoms, self.n_charges_atoms = np.histogram(np.abs(charges), bins=[0, 0.3, 0.7, 10])[0]
        self.hydrophobicity = self.n_apolar_atoms / self.atoms.shape[0]
        self.score = None

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
    def accessible(self):
        return self.accessible_points.shape[0] * self.volume_unit**(2/3) * 64
    
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
    def accessible_positions(self):
        '''
        Points of the pocket as floating point coordinate 
        '''
        return self.origin + self.accessible_points @ self.spacing.T

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
        return [self.volume, self.enclosure, *self.energy_quants, *self.charges_histogram]

    def props(self):
        return {
            k: getattr(self, v[0]) if len(v) == 2 else getattr(self, v[0])[v[2]]
            for k, v in self.properties.items()
        }

    def __repr__(self) -> str:
        return f'<Pocket volume={self.volume:.1f} #atoms={len(self.atoms)}>'
    
    properties = {
        'Score': ('score','#f74d75'),
        'Volume': ('volume','#90be6d'),
        'Surface': ('surface','#577590'), 
        'Accessible': ('accessible','#f8961e'),
        'Length1': ('lengths', '#43aa8b', 0),
        'Length2': ('lengths', '#43aa8b', 1),
        'Length3': ('lengths', '#43aa8b', 2),
        'Hydrophobicity': ('hydrophobicity', '#f9c74f'),
        'Enclosure': ('enclosure', '#2a6f97'),
        '#PolarAtoms': ('n_polar_atoms', '#f8961e'),
        '#ChargedAtoms': ('n_charges_atoms', '#9b5de5'),
    }