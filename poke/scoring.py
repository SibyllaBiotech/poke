import copy
import io
import json
import gzip
import warnings
import itertools
import tempfile
import numpy as np
import MDAnalysis
import pickle

import torch
import torch.nn as nn
import torch.optim as optim
from torcheval.metrics import R2Score

from pathlib import Path
from MDAnalysis import Universe, AtomGroup, Writer
from .pocket import Pocket
from .writer import PocketWriter

from .docking.gnina import Gnina
from .docking.molecule import SDFWriter, SDFReader, Molecule

def chunk_generator(iterable, size=10):
    iterator = iter(iterable)
    while chunk := tuple(itertools.islice(iterator, size)):
        yield chunk

class Scorer(nn.Module):

    def __init__(self, in_features: int):
        super(Scorer, self).__init__()
        n_hidden_1 = 64

        self.layers = nn.Sequential(
            nn.Linear(in_features, n_hidden_1, bias=False),
            nn.BatchNorm1d(n_hidden_1),
            nn.Sigmoid(),
            nn.Dropout(0.1),
            nn.Linear(n_hidden_1, 1),
        )

    def forward(self, x):
        return self.layers(x)
    
    @classmethod
    def train_with(cls, train_pockets: list[Pocket], train_scores: list[float],
                        test_pockets: list[Pocket], test_scores: list[float],   
                        epochs: int=50, learning_rate: float=0.0005, report_every: int=None):

        if torch.cuda.is_available(): 
            device = torch.device("cuda:0") 
        else: 
            device = torch.device("cpu") 
            
        # Compute the features of each pocket and targets
        train_features = torch.tensor(np.array([p.features() for p in train_pockets]), dtype=torch.float32).to(device)
        train_targets = torch.tensor(train_scores, dtype=torch.float32).unsqueeze(1).to(device)

        test_features = torch.tensor(np.array([p.features() for p in test_pockets]), dtype=torch.float32).to(device)
        test_targets = torch.tensor(test_scores, dtype=torch.float32).unsqueeze(1).to(device)

        # Initialize the model
        model = cls(in_features=train_features.shape[1]).to(device)

        optimizer = optim.Adamax(model.parameters(), lr=learning_rate)
        criterion = nn.MSELoss()

        # set variables to keep track of the model that peforms best on the test set, to avoid overfitting
        best_test_r2 = 0
        best_model = None
        
        consecutive_train_better = 0

        try:
            for i in range(epochs):

                # Reset the gradients for the model parameters
                optimizer.zero_grad()

                # Make predictions
                outputs = model(train_features)

                # Compute the loss and its gradients
                loss = criterion(outputs, train_targets)
                loss.backward()

                # Adjust learning weights
                optimizer.step()

                if report_every and i % report_every == 0:
                    model.eval()
                    test_predictions = model(test_features)
                    test_loss = criterion(test_predictions, test_targets)

                    metric = R2Score()
                    metric.update(test_predictions, test_targets)
                    test_r2 = metric.compute()

                    metric = R2Score()
                    metric.update(outputs, train_targets)
                    train_r2 = metric.compute()

                    if loss < test_loss:
                        consecutive_train_better += 1
                    else:
                        consecutive_train_better = 0

                    if test_r2 > best_test_r2:
                        best_test_r2 = test_r2
                        best_model = copy.deepcopy(model)

                    model.train()
                    print(f'Epoch {i:<4} Train Loss: {loss:.4g} Test Loss: {test_loss:.4g} Train r2: {train_r2:.4g} Test r2: {test_r2:.4g}', end='\r')
        except KeyboardInterrupt:
            print('')
            pass

        return best_model.cpu()

    def predict(self, pockets: list[Pocket]):
        self.eval()
        features = torch.tensor(np.array([p.features() for p in pockets]), dtype=torch.float32)
        return self(features).detach().numpy().ravel()
        

def dock_molecules_on_pocket(output_dir: Path, atoms: AtomGroup, pocket: Pocket, smiles: list[str], 
                             n_cpus: int, docking: bool=True, 
                             only_cached: bool=False, compute_lle: bool=False):
    '''
        Dock a set of molecules inside a pocket.
        A folder will be created inside output_dir containing the results of the docking.
        If a similar pocket is found in an already existing folder, it will stop.
        If some molecule is missing it will dock only the missing ligands.
    '''
    output_dir.mkdir(parents=True, exist_ok=True)
    protein_path = output_dir / 'protein.gro'

    # Initialize a fake text file
    string_io = io.StringIO()

    # Write the atoms group to the string_io    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)    
        with Writer(MDAnalysis.lib.util.NamedStream(string_io, filename=f'name.gro')) as writer:
            writer.write(atoms)
    protein_gro = string_io.getvalue().encode()

    # Save the protein for the pymol scripts of this structure.
    if not protein_path.is_file():
        protein_path.write_bytes(protein_gro)
    else:
        written_protein_gro = protein_path.read_bytes()
        if protein_gro != written_protein_gro:
            written_atoms = Universe(protein_path).atoms

            if written_atoms.n_atoms != atoms.n_atoms or np.any(np.abs(written_atoms.positions - atoms.positions) > 0.1):
                raise RuntimeError(f'The protein inside {output_dir} is different to the one passed by argument.')

    # Find the oriented box with the minimal volume
    origin, box_vectors = pocket.minimal_bounding_box()
    pocket_dir = search_similary_box(output_dir, origin, box_vectors)

    if pocket_dir:
        info_path = pocket_dir / 'info.json'
        info = json.loads(info_path.read_text())
    else:

        # No pockets found? Create a new one
        info = {
            'origin': origin.tolist(),
            'box_vectors': box_vectors.tolist(),
        }

        pocket_dir = Path(tempfile.mkdtemp(prefix='pocket_', dir=output_dir))
        info_path = pocket_dir / 'info.json'
        info_path.write_text(json.dumps(info, indent=1))

        # Write the pocket for the pymol script
        with gzip.open(pocket_dir / 'volume.pdb.gz', 'wt') as file:
            with PocketWriter(file) as pocket_writer:
                pocket_writer.write_frame([pocket])
                pocket_writer.write_pymol_script(pocket_dir / 'pymol.pml', 
                                                 structures=[protein_path], 
                                                 additional_strucutres={'molecules': pocket_dir / 'molecules.sdf.gz'})

    sdf_path = pocket_dir / 'molecules.sdf.gz'
    energies_path = Path(pocket_dir) / 'energies.pickle'


    if not energies_path.is_file():
        if only_cached:
            return None

        # Read the list of docked molecules
        try:
            with gzip.open(sdf_path, 'rb') as file:
                with SDFReader(file, sanitize=False) as sdf:
                    already_docked = list(sdf)
        except OSError:
            already_docked = set()
        except SystemError:
            return

        to_dock = smiles - set(m['smiles'] for m in already_docked)
    
        # Find remaining molecules to be docked
        if docking:
            gnina = Gnina(atoms, padding=3, n_poses=5, verbose=True, n_cpus=n_cpus, exhaustiveness=8,
                        box_origin=np.array(info['origin']), 
                        box_vectors=np.array(info['box_vectors']))

            # Write docked molecules    
            with SDFWriter(sdf_path) as sdf:

                # Write already docked molecules
                for m in already_docked:
                    sdf.write(m)

                # Dock molecules in chunks
                for chunk in chunk_generator(to_dock, 5):
                    docked = gnina.dock_list([Molecule.from_smiles(s, props={'smiles': s}) for s in chunk])

                    # Write the docked chunk
                    for poses in docked:
                        for m in poses:
                            sdf.write(m)

        if len(to_dock) < 20:
            energies = {}

            for m in already_docked:
                energies.setdefault(m['smiles'], []).append(float(m['CNNaffinity']))

            # import matplotlib.pyplot as plt
            # plt.hist([np.max(e) for e in energies.values()], bins=50)
            # plt.savefig(pocket_dir / 'energies.png', dpi=300)
            # plt.close()

            energies_path = Path(pocket_dir) / 'energies.pickle'
            energies_path.write_bytes(pickle.dumps(energies))
        else:
            return None
    else:
        energies = pickle.loads(energies_path.read_bytes())

    return energies

def search_similary_box(folder: Path, origin: np.ndarray, box_vectors: np.ndarray):
    
    # Create a grid of points inside the bounding_box 
    linspace = np.linspace(0.001, 0.999, 8, dtype=np.float32)
    unit_points = np.array(list(itertools.product(linspace, linspace, linspace)))
    box_points = origin + (unit_points @ box_vectors)
    box_volume = np.linalg.det(box_vectors)

    # Search for pockets already processed
    for pocket_dir in folder.iterdir():
        info_path = pocket_dir / 'info.json'

        if info_path.is_file():
            info = json.loads(info_path.read_text())

            # Check the similarity between the volume of  
            # this pocket and the one already performed
            if jaccard_similarity_boxes(box_points, box_volume, info) > 0.75:
                return pocket_dir
            
    return None

       
def jaccard_similarity_boxes(box_points: np.ndarray, box_volume: float, info: dict):
    '''
        Compute the jaccard similarity between the volumes of two pockets.
    '''
    pocket_origin      = np.array(info['origin'])
    pocket_box_vectors = np.array(info['box_vectors'])

    transformed = (box_points - pocket_origin) @ np.linalg.inv(pocket_box_vectors)

    points_inside = np.all((transformed >= 0) & (transformed <= 1), axis=1)
    intersection = box_volume * np.count_nonzero(points_inside) / box_points.shape[0]

    union = box_volume + np.linalg.det(pocket_box_vectors) - intersection
    return intersection / union