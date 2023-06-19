import os
import io
import gzip
import numpy as np

from pathlib import Path

from meeko import MoleculePreparation
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

class Molecule:

    def __init__(self, rdkit: Chem.Mol, props: dict[str,str]=None):
        self.rdkit = rdkit

        if props is not None:
            self.update(props)
    
    @classmethod
    def from_rdkit(cls, rdkit: Chem.Mol, props: dict[str,str]=None):
        return cls(Chem.Mol(rdkit), props=props)

    @classmethod
    def from_smiles(cls, smiles: str, props: dict[str,str]=None):
        return cls(Chem.MolFromSmiles(smiles), props=props)
    
    @classmethod
    def from_bytes(cls, data: bytes, props: dict[str,str]=None):
        return cls(Chem.Mol(data), props=props)

    @classmethod
    def from_sdf(cls, sdf: os.PathLike, props: dict[str,str]=None):
        path = Path(sdf)

        if path.suffix == '.gz':
            supplier = Chem.ForwardSDMolSupplier(gzip.open(path), removeHs=False)
        else:
            supplier = Chem.SDMolSupplier(path.as_posix(), removeHs=False)

        for rdkit_mol in supplier:
            if rdkit_mol:
                yield cls(rdkit_mol, props=props)



    # Props managment ----------------------
    
    @property
    def id(self):
        return self.rdkit.GetProp('_Name')

    @id.setter
    def id(self, identifier):
        self['_Name'] = identifier

    @property
    def props(self) -> dict[str,str]:
        return self.rdkit.GetPropsAsDict()
    
    @props.setter
    def props(self, new: dict[str,str]):
        for key in self.rdkit.GetPropNames():
            if key not in new:
                self.rdkit.ClearProp(key)
        self.update(new)

    def update(self, props: dict[str,str]):
        for key, value in props.items():
            self.rdkit.SetProp(key, value)

    def __contains__(self, key: str) -> bool:
        return self.rdkit.HasProp(key)

    def __getitem__(self, key: str) -> str:
        return self.rdkit.GetProp(key)
  
    def __setitem__(self, key: str, value: str):
        assert isinstance(value, str), f'The property {key} must be a string, not a {type(value).__name__}'
        self.rdkit.SetProp(key, value)

    def __delitem__(self, key: str):
        self.rdkit.ClearProp(key)

    def get(self, key: str, default):
        return self.rdkit.GetProp(key) if self.rdkit.HasProp(key) else default

    # Serializations ----------------------

    @property
    def bytes(self) -> bytes:
        return self.rdkit.ToBinary(propertyFlags=True)

    @property
    def smiles(self) -> str:
        return Chem.MolToSmiles(self.rdkit)

    @property
    def inchi(self) -> str:
        return Chem.MolToInchi(self.rdkit)

    @property
    def inchi_key(self) -> str:
        return Chem.MolToInchiKey(self.rdkit)

    @property
    def pdbqt(self) -> str:
        # macrocycles flexible by default since v0.3.0
        preparator = MoleculePreparation(hydrate=False, flexible_amides=True) 
        preparator.prepare(self.rdkit)
        return preparator.write_pdbqt_string()
    
    @property
    def sdf(self) -> str:
        buffer = io.StringIO()
        with SDFWriter(buffer) as writer:
            writer.write(self)
        return buffer.getvalue()

    def __str__(self) -> str:
        return f'{self["id"]}' if 'id' in self else self.smiles

    def __repr__(self) -> str:
        return f'Molecule<{self.__str__()}>'
        
    def __reduce__(self) -> tuple:
        return self.__class__.from_bytes, (self.bytes,)
    

    # Chemical modifications--------------------

    def generate_3d_coordinates(self, mmffVariant='MMFF94s', maxIters=1500) -> 'Molecule':
        rdkit = Chem.Mol(self.rdkit)
        AllChem.MMFFOptimizeMolecule(rdkit, mmffVariant=mmffVariant, maxIters=maxIters)
        return Molecule(rdkit)

    def add_hydrogens(self) -> 'Molecule':
        return Molecule(Chem.AddHs(self.rdkit))

    def apply_matrix(self, matrix: np.ndarray):
        
        # If the matrix is 3x3 create a matrix 4x4
        # compatible with rdkit
        if matrix.shape[0] == 3:
            new_matrix = np.eye(4)
            new_matrix[:3,:3] = matrix
            matrix = new_matrix
        
        rdMolTransforms.TransformConformer(self.rdkit.GetConformer(0), matrix)
        return self

class SDFWriter:
    def __init__(self, out: os.PathLike) -> None:
        if isinstance(out, Path):
            out = out.as_posix()
            
        if out.endswith('.gz'):
            self.stream = gzip.open(out, 'wt')
            self.writer = Chem.SDWriter(self.stream)
        else:
            self.stream = None
            self.writer = Chem.SDWriter(out)

        self.context = self.writer

    def __enter__(self):
        self.context = self.writer.__enter__()
        return self
    
    def __exit__(self, type, value, traceback):
        self.writer.__exit__(type, value, traceback)
        
        if self.stream is not None:
            self.stream.__exit__(type, value, traceback)

    def write(self, molecule: Molecule):
        return self.context.write(molecule.rdkit)

    def flush(self):
        self.context.flush()

class SDFReader:
    def __init__(self, input: os.PathLike, sanitize=True, removeHs=False) -> None:
        if isinstance(input, Path):
            input = input.as_posix()
        if isinstance(input, io.IOBase):
            self.reader = Chem.ForwardSDMolSupplier(input, removeHs=removeHs, sanitize=sanitize)
        else:
            self.reader = Chem.SDMolSupplier(input, removeHs=removeHs, sanitize=sanitize)
        self.context = self.reader

    def __enter__(self):
        self.context = self.reader.__enter__()
        return self
    
    def __exit__(self, type, value, traceback):
        self.reader.__exit__(type, value, traceback)

    def __iter__(self):
        for rdmol in self.context:
            yield Molecule(rdmol)