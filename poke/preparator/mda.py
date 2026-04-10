import inspect
import numpy as np
import warnings
from MDAnalysis import Universe
from MDAnalysis.core.topology import Topology
from MDAnalysis.core.topologyattrs import AtomAttr
from MDAnalysis.coordinates.base import ProtoReader
from MDAnalysis.topology.guessers import guess_types, guess_masses

from . import Preparator

class PreparatorReader(ProtoReader):
    """A wrapper around an MDAnalysis trajectory reader that applies a Preparator to the coordinates."""

    def __init__(
        self, 
        reader, 
        preparator, 
        n_atoms, 
        input_preserved_ix, 
        input_transform_ix, 
        out_preserved_ix, 
        out_transform_ix, 
        **kwargs
    ):
        """Initialize the PreparatorReader.

        Args:
            reader (ProtoReader): The original MDAnalysis trajectory reader.
            preparator (Preparator): The Preparator to apply to the coordinates.
            n_atoms (int): The number of atoms in the output trajectory.
            input_preserved_ix (np.ndarray): Indices of atoms to preserve from the input.
            input_transform_ix (np.ndarray): Indices of atoms to transform from the input.
            out_preserved_ix (np.ndarray): Indices of atoms to preserve in the output.
            out_transform_ix (np.ndarray): Indices of atoms to transform in the output.
            **kwargs: Additional keyword arguments to pass to the Timestep constructor.
        """
              
        super().__init__()
        self._reader = reader
        self._preparator = preparator
        self._input_preserved_ix = input_preserved_ix
        self._input_transform_ix = input_transform_ix
        self._out_preserved_ix = out_preserved_ix
        self._out_transform_ix = out_transform_ix
        self._kwargs = kwargs
        self.n_atoms = n_atoms
        self.ts = self._Timestep(self.n_atoms, **kwargs)
        self.ts.has_velocities = self._reader.ts.has_velocities
        self.ts.has_forces = self._reader.ts.has_forces
        self.ts.data = reader.ts.data
        self._read_frame(reader.ts.frame)

    def __getattr__(self, name):
        """Delegate attribute access to the underlying reader."""
        return getattr(self._reader, name)
    
    def _update_ts(self, orig_ts):
        """Populate own Timestep with data from the original ts, applying the preparator as needed."""
        self.ts.frame = orig_ts.frame
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            self.ts.dt = orig_ts.dt
            self.ts.time = orig_ts.time
        self.ts.dimensions = orig_ts.dimensions
        self.ts._pos[self._out_preserved_ix] = orig_ts._pos[self._input_preserved_ix]
        self.ts._pos[self._out_transform_ix] = self._preparator(orig_ts._pos[self._input_transform_ix])
        if orig_ts.has_velocities:
            self.ts._velocities[self._out_preserved_ix] = orig_ts._velocities[self._input_preserved_ix]
            self.ts._velocities[self._out_transform_ix] = self._preparator.reshape(
                orig_ts._velocities[self._input_transform_ix], fill_value=0.0
            )
        if orig_ts.has_forces:
            self.ts._forces[self._out_preserved_ix] = orig_ts._forces[self._input_preserved_ix]
            self.ts._forces[self._out_transform_ix] = self._preparator.reshape(
                orig_ts._forces[self._input_transform_ix], fill_value=0.0
            )
        return self.ts
        
    
    def _read_next_timestep(self, ts=None):
        orig_ts = self._reader._read_next_timestep(ts)
        return self._update_ts(orig_ts)
    
    def _read_frame(self, frame):
        orig_ts = self._reader._read_frame(frame)
        return self._update_ts(orig_ts)

    def _reopen(self):
        self.ts.frame = -1
        self.ts.time = -1
        self._reader._reopen()

    def __del__(self):
        if hasattr(self._reader, '__del__'):
            self._reader.__del__()

    def copy(self):
        """Return a copy of the reader."""
        return self.__class__(self._reader.copy(), self._preparator, self.n_atoms, self._input_preserved_ix, 
                              self._input_transform_ix, self._out_preserved_ix, self._out_transform_ix, **self._kwargs)
        
        
class PreparatorTopology(Topology):
    """A subclass of MDAnalysis Topology to indicate it has been processed by a Preparator."""
    pass


class MDAUniverse(Universe):
    """
    A subclass of MDAnalysis Universe that includes a Preparator.
    Modifications are immediately reflected in the universe atoms at object creation. 
    Coordinates are transformed on-the-fly when reading trajectory frames.
    """

    # Deafults to populate the attributes of new atoms created by the Preparator
    _attr_default = {
        'elements': guess_types,
        'types': guess_types,
        'masses': lambda names: guess_masses(guess_types(names)),
    }

    def __init__(
        self, 
        topology=None, 
        *coordinates, 
        selection=None, 
        forcefield=None, 
        **kwargs
    ):
        """Initialize the MDAUniverse.
        Args:
            topology: The topology file or object.
            *coordinates: The coordinate files.
            selection (str, optional): The atom selection string. If None, a default scoutable selection is used.
            forcefield (str, optional): The forcefield to use for the Preparator. 
                Required if topology is not provided as a PreparatorTopology.
            **kwargs: Additional keyword arguments to pass to the Universe and Preparator constructors.
        """
        
        if isinstance(topology, PreparatorTopology):
            super().__init__(topology, **kwargs)
        
        else:
            if forcefield is None:
                raise ValueError("Forcefield must be provided for MDAUniverse initialization.")
            
            # Separate kwargs for different components
            prep_kwargs = {k: v for k, v in kwargs.items() if k in inspect.signature(Preparator).parameters}
            caps_kwargs = {k: v for k, v in kwargs.items() if k in inspect.signature(Preparator.get_caps).parameters}
            sres_kwargs = {k: v for k, v in kwargs.items() if k in inspect.signature(Preparator.select_restypes).parameters}
            skip_kwargs = set(prep_kwargs) | set(caps_kwargs) | set(sres_kwargs)
            mdau_kwargs = {k: v for k, v in kwargs.items() if k not in skip_kwargs}

            original_universe = Universe(topology, *coordinates, **mdau_kwargs)
            if selection is None:
                selection = Preparator.scoutable_selection(original_universe.residues.resnames)
            
            selected = original_universe.select_atoms(selection)
            input_transform_ix = selected.atoms.ix
            input_preserved_ix = np.setdiff1d(original_universe.atoms.ix, input_transform_ix)

            resnames = selected.residues.resnames
            resids = selected.residues.resids
            segids = selected.residues.segids
            atom_names = selected.residues.names
            
            if hasattr(selected.residues, 'altLocs'):
                alt_locs = selected.residues.altLocs
                tempfactors = getattr(selected.residues, 'tempfactors', [[0] * len(al) for al in alt_locs])
                occupancies = getattr(selected.residues, 'occupancies', [[0] * len(al) for al in alt_locs])
                
                atom_names, _, _, to_keep = Preparator.select_altlocations(atom_names, alt_locs, tempfactors, occupancies)
                input_transform_ix = input_transform_ix[to_keep]
            
            positions = original_universe.atoms.positions.astype(np.float32)[input_transform_ix]
            caps = Preparator.get_caps(resnames, resids, segids, **caps_kwargs)
            resnames = Preparator.select_restypes(resnames, atom_names, positions=positions, **sres_kwargs)
            preparator = Preparator(forcefield, resnames, resids, atom_names, caps, **prep_kwargs)

            # Determine the location of preserved and transformed atoms in the output array of positions
            # The original residue order is preserved, with transformed residues expanded in place
            out_ix = 0
            transformed_res_ix = 0
            out_transform_ix = []
            out_preserved_ix = []
            for res in original_universe.residues:
                n_preserved = np.sum(np.isin(res.atoms.ix, input_preserved_ix))
                out_preserved_ix.extend(range(out_ix, out_ix + n_preserved))
                out_ix += n_preserved
                if n_preserved < res.atoms.n_atoms:
                    n_transformed = len(preparator.atom_names[transformed_res_ix])
                    out_transform_ix.extend(range(out_ix, out_ix + n_transformed))
                    out_ix += n_transformed
                    transformed_res_ix += 1
            
            out_transform_ix = np.array(out_transform_ix, dtype=int)
            out_preserved_ix = np.array(out_preserved_ix, dtype=int)
            output_len = out_ix
    
            # Copy over attributes from the original topology, texpanding the atom attributes as needed
            base_attrs = set(dir(Topology(attrs=[])))
            conn_attrs = set(['bonds', 'angles', 'dihedrals', 'impropers']) # skip connectivity attributes for now
            curr_attrs = set(dir(original_universe._topology))

            attrs = []
            for attr_name in curr_attrs - base_attrs - conn_attrs:
                attribute = getattr(original_universe._topology, attr_name)
                attr_type = type(attribute)
                if issubclass(attr_type, AtomAttr):
                    new_values = np.empty(output_len, dtype=attribute.values.dtype)
                    new_values[out_preserved_ix] = attribute.values[input_preserved_ix]
                    if attr_name == 'names':
                        new_values[out_transform_ix] = preparator.atom_names_flat
                    else:
                        # If a default function is provided for this attribute, use it to populate new values
                        # Otherwise, pass None to fill with value of previous original atom
                        new_values[out_transform_ix] = preparator.reshape(
                            attribute.values[input_transform_ix], fill_value=self._attr_default.get(attr_name, None)
                        )
                    attribute = attr_type(new_values)
                attrs.append(attribute)

            # Assign atoms to residues based on original universe. For the new atoms, passing None to reshape uses the 
            # residue index of the previous original atom, to which the new atoms likely belong.
            new_residx = np.empty(output_len, dtype=original_universe.atoms.resindices.dtype)
            new_residx[out_preserved_ix] = original_universe.atoms.resindices[input_preserved_ix]
            new_residx[out_transform_ix] = preparator.reshape(original_universe.atoms.resindices[input_transform_ix])
            top = PreparatorTopology(
                output_len, 
                original_universe.residues.n_residues, 
                original_universe.segments.n_segments,
                attrs=attrs,
                atom_resindex=new_residx,
                residue_segindex=original_universe.residues.segindices
            )
  
            super().__init__(top, **kwargs)
            self.trajectory = PreparatorReader(
                original_universe.trajectory, 
                preparator, 
                output_len, 
                input_preserved_ix, 
                input_transform_ix,
                out_preserved_ix, 
                out_transform_ix, 
                **self.kwargs
            )

        self.selection = selection
        self._kwargs.update({
            'forcefield': forcefield,
            'selection': selection,
        })

        

