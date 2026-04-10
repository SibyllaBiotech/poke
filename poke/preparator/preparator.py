import json
import warnings
import numpy as np
from numba import njit
from pathlib import Path
from functools import cache
from typing import Any, Callable

from .newpos import sp2_1h, sp3_lone_1h, sp2_2h, sp3_3h, sp3_1h, sp3_2h, sp3_2o, sp3_ooh
from .special import is_hbonded, is_covbonded
from .species import get_species, COFACTOR_RESIDUES, STRUCTURAL_ION_RESIDUES


class Preparator:
    """Object used to prepare a structure for pocket scouting.

    For the moment this class adds missing hydrogens, caps protein and nucleic chains and selects best alt locations.
    
    The class computes the modifications needed to be applied to the structure during initialization, 
    then provides a callable interface to apply these modifications to the positions of the atoms.
    The application of modifications is fast and can be used in loops over trajectories.
    """

    class StrictLevel:
        NONE = 0
        WARNING = 1
        ERROR = 2


    conversions = {
        'ILE': {
            'CD1': 'CD',
        }
    }

    altnames = {
        'HEM': 'HEME',
        'NMA': 'NME',
    }

    

    def __init__(
        self, 
        forcefield: str | Path,
        resnames: list[str], 
        resids: list[int], 
        atom_names: list[list[str]], 
        caps: dict[int, str | list[str]], 
        keep_existing: bool=True, 
        strictness: int=StrictLevel.WARNING
    ) -> None:
        """Initialize the Preparator object.
        
        Args:
            forcefield: Path to the forcefield database.
            resnames: List of residue names.
            resids: List of residue IDs.
            atom_names: List of lists of atom names for each residue.
            caps: Dictionary mapping residue indices to cap types.
            keep_existing: Whether to keep existing hydrogens or not.
            strictness: Strictness level for missing data handling.
        """
        
        self.hydrogens_db, self.termini_db = cached_loader(forcefield)
        
        self.resnames = resnames.copy()
        self.atom_names = []
        atom_names_flat = []
        mods = []
        is_orig = [False]
        load_ix = [0]
        outp_ix = []
        newn_ix = []
        orig_ix = []
        res_atom_names_ix_map = []
        res_mods_to_do = []
        input_absolute_ix = 0
        absolute_ix = 1  # 0 is reserved for discarding generated hydrogens

        # Collect all modifications to be applied. We perform a first pass to collect all modifications before 
        # evaluating each one, since some modifications can depend on others.
        for resix, (resname, res_atom_names) in enumerate(zip(resnames, atom_names)):
            if resname in Preparator.altnames:
                resname = Preparator.altnames[resname]
            
            if resname not in self.hydrogens_db:
                # Residue not found in hydrogens database, skip or raise error
                if strictness == Preparator.StrictLevel.ERROR:
                    raise RuntimeError(f'No hydrogens data for {resname}. Add it to {forcefield}/hydrogens.json')
                
                if strictness == Preparator.StrictLevel.WARNING:
                    warnings.warn(f'No hydrogens data for {resname}. Add it to {forcefield}/hydrogens.json')
                
                additions = []
                deletions = set()

            else:
                # Residue found in hydrogens database, get modifications
                additions = self.hydrogens_db[resname]
                deletions = set()
                
                if resix in caps:
                    cap = caps[resix]
                    if isinstance(cap, str):
                        cap = [cap]
                    for c in cap:
                        if c not in self.termini_db:
                            raise RuntimeError(f'No termini data for {c}. Add it to {forcefield}/termini.json')
                        cap_mod = self.termini_db[c].get(resname, self.termini_db[c]['default'])
                        additions = [*additions, *cap_mod['add']]
                        deletions.update(cap_mod['delete'])

            # Create a mapping from atom names to indeces of the internal memory array that will contain
            # positions. This mapping is used to apply modifications later. Also populate the load_ix array
            # used to insert positions into the memory.
            atom_names_ix_map = {}
            for atom_name in res_atom_names:
                if atom_name[0].isdigit():
                    atom_name = atom_name[1:] + atom_name[0]
                if resname in Preparator.conversions and atom_name in Preparator.conversions[resname]:
                    atom_name = Preparator.conversions[resname][atom_name]
                if keep_existing or not atom_name.lower().startswith('h'):
                    atom_names_ix_map[atom_name] = absolute_ix
                    load_ix.append(input_absolute_ix)
                    is_orig.append(True)
                    absolute_ix += 1
                input_absolute_ix += 1    
            
            # Collect addtions to do, assigning a mapping to the new atom names 
            additions_to_do = []
            for add in additions:
                new_names = add['names']
                if not set(new_names) <= atom_names_ix_map.keys():
                    for new_name in new_names:
                        load_ix.append(input_absolute_ix-1)
                        is_orig.append(False)
                        atom_names_ix_map[new_name] = absolute_ix
                        absolute_ix += 1
                    
                    additions_to_do.append(add)
            
            res_atom_names_ix_map.append(atom_names_ix_map)
            res_mods_to_do.append((additions_to_do, deletions))

        
        # Perform a second pass to create the final array of modifications to be applied.
        # Also populate the final array of output atom names
        # TODO: The modifications array should be ordered to correclty account for dependencies. For now, it is
        # correct only for modifications depending on previous residues modifications, and on current residue if the
        # correct ordering is given in the input force-field modifications database file.
        for resix, (additions_to_do, deletions) in enumerate(res_mods_to_do):
            for addition in additions_to_do:
                mod = [addition['type'], 0, 0, 0, 0, 0, 0, 0]
                for ai, (anchor_name, anchor_rel_resix) in enumerate(zip(addition['anchors'], addition['resindex'])):
                    anchor_resix = resix + anchor_rel_resix
                    if (anchor_resix < 0 or anchor_resix >= len(res_atom_names_ix_map) 
                            or (resids[anchor_resix] - resids[resix]) != anchor_rel_resix 
                            or anchor_name not in res_atom_names_ix_map[anchor_resix]):
                        # Anchor atom not found, skip this modification and remove associated new names
                        deletions.update(addition['names'])
                        break
                    # Translate anchor name to the index in the memory array where it is stored
                    mod[ai + 1] = res_atom_names_ix_map[anchor_resix][anchor_name]
                else:
                    # All anchors were found, so we can add the modification
                    for ni, new_name in enumerate(addition['names'], start=5):
                        mod[ni] = res_atom_names_ix_map[resix][new_name]
                    mods.append(mod)

            # Collect output atom names for the residue, excluding deleted atoms
            out_atom_names = []
            for atom_name, atom_ix in res_atom_names_ix_map[resix].items():
                if atom_name not in deletions:
                    out_atom_names.append(atom_name)
                    atom_names_flat.append(atom_name)
                    outp_ix.append(atom_ix)
                    if is_orig[atom_ix]:
                        orig_ix.append(atom_ix)
                    else:
                        newn_ix.append(atom_ix)
            self.atom_names.append(out_atom_names)
        
        self.mods = np.array(mods, dtype=np.int32, ndmin=2)
        self.atom_names_flat = np.array(atom_names_flat)
        self.memory = np.empty((len(load_ix), 3), dtype=np.float32)
        self.load_ix = np.array(load_ix, dtype=np.int32)
        self.outp_ix = np.array(outp_ix, dtype=np.int32)
        self.newn_ix = np.array(newn_ix, dtype=np.int32)
        self.orig_ix = np.array(orig_ix, dtype=np.int32)
            
    
    def reshape(self, array: np.ndarray | list, fill_value: Any | Callable[[np.ndarray], Any] = None) -> np.ndarray:
        """Reshape the input array to match the output positions length.
        
        Args:
            array: Input array to reshape.
            fill_value: Value to fill in the new positions. 
                If a callable is provided, it will be called with the new atom names to generate the fill value.
        
        Returns:
            Reshaped array with the same length as the output positions.
        """

        array = np.array(array, copy=False)
        tmp = array[self.load_ix]
        if fill_value is not None:
            if callable(fill_value):
                new_names = self.atom_names_flat[np.where(np.isin(self.outp_ix, self.newn_ix))[0]]
                fill_value = fill_value(new_names)
            tmp[self.newn_ix] = fill_value
        return tmp[self.outp_ix]


    @staticmethod
    def scoutable_selection(resnames: list[str] | set[str]) -> str:
        """Generate an MDAnalysis selection string for scouting. 
        The selection contains protein, nucleic, cofactors and structural ions.
        
        Args:
            resnames: List of residue names.
        
        Returns:
            Selection string for scouting.
        """
        resnames = set(resnames)
        selection = 'protein or resname CSO or nucleic'
        cofactor_res = resnames & COFACTOR_RESIDUES
        structural_ions = resnames & STRUCTURAL_ION_RESIDUES
        if cofactor_res:
            selection += ' or resname ' + ' or resname '.join(cofactor_res)
        if structural_ions:
            selection += ' or resname ' + ' or resname '.join(structural_ions)
        return selection

    
    @staticmethod
    def select_altlocations(
        atom_names: list[list[str]], 
        alt_locs: list[list[str]], 
        temp_factors: list[list[float]], 
        occupancies: list[list[float]]
    ) -> tuple[list[list[str]], list[list[float]], list[list[float]], np.ndarray]:
        """Select the best alt locations for each residue. 

        Selection is based on highest occupancy and, in a tie, lowest temperature factor of the atoms.
        
        Args:
            atom_names: List of lists of atom names for each residue.
            alt_locs: List of lists of alternative location letter for each atom.
            temp_factors: List of lists of temperature factors for each atom.
            occupancies: List of lists of occupancies for each atom.
        
        Returns:
            Tuple containing:
                - List of lists of selected atom names for each residue.
                - List of lists of selected temperature factors for each residue.
                - List of lists of selected occupancies for each residue.
                - Numpy array of indices to keep in the original flattened atom names list.
        """
        out_atom_names = []
        out_temp_factors = []
        out_occupancies = []
        to_keep = []
        ai = 0
        for res_atom_names, res_alt_locs, res_temp_factors, res_occupancies in zip(atom_names, alt_locs, temp_factors, 
                                                                                   occupancies):
            at_dict = {}
            br_dict = {}
            for atom_name, alt_loc, temp_factor, occupancy in zip(res_atom_names, res_alt_locs, res_temp_factors, 
                                                                  res_occupancies):
                at_dict.setdefault(atom_name, {})[alt_loc] = (ai, temp_factor, occupancy)
                if alt_loc != '':
                    cur = br_dict.get(alt_loc, (0,0))
                    br_dict[alt_loc] = (cur[0]+occupancy, cur[1]-temp_factor)
                ai += 1
            best_branch = max(br_dict, key=br_dict.__getitem__) if br_dict else None
            
            res_out_atom_names = []
            res_out_temp_factors = []
            res_out_occupancies = []
            for atom_name, alt in at_dict.items():
                ai_sel, temp_factor, occupancy = alt.popitem()[1] if len(alt) == 1 else alt[best_branch]
                to_keep.append(ai_sel)
                res_out_atom_names.append(atom_name)
                res_out_temp_factors.append(temp_factor)
                res_out_occupancies.append(occupancy)
            out_atom_names.append(res_out_atom_names)
            out_temp_factors.append(res_out_temp_factors)
            out_occupancies.append(res_out_occupancies)
     
        return out_atom_names, out_temp_factors, out_occupancies, np.array(to_keep, dtype=np.int32)


    @staticmethod
    def get_caps(
        resnames: list[str], 
        resids: list[int], 
        segids: list[str], 
        caps_types: dict[str,tuple[str,str]] = {'protein': ('NH3+', 'COO-'), 'nucleic_acid': ('5TER', '3TER')},
        max_interruption_size: int = 5,
        on_segid: bool = True, 
        resid_threshold: int = 50, 
        skip: set[str] | list[str] = {'ACE', 'NME'}, 
        ter: list[int] | set[int] = [] 
    ) -> dict[int, str | list[str]]:
        """Select the caps for the different chains.

        The caps are added for chains of the type specified in caps_types. This method tries to determine
        the chains separations based on different criteria. Separations occurs when:
            - The species type changes. In this case a maximum number of residues of type 'other' is allowed, given by
            max_interruption_size, before the chain is considered interrupted.
            - The segid changes, if on_segid is True.
            - There is jump in residue IDs of at least resid_threshold, only if segid is empty for residues across
            the jump.
            - At specific residue indices, given by ter.
        
        Args:
            resnames: List of residue names.
            resids: List of residue IDs.
            segids: List of segment IDs.
            caps_types: Dictionary mapping species types to cap types.
            max_interruption_size: Maximum number of 'other' residues allowed before interruption.
            on_segid: Whether to consider segid changes as chain separations.
            resid_threshold: Minimum jump in residue IDs to consider as chain separation.
            skip: Residue names to skip when adding caps.
            ter: List or set of residue indices where explicit termini should be added.
        
        Returns:
            Dictionary mapping residue indices to cap types or lists of cap types.
        """

        ter = set(ter)
        skip = set(skip)
        caps = {}
        chain_start = None
        chain_specs = None
        chain_segid = None
        inter_start = None
        for i in range(len(resnames)+1):
            species = get_species(resnames[i]) if i < len(resnames) else 'other'
            segid = segids[i] if i < len(resnames) else ''
            same_segid = (not on_segid or chain_segid == segid)
            explicit_term = i - 1 in ter
            resid_jump = (0 < i < len(resnames) and abs(resids[i] - resids[i-1]) > resid_threshold 
                          and segids[i].strip() == '' and segids[i-1].strip() == '')

            if chain_specs == species and same_segid and not explicit_term and not resid_jump:
                inter_start = None
            
            else:
                if chain_specs is not None:
                    assert chain_start is not None, 'chain_start should not be None if chain_specs is not None'
                    
                    if inter_start is None:
                        inter_start = i - 1

                    if (species != 'other' or i == len(resnames) or not same_segid or explicit_term or resid_jump 
                        or i - inter_start > max_interruption_size):
                        for teri, pos in enumerate((chain_start, inter_start)):
                            if resnames[pos] not in skip:
                                if pos in caps:
                                    if isinstance(caps[pos], list):
                                        caps[pos].append(caps_types[chain_specs][teri])
                                    else:
                                        caps[pos] = [caps[pos], caps_types[chain_specs][teri]]
                                    
                                else:
                                    caps[pos] = caps_types[chain_specs][teri]
                    
                        chain_start = None
                        chain_specs = None
                
                
                if species in caps_types:
                    chain_start = i
                    chain_specs = species   
                    chain_segid = segids[i]
                    inter_start = None 
        return caps
            
    
    @staticmethod
    def select_restypes(
        resnames: list[str], 
        atom_names: list[list[str]], 
        his: str='auto', 
        lys: str='auto', 
        arg: str='auto', 
        glu: str='auto', 
        asp: str='auto', 
        cys: str='auto', 
        specific: dict[int, str] | None = None, 
        positions: np.ndarray | None = None
    ) -> list[str]:
        """Select the best residue types for each residue, considering possible protonation states.

        For the standard residues with alternative protonation states, specific residue types can be provided, 
        or 'auto' for automatic detection. The automatic detection depends on the specific residue and 
        the presence of certain hydrogens in the residue.
        
        For histidine, the protonation state is determined based on the presence of HD1 and HE2 hydrogens. 
        If these cannot be found, it is determined based on hydrogen bonding detection. In this case,
        the positions of the atoms are required. To only use the hydrogen bonding detection, set his to 'bond'.

        For cysteine, the protonation state is determined based on the presence of HG1 hydrogen. 
        If this hydrogen is not found, it is determined based on covalent bonding detection. 
        In this case, the positions of the atoms are required. To only use the covalent bonding detection, 
        set cys to 'bond'.

        For lysine, arginine, glutamate, and aspartate, the protonation state is determined based on the presence of
        specific hydrogens. 

        Additionally, specific residue names can be provided for residues at specific indices.

        Args:
            resnames: List of residue names.
            atom_names: List of lists of atom names for each residue.
            his: Protonation state for histidine ('auto', 'bond', 'HSP', 'HSE', 'HSD').
            lys: Protonation state for lysine ('auto', 'LYS', 'LSN').
            arg: Protonation state for arginine ('auto', 'ARG', 'ARGN').
            glu: Protonation state for glutamate ('auto', 'GLU', 'GLUP').
            asp: Protonation state for aspartate ('auto', 'ASP', 'ASPP').
            cys: Protonation state for cysteine ('auto', 'bond', 'CYS', 'CYS2').
            specific: Dictionary mapping residue indices to specific residue names.
            positions: Numpy array of positions, required for some automatic protonation state detections.
        
        Returns:
            List of residue names with the best protonation states selected.
        """
        

        if specific is None:
            specific = {}

        out_resnames = []

        acceptor_positions = None
        covalent_bonding = None
        
        ai = 0
        for resi, (resname, res_atom_names) in enumerate(zip(resnames, atom_names)):
            if resi not in specific:
                if resname == 'HIS' or resname == 'HSP' or resname == 'HSE' or resname == 'HSD':
                    if his in ('auto', 'bond'):
                        res_a = {}
                        for atom_name in res_atom_names:
                            res_a[atom_name] = ai
                            ai += 1
                        if his == 'auto':
                            if 'HD1' in res_a and 'HE2' in res_a:
                                out_resnames.append('HSP')
                                continue
                            elif 'HD1' in res_a:
                                out_resnames.append('HSD')
                                continue
                            elif 'HE2' in res_a:
                                out_resnames.append('HSE')
                                continue
                        
                        # Check if we have covalent bonding (with HEME)
                        if covalent_bonding is None:
                            covalent_bonding = Preparator._get_covalent_bonding(resnames, atom_names, positions)
                        
                        if resi in covalent_bonding:
                            out_resnames.append('HSD')
                        else:

                            # Time for auto detection with h bond
                            if 'ND1' in res_a and 'CG' in res_a and 'CE1' in res_a:
                                if 'NE2' in res_a and 'CD2' in res_a:
                                    if acceptor_positions is None:
                                        acceptor_positions = Preparator._get_acceptor_positions(resnames, 
                                                                                                atom_names, 
                                                                                                positions)
                                    
                                    assert positions is not None, 'Positions must be provided for HIS detection'
                                    hd1 = sp2_1h(positions[res_a['ND1']], 
                                                 positions[res_a['CG']], 
                                                 positions[res_a['CE1']])
                                    he2 = sp2_1h(positions[res_a['NE2']], 
                                                 positions[res_a['CE1']], 
                                                 positions[res_a['CD2']])
                                    
                                    if is_hbonded(positions[res_a['ND1']], hd1, acceptor_positions):
                                        if is_hbonded(positions[res_a['NE2']], he2, acceptor_positions):
                                            out_resnames.append('HSP')
                                        else:
                                            out_resnames.append('HSD')
                                    else:
                                        out_resnames.append('HSE')
                                else:
                                    out_resnames.append('HSD')
                            else:
                                out_resnames.append('HSE')
                        
                    else:
                        out_resnames.append(his.upper())    
                        ai += len(res_atom_names)          
                
                elif resname == 'CYS' or resname == 'CYM' or resname == 'CYS2':
                    if cys in ('auto', 'bond'):
                        if cys == 'auto' and 'HG1' in res_atom_names:
                                out_resnames.append('CYS')
                         
                        else:
                            if covalent_bonding is None:
                                covalent_bonding = Preparator._get_covalent_bonding(resnames, atom_names, positions)
                            if resi in covalent_bonding:
                                out_resnames.append('CYS2')
                            else:
                                out_resnames.append('CYS')
                    
                    else:
                        out_resnames.append(cys.upper())
                
                    ai += len(res_atom_names)
                            

                elif resname == 'LYS' or resname == 'LSN':
                    if lys == 'auto':
                        if ('HZ1' in res_atom_names 
                            and 'HZ2' in res_atom_names 
                            and 'HZ3' not in res_atom_names):
                            out_resnames.append('LSN')
                        else:
                            out_resnames.append('LYS')
                    
                    else:
                        out_resnames.append(lys.upper())

                    ai += len(res_atom_names)

                elif resname == 'ARG' or resname == 'ARGN':
                    if arg == 'auto':
                        if 'HH11' in res_atom_names and 'HH12' not in res_atom_names:
                            out_resnames.append('ARGN')
                        else:
                            out_resnames.append('ARG')
                    
                    else:
                        out_resnames.append(arg.upper())

                    ai += len(res_atom_names)
                    
                elif resname == 'GLU' or resname == 'GLUP':
                    if glu == 'auto':
                        if 'HE2' in res_atom_names:
                            out_resnames.append('GLUP')
                        else:
                            out_resnames.append('GLU')
                    
                    else:
                        out_resnames.append(glu.upper())

                    ai += len(res_atom_names)

                elif resname == 'ASP' or resname == 'ASPP':
                    if asp == 'auto':
                        if 'HD2' in res_atom_names:
                            out_resnames.append('ASPP')
                        else:
                            out_resnames.append('ASP')
                    
                    else:
                        out_resnames.append(asp.upper())

                    ai += len(res_atom_names)
                
                else:
                    out_resnames.append(resname.upper())
                    ai += len(res_atom_names)

            else:
                out_resnames.append(specific[resi].upper())
                ai += len(res_atom_names)

        return out_resnames  
    
    
    def __call__(self, positions: np.ndarray) -> np.ndarray:
        """Apply the modifications to the positions and return the modified positions.
        
        Args:
            positions: Numpy array of positions to modify.
        
        Returns:
            Numpy array of modified positions.
        """
        return Preparator._apply(positions, self.mods, self.memory, self.load_ix, self.outp_ix)
        

    @staticmethod
    @njit(cache=True, fastmath=True)
    def _apply(
        positions: np.ndarray, 
        mods: np.ndarray, 
        memory: np.ndarray, 
        load_ix: np.ndarray, 
        outp_ix: np.ndarray
    ) -> np.ndarray:
        """Apply the modifications to the positions and return the modified positions. Numba-optimized version.
        
        Args:
            positions: Numpy array of positions to modify.
            mods: Numpy array of modifications to apply.
            memory: Numpy array to store intermediate results.
            load_ix: Numpy array of indices to load positions into memory.
            outp_ix: Numpy array of indices to output modified positions.
        
        Returns:
            Numpy array of modified positions.
        """
        memory[:] = positions[load_ix]
        if (mods.size): 
            for mod in mods:
                if mod[0] == 1:
                    memory[mod[5]] = sp2_1h(memory[mod[1]], memory[mod[2]], memory[mod[3]])
                elif mod[0] == 2:
                    memory[mod[5]] = sp3_lone_1h(memory[mod[1]], memory[mod[2]], memory[mod[3]])
                elif mod[0] == 3:
                    memory[mod[5:7]] = sp2_2h(memory[mod[1]], memory[mod[2]], memory[mod[3]])
                elif mod[0] == 4:
                    memory[mod[5:8]] = sp3_3h(memory[mod[1]], memory[mod[2]], memory[mod[3]])
                elif mod[0] == 5:
                    memory[mod[5]] = sp3_1h(memory[mod[1]], memory[mod[2]], memory[mod[3]], memory[mod[4]])
                elif mod[0] == 6:
                    memory[mod[5:7]] = sp3_2h(memory[mod[1]], memory[mod[2]], memory[mod[3]])
                elif mod[0] == 8:
                    memory[mod[5:7]] = sp3_2o(memory[mod[1]], memory[mod[2]], memory[mod[3]])
                elif mod[0] == 9:
                    memory[mod[5:8]] = sp3_ooh(memory[mod[1]], memory[mod[2]], memory[mod[3]])
                else:
                    raise RuntimeError(f'Unknown modification type {mod[0]}')
            
        return memory[outp_ix]

    @staticmethod
    def _query_atoms(
        resnames: list[str], 
        atom_names: list[list[str]], 
        query: list[str] | set[str] | dict[str, list[str] | set[str]]
    ) -> tuple[list[int],list[int]]:
        """
        Get indices of query atom names in the flattened atom_names list, and the indeces of the residues in which 
        they are located.
        
        Args:
            resnames: List of residue names.
            atom_names: List of lists of atom names for each residue.
            query: List or set of atom names to query, or a dictionary mapping atom names to query to sets of residue 
            names to restrict the search.
        
        Returns:
            Tuple containing:
                - List of residue indices where the query atom names are found.
                - List of indices in the flattened atom_names list where the query atom names are found.
        """
        
        res_ixs = []
        ixs = []
        ix = 0
        for resi, (resname, res_atom_names) in enumerate(zip(resnames, atom_names)):
            for atom_name in res_atom_names:
                if atom_name in query and (not isinstance(query, dict) or resname in query[atom_name]):
                    res_ixs.append(resi)
                    ixs.append(ix)
                ix += 1
        
        return res_ixs, ixs
    

    @staticmethod
    def _get_covalent_bonding(
        resnames: list[str], 
        atom_names: list[list[str]], 
        positions: np.ndarray | None
    ) -> dict[int, int]:
        """ Determine covalent bonding between residues, specifically for CYS and HIS residues.

        Args:
            resnames: List of residue names.
            atom_names: List of lists of atom names for each residue.
            positions: Numpy array of positions, required for covalent bonding detection.
        
        Returns:
            Dictionary mapping residue indices to indices of residues they are covalently bonded to.
        """
        
        if positions is None:
            raise RuntimeError(
                'Positions must be provided to automatically determine special covalent bonding '
                '(e.g., for CYS determination).'
            )
        
        bonded = {}
        cys_res_ixs, cys_ixs = Preparator._query_atoms(resnames, atom_names, {'SG': {'CYS', 'CYM', 'CYS2'}})
        his_res_ixs, his_ixs = Preparator._query_atoms(resnames, atom_names, {'NE2': {'HIS', 'HSP', 'HSD', 'HSE'}})
        hemef_res_ixs, hemef_ixs = Preparator._query_atoms(resnames, atom_names, {'FE': {'HEM', 'HEME'}})
        hemec_res_ixs, hemec_ixs = Preparator._query_atoms(resnames, atom_names, {'CAB': {'HEM', 'HEME'}, 
                                                                                  'CAC': {'HEM', 'HEME'}})

        for i in range(len(his_ixs)):
            if his_res_ixs[i] in bonded:
                continue
            j = is_covbonded(positions[his_ixs[i]], positions[hemef_ixs], cov_dist=2.0)
            if j >= 0 and hemef_res_ixs[j] not in bonded:
                bonded[his_res_ixs[i]] = hemef_res_ixs[j]
                bonded[hemef_res_ixs[j]] = his_res_ixs[i]
            

        for i in range(len(cys_ixs)):
            if cys_res_ixs[i] in bonded:
                continue
            j = is_covbonded(positions[cys_ixs[i]], positions[hemef_ixs], cov_dist=2.5)
            if j >= 0 and hemef_res_ixs[j] not in bonded:
                bonded[cys_res_ixs[i]] = hemef_res_ixs[j]
                bonded[hemef_res_ixs[j]] = cys_res_ixs[i]
                continue
            j = is_covbonded(positions[cys_ixs[i]], positions[hemec_ixs], cov_dist=1.8)
            if j >= 0 and hemec_res_ixs[j] not in bonded:
                bonded[cys_res_ixs[i]] = hemec_res_ixs[j]
                bonded[hemec_res_ixs[j]] = cys_res_ixs[i]
                continue
            j_rel = is_covbonded(positions[cys_ixs[i]], positions[cys_ixs[i+1:]], cov_dist=2.0)
            if j_rel >= 0:
                j = j_rel + i + 1
                bonded[cys_res_ixs[i]] = cys_res_ixs[j]
                bonded[cys_res_ixs[j]] = cys_res_ixs[i]
        return bonded   
    

    @staticmethod
    def _get_acceptor_positions(
        resnames: list[str], 
        atom_names: list[list[str]], 
        positions: np.ndarray | None
    ) -> np.ndarray:
        """
        Get positions of hydrogen acceptor atoms.

        Args:
            resnames: List of residue names.
            atom_names: List of lists of atom names for each residue.
            positions: Numpy array of positions, required for hydrogen bonding detection.
        
        Returns:
            Numpy array of positions of hydrogen acceptor atoms.
        """
        
        if positions is None:
            raise RuntimeError('Positions must be provided to automatically determine hydrogen bonding (e.g., for HIS determination).')

        acceptor = { 'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'OW' }
        _, ixs = Preparator._query_atoms(resnames, atom_names, acceptor)
        return positions[ixs]



@cache
def cached_loader(ff: Path):
    hydrogens_db = json.loads(Path(ff / 'hydrogens.json').read_bytes())
    termini_db = json.loads(Path(ff / 'termini.json').read_bytes())
    return hydrogens_db, termini_db


