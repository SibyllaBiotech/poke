import os
from typing import TextIO
from pathlib import Path

from .pocket import Pocket

class PocketWriter:
    def __init__(self, output: TextIO | os.PathLike, resid_start=1) -> None:
        if isinstance(output, (str, Path)):
            self.output = open(output, 'wt')
            self.close_stream = True
        else:
            self.output = output
            self.close_stream = False

        self.output.write('TITLE     POKE POCKETS\n')
        self.output.write('CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n')
        self.model_idx = 1
        self.atom_idx = 1
        self.resid_idx = resid_start

    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.close_stream:   
            self.output.__exit__(exc_type, exc_val, exc_tb)

    def write_frame(self, pockets: list[Pocket]):
        pockets = [pockets] if isinstance(pockets, Pocket) else pockets
    
        self.output.write(f'MODEL {self.model_idx:>8}\n') 
        atom_idx = 0
        for p in pockets:
            for (x,y,z), e in zip(p.volume_positions, p.electro):
                self.output.write(f'ATOM  {atom_idx:5d}  A   PCK  {self.resid_idx:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00{e:6.2f}\n')
                atom_idx += 1
            self.resid_idx +=1
        self.output.write(f'ENDMDL\n')
        self.model_idx += 1
            
    def write_pymol_script(self, output: os.PathLike, structures: list[os.PathLike]=[], additional_strucutres: dict[str, Path]={}, pockets: os.PathLike=None,  scouted_output_indexes: list[int]=[],  extra_output_indexes: list[int]=[]) :
        output = Path(output)

        if pockets is None:
            try:
                pockets = self.output.name
            except:
                raise ValueError('The output stream has not a name. You must specify the "pockets" '
                                 'parameter with the path of the written pockets.')

        structure_lines = []

        for fn in structures:
            if len(structure_lines) == 0:
                structure_lines.append(f'cmd.load(os.path.join(script_path, "{os.path.relpath(fn, output.parent)}"), "protein")')
            elif len(structure_lines) == 1:
                structure_lines.append(f'cmd.load_traj(os.path.join(script_path, "{os.path.relpath(fn, output.parent)}"), "protein", state=1)')
            else:
                structure_lines.append(f'cmd.load_traj(os.path.join(script_path, "{os.path.relpath(fn, output.parent)}"), "protein", state=0)')

        for name, fn in additional_strucutres.items():
            structure_lines.append(f'cmd.load(os.path.join(script_path, "{os.path.relpath(fn, output.parent)}"), "{name}")')
            
        output.write_text('\n'.join([
            'import os',
            '',
            'script_path = os.path.abspath(os.path.dirname(pymol.__script__))',
            '',
            'set async_builds, 1',
            'set defer_builds_mode, 3',
            '',
            *structure_lines,
            '',
            f'create scouted, model protein and rank {condense_list(scouted_output_indexes)}' if len(scouted_output_indexes) else '',
            'create surf, scouted',
            'disable surf',
            'hide everything, surf',
            'util.cbaw surf',
            'show surface, surf',
            '',
            'set sphere_scale,0.1',
            'set connect_mode,1',
            f'cmd.load(os.path.join(script_path, "{os.path.relpath(pockets, output.parent)}"), "pockets")',
            'spectrum b, blue_white_red, pocket, minimum=-1, maximum=+1',
            'set connect_mode,0',
            '',
            f'create extra, model protein and rank {condense_list(extra_output_indexes)}' if len(extra_output_indexes) else '',
            'delete protein'
        ]))


def condense_list(indexes: list[int]):
    sorted_i = sorted(indexes)
    ranges = []
    first = sorted_i[0]
    last = sorted_i[0]
    for i in sorted_i[1:]:
        if i == last + 1:
            last = i
        else:
            if first == last:
                ranges.append(str(first))
            else:
                ranges.append(f'{first}-{last}')
            first = i
            last = i
    if first == last:
        ranges.append(str(first))
    else:
        ranges.append(f'{first}-{last}')
    
    return '+'.join(ranges)