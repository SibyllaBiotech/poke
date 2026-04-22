![poke](./img/poke-banner.svg)
[![PyPI](https://img.shields.io/pypi/v/poke)](https://pypi.org/project/poke/)
[![License](https://img.shields.io/pypi/l/poke)](https://github.com/SibyllaBiotech/poke/blob/main/LICENSE)
[![Python](https://img.shields.io/pypi/pyversions/poke)](https://pypi.org/project/poke/)

# Overview
Poké is a high-performance pocket scouting engine. It uses an energy-based approach to identify low-energy regions on protein surfaces. These regions are characterized by a set of physical descriptors, which are used by a neural network to assign a ligandability score to the pocket.

It is designed to be used on long molecular dynamics trajectories, but it works equally well on single frames.


## Installation
Poké is available as a wheel package for macOS, Windows and Linux on PyPI. Install it using `pip`:
```console
python -m pip install poke
```
Alternatively, the package can be installed directly from GitHub:
```console
git clone https://github.com/SibyllaBiotech/poke
pip install ./poke
```

## Getting started
> [!CAUTION]
> All input structures passed to *Poké* must be protonated and the atom names should match the CHARMM or Amber convention. Alternatively, use the `--prepare` flag to let *Poké* handle protonation automatically (see [Structure Preparation](#structure-preparation)).

*Poké* can be used on a single structure with the following command:
```console
poke -dvp protein.pdb
```

The result will appear in a folder named after the input file (in this example `protein`) in the current working directory.
The flags `-dvp` specify the desired output files:
- `-d` for the descriptors CSV
- `-v` for the volumes (a PDB file with the grid points of each pocket)
- `-p` for a PyMOL script to visualize the pockets

Use `-h` to see all available options, or see the [CLI Reference](CLI.md) for the complete list.

*Poké* can be used on multiple files:
```console
poke -dvp protein1.pdb protein2.pdb
```
In this case it will create two output folders called `protein1` and `protein2`.

To run the scouter on all frames of a trajectory, a reference structure must be provided to define the atom types (its positions will not be used):
```console
poke -dvp reference.pdb trajectory.xtc
```
Multiple trajectories can be concatenated:
```console
poke -dvp reference.pdb trajectory1.xtc trajectory2.xtc
```

If the input is a long list of different files they can be listed in a text file (`.txt`, one path per line) and passed to *Poké*:
```console
poke -dvp list.txt
```


## Structure Preparation

The `--prepare` flag enables automatic structure preparation before pocket scouting. This includes selecting the best alternate locations, adding missing hydrogen atoms according to the CHARMM force field, and assigning terminal caps to chains. When preparation is enabled, the extraction flag (`-e`) is automatically turned on so that the protonated structure is written to the output.

```console
poke -dvp --prepare protein.pdb
```

> [!TIP]
> When `--prepare` is used without `-ss`, the scouter selection is automatically set to the prepared atoms (proteins, nucleic acids, cofactors, and structural ions) instead of all atoms.

### Alt-location selection

When atoms have alternate locations (alt-locs), the preparator selects the best conformer for each residue based on highest occupancy. In case of a tie, the conformer with the lowest temperature factor is preferred.

### Hydrogen addition

Missing hydrogen atoms are added based on the CHARMM force field definitions. By default all existing hydrogens are removed and re-added. Use `--keep-hydrogens` to retain existing hydrogens and only add the missing ones:
```console
poke -dvp --prepare --keep-hydrogens protein.pdb
```

### Capping

Terminal caps are automatically added to chain ends. The preparator detects chain separations based on several criteria: changes in species type (protein, nucleic acid), segment ID (chain) changes, large jumps in residue IDs, and explicit termini indices.

The default caps are `NH3+`/`COO-` for proteins and `5TER`/`3TER` for nucleic acids. These can be customized:

| Flag | Description |
|---|---|
| `--caps-types TYPE N_TERM C_TERM ...` | Cap types for each species (e.g., `protein NH3+ COO-`) |
| `--caps-exclude NAME ...` | Residue names to skip when capping (default: `ACE NME`) |
| `--ignore-segid` | Do not use segment IDs to detect chain separations |
| `--max-other-len LEN` | Max "other" residues before a chain is considered interrupted (default: `5`) |
| `--resid-threshold LEN` | Residue ID jump to consider a chain interrupted (default: `50`, used only if segid is absent) |
| `--ter-indices IX ...` | Explicit 0-based residue indices to treat as terminal |

### Protonation states

Residues with titratable side chains can be assigned specific protonation states. The default is `auto`, which uses the presence of specific hydrogen atoms or geometric criteria (covalent/hydrogen bond detection) to determine the correct state.

| Flag | Options | Auto-detection method |
|---|---|---|
| `--his TYPE` | `auto`, `bond`, `HSP`, `HSE`, `HSD` | Presence of HD1/HE2 or hydrogen bonding |
| `--cys TYPE` | `auto`, `bond`, `CYS`, `CYS2` | Presence of HG1 or covalent bonding (disulfide) |
| `--lys TYPE` | `auto`, `LYS`, `LSN` | — |
| `--arg TYPE` | `auto`, `ARG`, `ARGN` | — |
| `--glu TYPE` | `auto`, `GLU`, `GLUP` | — |
| `--asp TYPE` | `auto`, `ASP`, `ASPP` | — |

For fine-grained control, `--custom-restype IX TYPE` can be used to override the residue type at specific 0-based indices:
```console
poke -dvp --prepare --custom-restype 42 HSE 105 CYS2 protein.pdb
```

### Preparation selection

By default, only proteins, nucleic acids, and common cofactors and structural ions are prepared. Use `--selection-prepare` with an [MDAnalysis selection string](https://docs.mdanalysis.org/stable/documentation_pages/selections.html) to control which atoms undergo preparation:
```console
poke -dvp --prepare --selection-prepare "protein" protein.pdb
```


## Output Descriptors

Use the `-d` flag to save pocket descriptors in a CSV file (`descriptors.csv`). The file contains one row per pocket with the following columns, in order:

| Column | Description |
|---|---|
| `OutputFrame` | Index of the frame in the output trajectory. Equal to the input frame index when no extraction or sorting is performed. |
| `OriginalFrame` | Index of the frame in the input trajectory. When multiple trajectories are concatenated this is the index in the concatenated trajectory (the reference frame is not counted). |
| `PocketResid` | Residue ID assigned to the pocket in the volume PDB file (starting at 800). |
| `Score` | Druggability score assigned by the neural network model. |
| `Volume` | Low-energy space available to the ligand, in Å³. Proportional to the number of grid points and the grid spacing. |
| `Surface` | Surface area of the pocket, in Å². A point is considered superficial when it has fewer than 5 adjacent points. The surface-to-volume ratio provides insight into pocket shape. |
| `Accessible` | Solvent-accessible surface area of the pocket, in Å². |
| `Length1` | First principal length of the pocket, computed as the square root of the largest eigenvalue of the covariance matrix of the pocket points. |
| `Length2` | Second principal length. |
| `Length3` | Third principal length (smallest). |
| `Hydrophobicity` | Fraction of pocket atoms with an absolute partial charge less than 0.3e. |
| `Enclosure` | Enclosure indicator, computed by comparing the number of neighboring atoms to the pocket surface. |
| `#PolarAtoms` | Number of pocket atoms with an absolute partial charge between 0.3e and 0.7e. |
| `#ChargedAtoms` | Number of pocket atoms with an absolute partial charge greater than 0.7e. |
| `Residues` | Space-separated list of residue IDs that form the pocket. |
| `Atoms` | *(Optional)* Space-separated list of atom indices belonging to the pocket. Only included when the `--write-atoms` flag is used. |


## Other Outputs

### Extracted trajectory

Normally *Poké* does not produce a new trajectory. Use the `-e` flag to extract frames that contain at least one pocket into a new trajectory file (`frames.xtc` + `structure.pdb`).

Frames can be sorted by any pocket property using the `-s` flag (descending order). Combined with `-n`, only the top *n* frames are kept:
```console
poke -dvp -e -s score -n 10 protein.pdb
```
The `OutputFrame` column in the descriptors CSV reflects the new frame ordering after sorting.

### Volume points

The `-v` flag writes a multi-model PDB file (`volumes.pdb`) containing the grid points of each pocket. Each pocket is assigned a unique residue ID (starting at 800) with residue name `PCK`. The B-factor column stores the electrostatic potential value at each point.

### Density maps

The `-dx` flag produces an OpenDX density map (`density.dx`) showing pocket occupancy across trajectory frames. The grid uses a 1.0 Å spacing with a 5.0 Å margin around the structure bounding box. If the `mrcfile` package is installed, an MRC format map (`density.mrc`) is also generated.

### PyMOL script

The `-p` flag generates a PyMOL script (`pymol.pml`) that loads the structure, trajectory, and volume files and sets up a ready-to-use visualization.

### Pickle

The `-k` flag saves the raw `Pocket` objects in a pickle file (`pockets.pickle`), one list per frame. This is useful for custom analysis in Python.

### Histograms

The `-g` flag generates histogram plots (`histograms/*.png`) showing the distribution of each descriptor across all pockets found. The appearance can be customized with environment variables:

| Variable | Description |
|---|---|
| `POKE_FONT` | Font family for plots |
| `POKE_FIGSIZE` | Figure size as comma-separated values (e.g., `8,4`) |
| `POKE_COLOR` | Override bar color |
| `POKE_QUANTILES` | Set to `no`, `false`, or `0` to hide quantile lines |

### Launch script

The `-ls` flag saves the full *Poké* command to a bash script (`launch.sh`) for reproducibility.

### Compression

The `-z` flag compresses the volume PDB and pickle outputs using gzip (`.pdb.gz`, `.pickle.gz`). Most downstream programs can read gzip-compressed files directly.


## Input Atoms Filtering
It is possible to use only a subset of the atoms to compute energy maps. The atoms can be selected using the `-ss` flag with an [MDAnalysis selection string](https://docs.mdanalysis.org/stable/documentation_pages/selections.html). For example, to select only protein atoms when the PDB file contains water molecules:
```console
poke -dvp -ss "protein" protein.pdb
```
The `-so` flag controls which atoms are written in the output trajectory. Atoms not in the `-so` selection are still used to compute energy maps unless removed with `-ss`.

## Pocket Filtering

Pockets can be filtered by any property or by the residues involved. The `-f` flag takes a Python expression; a pocket is kept only if the expression evaluates to `True`.

To select pockets with volume greater than 100 Å³:
```console
poke -dvp -f "volume > 100" protein.pdb
```

To combine conditions:
```console
poke -dvp -f "score > 6 and enclosure > 0.5" protein.pdb
```

### Available properties

All pocket descriptors are available in filter and sort expressions: `score`, `volume`, `surface`, `accessible`, `hydrophobicity`, `enclosure`, `lengths` (array of 3), `n_polar_atoms`, `n_charged_atoms`, `n_apolar_atoms`.

### Filter functions

| Function | Description |
|---|---|
| `residues(*resids, min_common=, max_extra=)` | Check if the pocket contains at least `min_common` residues from the given list and at most `max_extra` additional residues. |
| `jaccard(*resids)` | Jaccard similarity between the pocket's residues and the given residue IDs. |
| `tversky(*resids, a=, b=)` | Tversky similarity index with tunable weights. |
| `rmsd(selection="all")` | RMSD of the pocket atoms against the reference. Use `selection="CA"` for Cα-only RMSD. |
| `resids` | Set of residue IDs belonging to the pocket. |

Example — keep only pockets that overlap with residues 42, 58, and 102:
```console
poke -dvp -f "residues(42, 58, 102, min_common=2)" protein.pdb
```

### Sorting

The `-s` flag sorts pockets using the same expression syntax (descending order). After sorting, `-n` limits the output to the first *n* pockets:
```console
poke -dvp -s score -n 5 protein.pdb
```

## Python Package
The scouter can be used through Python for maximum flexibility:
```python
from poke import Poke
from MDAnalysis import Universe

u = Universe('1gkc.pdb') 

scouter = Poke()
scouter.from_mdanalysis(u.atoms)
pockets = scouter.get_pockets(u.atoms.positions)  

for p in pockets:
    print(p)
```

> [!NOTE] 
> When using the CLI, the `-k` flag can be used to save the `Pocket` objects in a pickle file for later use in Python.

## Plugins

### Chimera plugin

> [!WARNING]
> The Chimera plugin is under development and could be currently broken.

The plugin allows interactive pocket scouting directly inside UCSF Chimera. To install, copy the plugin directory to the Chimera share folder:
```console
cp -r plugins/chimera/Poke /path-to-chimera/share/
```