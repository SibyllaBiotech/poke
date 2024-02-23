![poke](./img/poke-banner.svg)

# Overview
Poke scouter is a blazing fast pocket scouter based on charmm forcefield. It uses an energy-based approach to identify the low-energy regions of the protein and it is able to compute a set of descriptors for each pocket useful to characterize it. The descriptors are used then by a neural network model to assign a score to the pocket. 

Designed to be used on long molecular dynamics trajectories, but it can be used also on a single frame.


## Installation
Poke is available as wheel package for macOS, Windows and Linux on PyPI. Install it using `pip`:
```console
python -m pip install poke
```
Alternatively, the package can be installed directly from `github`:
```console
git clone https://github.com/SibyllaBiotech/poke
pip install ./poke
```

## Getting started
> [!CAUTION]
> All the input stuctures passed to *Poke* must be protonated and The name of the atoms should match the charmm or amber conventions.

*Poke* can be used on a single structure with the following command:
```console
poke -dvp protein.pdb
```

The result will appear in a folder called the same as the input file (in this example `protein`) in the current working directory.
The flags `-dvp` are used to specify the desired output files:
- `-d` for the descriptors
- `-v` for the volumes (a pdb file with the points of the pocket)
- `-p` for a pymol script to visualize the pockets

use `-h` to see all the available options.

*Poke* can be used also on multiple files:
```console
poke -dvp protein1.pdb protein2.pdb
```
In this case it will create two output folders called `protein1` and `protein2`.

In order to run the scouter on all the frames of a trajectory a reference must be specified to know the atoms types (its positions will not be used):
```console
poke -dvp reference.pdb trajectory.xtc
```
Multiple trajectories can be concatenated with:
```console
poke -dvp reference.pdb trajectory1.xtc trajectory2.xtc
```

If the input is a long list of different files they can be put in a text file `.txt` with a path for each line and passed to *Poke*:
```console
poke -dvp list.txt
```

## Descriptors
Poke produce a csv file with all the pockets found. The file contains the following columns:

``OutputFrame``
  The index of the frame in the output trajectory. When no extraction or sorting is performed it is equal to the input frame.

``OriginalFrame``
  The index of the frame in the input trajectory.
  When multiple trajectories are concatenated it is the index in the concatenated trajectory. Notice that the reference frame is not included in the count.

``PocketResid``
  If there are multiple pockets in the same frame,
  the volume points are written with the resid specified in this column. Useless if sorting is performed.

``Score``
  The score assigned to the pocket by the neural network model. It is a measure of the pocket druggability. 

``Volume``
  The low-energy space available to the ligand measured in Å³. It is proportional to the number of points of the pocket and the grid spacing. 

``Surface``
  Proportional to the number of superficial points of the pockets measured in Å². A point is defined when the number of adjacent points is less than 5.
  The surface-to-volume ratio can provide insight into the shape of the pocket. 

``Lengths``
  Three numbers representing the average lengths of the pocket measured along orthogonal axis. They are computed as the squared root of the eigenvalues of the covariance matrix of the points.

``Hydrophobicity``
  Fraction of atoms belonging to the pocket with a absolute partial charge less than 0.3e .

``Enclosure``
  An indicator of the enclosure of the pocket. Computed comparing the number of atoms neighboring the pocket and its surface.

``#PolarAtoms``
  Number of atoms belonging to the pocket with a absolute partial charge between 0.3e and 0.7e.

``#ChargedAtoms``
  Number of atoms belonging to the pocket with a absolute partial charge greater than 0.7e.

## Python package
The scouter can be used through python to achieve maximum flexibility:
```python
from poke import Poke
from MDAnalysis import Universe

u = Universe('1gkc.pdb') 

scouter = Poke()
scouter.load_atoms_properties_from_mdanalysis(u.atoms)
pockets = scouter.get_pockets(u.atoms.positions)  

for p in pockets:
    print(p)
```

> [!NOTE] 
> When the scouter is used from command line, the option `-k` can be used to save the pockets in a pickle file.

## Chimera plugin
> [!CAUTION]
> The plugin is still under development.

The procedure to install the plugin is manual:
```console
cp -r plugins/chimera/Poke /path-to-chimera/share/
```