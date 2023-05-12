![poke](./img/poke-banner.svg)

# Overview
Poke scouter is a blazing fast pocket scouter based on charmm forcefield.


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
All the input stuctures passed to *Poke* must be protonated. The name of the atoms should match the charmm or amber conventions.
*Poke* can be used from command line with:
```console
poke -dvp protein.pdb
```

The result will appear in a folder called the same as the input file (in this example `protein`) in the current working directory.
*Poke* can be used also on multiple files:
```console
poke -dvp protein1.pdb protein2.pdb
```
In this case it will create two output folders called `protein1` and `protein2`.

In order to run the scouter on all the frames of a trajectory a reference must be specified to know the atoms types (its positions will not be used):
```console
poke -dvp reference.pdb trajectory.xtc
```
If the input is a long list of different files they can be put in a text file `.txt` with a path for each line and passed to *Poke*:
```console
poke -dvp list.txt
```

## Descriptors
Poke produce a set of descriptors for each pocket useful to characterize it. The descriptors are used then by a neural network model to assign a score to the pocket.

``volume``
  The low-energy space available to the ligand measured in Å³. It is proportional to the number of points of the pocket and the grid spacing.

``surface``
  Proportional to the number of superficial points of the pockets measured in Å². A point is defined when the number of adjacent points is less than 5.
  The surface-to-volume ratio can provide insight into the shape of the pocket. 

``n_regions``
  Number of disconnected low-energy regions. A ligand can reach the regions with small linkers.

``lengths``
  Average lengths of the pocket measured along orthogonal axis. They are computed as the squared root of the eigenvalues of the covariance matrix of the points.

``hydrophobicity``
  Fraction of atoms belonging to the pocket with a absolute partial charge less than 0.3e .

``monopole``
  A measure of the net charge of the pocket.

``dipole``
  A measure of the eletrical dipole of the pocket.

``n_polar_atoms``
  Number of atoms belonging to the pocket with a absolute partial charge between 0.3e and 0.7e.

``n_charged_atoms``
  Number of atoms belonging to the pocket with a absolute partial charge greater than 0.7e.

``enclousure``
  An indicator of the enclousure of the pocket. Computed comparing the number of atoms neighboring the pocket and its surface.

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

## Pymol plugin

## Chimera plugin
For now the procedure to install the plugin is manual:
```console
cp -r plugins/chimera/Poke /path-to-chimera/share/
```