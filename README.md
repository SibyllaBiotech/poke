![poke](./img/poke-banner.svg)

# Overview
Poke scouter is a blazing fast pocket scouter based on charmm forcefield. It uses an energy-based approach to identify the low-energy regions of the protein and it is able to compute a set of descriptors for each pocket useful to characterize it. The descriptors are used then by a neural network model to assign a druggability score to the pocket. 

Designed to be used on long molecular dynamics trajectories, but it can be used also on single frames.


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
> All input stuctures passed to *Poke* must be protonated and the atoms names should match the charmm or amber convention.

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
  The index of the frame in the output trajectory. When no extraction or sorting is performed it is equal to the input frame. See the section on the output trajectory for more details.

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


## Output trajectory
Normally poke does not produce a new trajectory, but it is possible to extract the frames with at least a pocket using the `-e` flag. In this case, the extracted frames are written in a new trajectory file. The frames can be sorted by a property using the `-s` flag. The sorting is performed in descending order. 
For example to extract the frames with at least a pocket and sort them by the score of the pocket use:

```console
poke -dvp -e -s score protein.pdb
```

## Atoms filtering
It is possible to use only a subset of the atoms in the input file to compute energy maps. The atoms can be selected using the `-ss` flag and a selection string. The selection string must be in the format of the `MDAnalysis` package. For example if the pdb file contains water molecule,  to select only the protein atoms use:
```console
poke -dvp -ss "protein" protein.pdb
```
Using the `-so` flag it is possible to select the atoms of the input structure that should be written in the output trajectory. Atoms not present in the `-so` selection are anyway used to compute the energy maps unless the `-ss` flag is used to remove them.

## Pocket filtering 
The pockets can be filtered by properties or by the residues involved. The filtering is performed using the `-f` flag and is a python expression. A pocket is kept only if the expression evaluates to `True`.

For example, to select only the pockets with a volume greater than 100 Å³:
```console
poke -dvp -f "volume > 100" protein.pdb
```

To select only the pockets with a score greater than 6 and with enclosure greater than 0.5 use:
```console
poke -dvp -f "score > 6 and enclosure > 0.5" protein.pdb
```

After the sorting, the `-n` flag can be used to select only the first `n` pockets. For example to select only the best 5 pockets by score use:
```console
poke -dvp -s score -n 5 protein.pdb
```

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