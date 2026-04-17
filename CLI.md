# CLI Reference

Complete list of command-line options for *Pokè*.

```
poke [options] structures [structures ...]
```

## Positional Arguments

| Argument | Description |
|---|---|
| `structures` | Input structures (`.pdb`, `.gro`, `.pdbqt`, `.xtc`, `.trr`, `.dcd`, ...). When a coordinate-only file is passed, the preceding structure is used to open it. A `.txt` file can be passed containing additional structure paths, one per line. |

## General Options

| Flag | Default | Description |
|---|---|---|
| `-o`, `--output PATH` | Current directory | Output path where result folders are created. |
| `-q`, `--quiet` | `False` | Suppress progress output for long trajectories. |
| `-nc`, `--ncores N` | Auto | Number of cores for parallel computation of trajectories or multi-model PDB files. |
| `-z`, `--compress` | `False` | Compress volume PDB and pickle outputs using gzip. |
| `-h`, `--help` | — | Show help message and exit. |

## Output Options

| Flag | Default | Description |
|---|---|---|
| `-d`, `--descriptors` | `False` | Compute descriptors and save them in a CSV file. |
| `--write-atoms` | `False` | Add the pocket atom indices to the descriptors CSV (requires `-d`). |
| `-v`, `--volumes` | `False` | Save pocket volumes as a multi-model PDB file. |
| `-dx`, `--density` | `False` | Save a density map of pocket occupancy in DX format. |
| `-e`, `--extract` | `False` | Extract frames with at least one pocket into a new trajectory. |
| `-p`, `--pymol` | `False` | Generate a PyMOL script to visualize the results. |
| `-g`, `--graph` | `False` | Plot descriptor distribution histograms. |
| `-k`, `--pickle` | `False` | Save pocket objects in a pickle file for Python analysis. |
| `-ls`, `--launch-script` | `False` | Save the launch command to a bash script for reproducibility. |

## Filtering and Sorting

| Flag | Default | Description |
|---|---|---|
| `-f`, `--filter EXPR` | None | Filter pockets with a Python expression. A pocket is kept only if the expression evaluates to `True`. All pocket properties and filter functions are available (see [Pocket Filtering](README.md#pocket-filtering)). |
| `-s`, `--sort EXPR` | None | Sort pockets by the value of a Python expression in descending order. Uses the same syntax as `--filter`. |
| `-n`, `--max-output-frames N` | None | Maximum number of pockets in the output. Only effective when `--sort` is used. |
| `-ks`, `--keep-subpockets` | `False` | Keep all pockets even if they share regions, instead of retaining only the highest-scoring one. |

## Alignment

| Flag | Default | Description |
|---|---|---|
| `-a`, `--align SEL` | None | Align the given atom selection to the reference before scouting. Uses MDAnalysis selection syntax. |
| `-r`, `--reference PATH` | None | Reference structure for alignment. Defaults to the first input structure. |

## Atom Selection

| Flag | Default | Description |
|---|---|---|
| `-ss`, `--selection-scouter SEL` | `all` (or auto with `--prepare`) | Selection of atoms used by the scouter. Uses MDAnalysis selection syntax. When `--prepare` is used without `-ss`, defaults to the prepared atoms selection. |
| `-so`, `--selection-output SEL` | `all` | Selection of atoms written in the output trajectory. Atoms outside this selection are still used for energy maps unless removed by `-ss`. |

## Scoring

| Flag | Default | Description |
|---|---|---|
| `-sf`, `--scoring-function NAME` | `model` | Scoring function to use for pocket druggability assessment. |

## Structure Preparation

These options control the automatic preparation of input structures. See [Structure Preparation](README.md#structure-preparation) for details.

| Flag | Default | Description |
|---|---|---|
| `--prepare` | `False` | Prepare structures before scouting (select alt-locs, add hydrogens, and caps). Default is to NOT prepare. |
| `--keep-hydrogens` | `False` | Keep existing hydrogens during preparation, adding only the missing ones. |
| `--selection-prepare SEL` | Auto | Selection of atoms to prepare. Others are kept as-is. If not provided, only proteins, nucleic acids, and common cofactors and structural ions are prepared. |

### Capping

| Flag | Default | Description |
|---|---|---|
| `--caps-types TYPE N_TERM C_TERM ...` | `protein NH3+ COO- nucleic_acid 5TER 3TER` | Types of caps to add, specified as triplets of species type, N-terminal cap, and C-terminal cap. |
| `--caps-exclude NAME ...` | `ACE NME` | Residue names to exclude from capping. |
| `--ignore-segid` | `False` | Do not use segment ID (chain) values to detect chain separations. |
| `--max-other-len LEN` | `5` | Maximum number of "other" residues allowed before a chain is considered interrupted. |
| `--resid-threshold LEN` | `50` | Residue ID jump required to consider a chain interrupted (used only if segment ID is absent). |
| `--ter-indices IX ...` | None | Explicit 0-based indices of residues to treat as terminal for chain splitting. |

### Protonation States

| Flag | Default | Options | Description |
|---|---|---|---|
| `--his TYPE` | `auto` | `auto`, `bond`, `HSP`, `HSE`, `HSD` | Histidine protonation state. `auto` uses HD1/HE2 presence or hydrogen bonding detection. |
| `--cys TYPE` | `auto` | `auto`, `bond`, `CYS`, `CYS2` | Cysteine protonation state. `auto` uses HG1 presence or covalent bonding (disulfide) detection. |
| `--lys TYPE` | `auto` | `auto`, `LYS`, `LSN` | Lysine protonation state. |
| `--arg TYPE` | `auto` | `auto`, `ARG`, `ARGN` | Arginine protonation state. |
| `--glu TYPE` | `auto` | `auto`, `GLU`, `GLUP` | Glutamate protonation state. |
| `--asp TYPE` | `auto` | `auto`, `ASP`, `ASPP` | Aspartate protonation state. |
| `--custom-restype IX TYPE ...` | None | — | Override residue type at specific 0-based indices (e.g., `--custom-restype 42 HSE 105 CYS2`). |
