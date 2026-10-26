# Changelog

All notable changes to Poké are recorded in this file, newest first. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the version numbers follow
[Semantic Versioning](https://semver.org/).

The user guide lives in a separate repository, [poke-docs](https://github.com/SibyllaBiotech/poke-docs).
Every entry here that changes what users see (a flag, an output file, a descriptor) must be mirrored there before the release.

## [0.2.0] - 2026-10-26

Reads mmCIF files. Requires Python 3.11 or later.

### Added
- mmCIF input files (`.cif`, `.mmcif`), the main format of the PDB and of the AlphaFold
  database. MDAnalysis reads them through the gemmi library, now a dependency of Poké.
  As with PDB files, the models of an NMR ensemble are read as frames.

### Changed
- Poké requires Python 3.11 or later and MDAnalysis 2.11 or later, the first version that
  reads mmCIF.

### Fixed
- `--prepare` crashed on single-frame files other than PDB, such as GRO and PDBQT
  ("GROReader only contains a single frame").

## [0.1.1] - 2026-10-20

Packaging and documentation for the first release on PyPI. The results are identical to 0.1.0.

### Changed
- The package is distributed on PyPI as `poke-scouter`; the command and the Python module
  are still called `poke`.
- The README is shipped inside the package metadata and serves as the PyPI project page.
- Command-line help: the options are grouped into sections, the values have named
  placeholders and the wording is corrected.

### Removed
- `CLI.md`: the command-line reference is now a page of the user guide, generated from the code.
- The developer-only code (force-field generator, docking benchmark, model training) is no
  longer in the repository. It was never part of the installed package.

### Fixed
- The Python example in the README did not run.
- Typos in the progress line and in an error message.

## [0.1.0] - 2026-08-20

The version used for the paper Colasurdo, Boldrini, Terruzzi et al., *ACS Omega* (2026), https://doi.org/10.1021/acsomega.6c04816, archived on Zenodo at https://doi.org/10.5281/zenodo.19634786.

Changes made during the peer review of the paper, compared with the internal version 0.0.1:

### Added
- `-mv`, `--min-volume` sets the minimum volume of the pocket core (default 15 Å³).
- `-ac`, `--accessible-cutoff` sets the distance around the pocket core within which low-energy
  points count as accessible volume (default 1.0 Å).
- `-ks`, `--keep-subpockets` keeps overlapping pockets instead of only the highest-scoring one.
- `--write-atoms` adds the pocket atom indices to the descriptors CSV file.
- Structure preparation handles structures with missing hydrogens.

### Changed
- The scoring model has a new architecture and was retrained in response to the peer review.
- Descriptors renamed: `Volume` is now `Core Volume`; `Accessible` is now `Accessible Volume`, in Å³.
- Overlapping pockets are pruned by default; only the highest-scoring one is kept.
- The pocket atom indices are no longer written to the descriptors CSV by default.
- Structure preparation is off by default; use `--prepare` to turn it on.

### Fixed
- Overflow when computing the pocket dimensions.
- Errors in the accessible volume computation.
- A bug in the structure preparator, and spurious MDAnalysis warnings during preparation.
