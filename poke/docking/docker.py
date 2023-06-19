from .molecule import Molecule

class Docker:

    def dock_list(self, mols: list[Molecule]) -> list[Molecule]:
        raise NotImplementedError('Method dock not implemented.')
    
    def dock_dict(self, mols: dict[str, Molecule]) -> dict[str, list[Molecule]]:
        keys = list(mols.keys())
        docked = self.dock_list([mols[k] for k in keys])
        return {k: docked[i] for i, k in enumerate(keys)}

