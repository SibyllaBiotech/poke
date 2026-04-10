"""
Species determination module for molecular structures.

This module contains residue type definitions and classification functions
to identify different molecular species (proteins, nucleic acids, solvents, etc.)
based on residue names.
"""
from collections import Counter

# Common protein residues (standard amino acids and variants)
PROTEIN_RESIDUES = {
    # CHARMM top_all27_prot_lipid.rtf
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HSD',
    'HSE', 'HSP', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR',
    'TRP', 'TYR', 'VAL', 'ALAD',
    # PDB
    'HIS', 'MSE',
    # from Gromacs 4.5.3 oplsaa.ff/aminoacids.rtp
    'ARGN', 'ASPH', 'CYS2', 'CYSH', 'QLN', 'PGLU', 'GLUH', 'HIS1', 'HISD',
    'HISE', 'HISH', 'LYSH',
    # from Gromacs 4.5.3 gromos53a6.ff/aminoacids.rtp
    'ASN1', 'CYS1', 'HISA', 'HISB', 'HIS2',
    # from Gromacs 4.5.3 amber03.ff/aminoacids.rtp
    'HID', 'HIE', 'HIP', 'ORN', 'DAB', 'LYN', 'HYP', 'CYM', 'CYX', 'ASH',
    'GLH', 'ACE', 'NME',
    # from Gromacs 2016.3 amber99sb-star-ildn.ff/aminoacids.rtp
    'NALA', 'NGLY', 'NSER', 'NTHR', 'NLEU', 'NILE', 'NVAL', 'NASN', 'NGLN',
    'NARG', 'NHID', 'NHIE', 'NHIP', 'NTRP', 'NPHE', 'NTYR', 'NGLU', 'NASP',
    'NLYS', 'NPRO', 'NCYS', 'NCYX', 'NMET', 'CALA', 'CGLY', 'CSER', 'CTHR',
    'CLEU', 'CILE', 'CVAL', 'CASF', 'CASN', 'CGLN', 'CARG', 'CHID', 'CHIE',
    'CHIP', 'CTRP', 'CPHE', 'CTYR', 'CGLU', 'CASP', 'CLYS', 'CPRO', 'CCYS',
    'CCYX', 'CMET', 'CME', 'ASF',
    # from GROMACS 2019
    'ABU', 'AIB', 'ASP1', 'ASPP', 'CASPH', 'CCYS2', 'CDAB', 'CGLUH', 'CHISD', 
    'CHISE', 'CHISH', 'CHYP', 'CLYSN', 'CORN', 'CT3', 'CYN', 'DALA', 'GLUP', 
    'HISP', 'LSN', 'LYP', 'MELEU', 'MEVAL', 'NAC', 'NASPH', 'NCYS2', 'NDAB', 
    'NGLUH', 'NH2', 'NHE', 'NHISD', 'NHISE', 'NHISH', 'NHYP', 'NLE', 'NLYSN', 
    'NORN', 'PHEH', 'PHEU', 'PHL', 'TRPH', 'TRPU', 'TYRH', 'TYRU', 'CSO'
}

# Common nucleic acid residues (DNA/RNA bases and variants)
NUCLEIC_RESIDUES = {
    'A', 'T', 'G', 'C', 'U', 'DA', 'DT', 'DG', 'DC', 'DU',
    'ADE', 'THY', 'GUA', 'CYT', 'URA', 'ATP', 'ADP', 'GTP', 'GDP',
    'AMP', 'GMP', 'CMP', 'TMP', 'UMP', 'RA', 'RU', 'RG', 'RC',
    'DA5', 'DC5', 'DG5', 'DT5',
    'DA3', 'DC3', 'DG3', 'DT3',
    'RA5', 'RU5', 'RG5', 'RC5',
    'RA3', 'RU3', 'RG3', 'RC3'
}

# Solvent residues (water and organic solvents)
SOLVENT_RESIDUES = {
    # Water molecules
    'SOL', 'WAT', 'TIP3', 'TIP4', 'TIP5', 'SPC', 'SPCE', 'HOH', 'H2O', 'OHH',
    
    # Alcohols
    'MEOH', 'ETOH', 'PROP', 'BUTO', 'PENT', 'HEXA', 'HEPT', 'OCTA',
    'MET', 'ETH', 'PRP', 'BUT', 'PEN', 'HEX', 'HEP', 'OCT',
    'IPA', 'ISOP', 'TBA', 'TERT', 'IPA', 'IOH', 'EDO', 'EGL',
    'MPD', 'MRD', 'EOH',
    
    # Ethers
    'DMET', 'DIET', 'DIPR', 'DIBU', 'THF', 'DIOXN', 'DIOX', 'ANISO',
    'DME', 'DEE', 'DPE', 'DBE',
    
    # Ketones
    'ACET', 'MEK', 'MIBK', 'CYCLOHEX', 'BUTAN', 'PENTA', 'ACETO',
    'ACE', 'BUT2', 'PEN2', 'HEX2',
    
    # Esters
    'MEAC', 'ETAC', 'PRAC', 'BUAC', 'MEFO', 'ETFO', 'PRFO',
    'EAC', 'PAC', 'BAC', 'MFO', 'EFO', 'PFO',
    
    # Aromatic solvents
    'BENZ', 'TOLU', 'XYLO', 'STYR', 'ANIL', 'PHEN', 'NAPH',
    'BEN', 'TOL', 'XYL', 'STY', 'ANI', 'PHE', 'NAP',
    
    # Halogenated solvents
    'CHCL3', 'CCL4', 'DCM', 'DCE', 'TCE', 'PERC', 'FREON',
    'CHL', 'DICL', 'TRICL', 'TETR', 'METH', 'ETCL',
    
    # Nitriles
    'ACEN', 'PRPN', 'BUTN', 'BENZ', 'ACN', 'PRN', 'BUN',
    
    # Amides
    'DMF', 'DMA', 'NMP', 'FORM', 'ACET', 'PROP',
    'DMFO', 'DMAC', 'NMPR', 'FOR', 'ACM', 'PRM',
    
    # Sulfur-containing solvents
    'DMSO', 'SULFO', 'CS2', 'THIO', 'MSO', 'SUL', 'DMS', 'BME'
    
    # Nitrogen-containing solvents
    'PYRI', 'QUIN', 'PIPER', 'MORPH', 'IMID', 'TRIETHYLAMINE',
    'PYR', 'QUI', 'PIP', 'MOR', 'IMI', 'TEA',
    
    # Other common organic solvents
    'DMPU', 'HMPA', 'TMPO', 'GLYME', 'DIGLY', 'TEGDME',
    'DMP', 'HMP', 'TMP', 'GLY', 'DIG', 'TEG', 'SIN', 'GOL',
}

# Structural and catalytic ions (typically metal centers in proteins)
STRUCTURAL_ION_RESIDUES = {
    'CU1', 'CU1P'
    
    # Divalent structural ions
    'CA2+', 'MG2+', 'ZN2+', 'FE2+', 'CU2+', 'MN2+', 'NI2+', 'CO2+',
    'CA2P', 'MG2P', 'ZN2P', 'FE2P', 'CU2P', 'MN2P', 'NI2P', 'CO2P',
    'CA2', 'MG2', 'ZN2', 'FE2', 'CU2', 'MN2', 'NI2', 'CO2',
    'CA', 'MG', 'ZN', 'FE', 'CU', 'MN', 'NI', 'CO', 'CAL', 'MAG',
    
    # Trivalent structural ions
    'FE3+', 'AL3+', 'CR3+', 'AL', 'CR','GD3+','GD3P'
    'FE3P', 'AL3P', 'CR3P', 
    'FE3', 'AL3', 'CR3'
    
    # Other metal ions (catalytic/structural)
    'AG', 'AU', 'PB', 'CD', 'HG', 'TL', 'V', 'W', 'MO',
}

# Solvent ions (typically buffer/salt ions)
SOLVENT_ION_RESIDUES = {
    # Common monovalent ions (buffer/salt)
    'Na+', 'K+', 'Cl-', 'F-', 'BR-', 'I-', 'NA', 'K', 'CL', 'F', 'BR', 'I',
    'SOD', 'POT', 'CLA', 'FLU', 'BRO', 'IOD',
    
    # Polyatomic anions (buffer components)
    'SO4', 'SUL', 'CO3', 'PO4', 'VO4'
}

# Combined ion residues (for backward compatibility)
ION_RESIDUES = STRUCTURAL_ION_RESIDUES | SOLVENT_ION_RESIDUES


# Heme groups and cofactors
COFACTOR_RESIDUES = {
    'HEM', 'HEME', 'HEMA', 'HEMB', 'HEMC', 'FEOX', 'FRED',
    'NAD', 'NADH', 'NADP', 'NADPH', 'FAD', 'FADH', 'FMN', 'FMNH',
    'COA', 'COASH', 'TPP', 'PLP', 'B12', 'FOL', 'BIO',
    'ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'GMP', 'CTP', 'CDP', 'CMP',
    'TTP', 'TDP', 'TMP', 'UTP', 'UDP', 'UMP', 
    'ACT', 'CSD', 'MO3', 'MYR', 'AGC',
    'P22', 'YF3', 'YF4', 'P25',
}

# Lipids and membrane components
LIPID_RESIDUES = {
    'POPC', 'POPE', 'POPS', 'DOPC', 'DOPE', 'DPPC', 'DLPC', 'DMPC',
    'CHOL', 'CHOLESTEROL', 'PALM', 'PALMITATE', 'STEA', 'STEARATE',
    'OLEA', 'OLEATE', 'LINO', 'LINOLEATE'
}

# Sugars and carbohydrates
SUGAR_RESIDUES = {
    'GLC', 'GAL', 'MAN', 'FUC', 'XYL', 'RIB', 'ARA', 'GlcNAc', 'GalNAc',
    'NEU', 'NEUA', 'NEUR', 'GLCA', 'IDOA', 'GLCN', 'GALA', 'NAG', 'BMA', 
    'NDG', 'BOG'
}

SPECIES_CATEGORIES = {
    'protein': PROTEIN_RESIDUES,
    'nucleic_acid': NUCLEIC_RESIDUES,
    'solvent': SOLVENT_RESIDUES,
    'structural_ions': STRUCTURAL_ION_RESIDUES,
    'solvent_ions': SOLVENT_ION_RESIDUES,
    'ions': ION_RESIDUES,  # kept for backward compatibility
    'cofactor': COFACTOR_RESIDUES,
    'lipid': LIPID_RESIDUES,
    'sugar': SUGAR_RESIDUES
}

def classify(resnames: list[str], coverage_th: float=0.8) -> tuple[str, float]:
    """
    Classify the species type based on residue names in the component.
    
    Args:
        resnames: List of residue names in the component
        
    Returns:
        Tuple containing the species type classification and the coverage of that species type in the component.
    """
    
    unique_resnames = dict(Counter(resnames))

    for species_type, residues in SPECIES_CATEGORIES.items():
        cov = sum(unique_resnames[res] for res in (residues & unique_resnames.keys())) / len(resnames)
        if cov >= coverage_th:
            return species_type, cov
    
    
    return 'other', 0.0


def is_species(resname: str, species_type: str) -> bool:
    """
    Check if a residue name belongs to a specific species type.
    
    Args:
        resname: Residue name to check
        species_type: Species type to check against
        
    Returns:
        True if residue belongs to the species type
    """
    if species_type not in SPECIES_CATEGORIES:
        raise ValueError(f"Unknown species type: {species_type}")
    
    return resname in SPECIES_CATEGORIES[species_type]


def get_species(resname: str) -> str:
    """
    Get the species type for a single residue.
    
    Args:
        resname: Residue name
        
    Returns:
        Species type, or 'other' if not found
    """
    
    for species_type, residues in SPECIES_CATEGORIES.items():
        if resname in residues:
            return species_type
    
    return 'other'