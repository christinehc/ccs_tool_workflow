# imports
from typing import Any, Callable
from rdkit import Chem


def mol_function(id_type: str) -> Callable[[str], Any]:
    """Parse identifier type into RDKit mol object generator.

    Parameters
    ----------
    id_type : str
        Type of molecular identifier ("smiles" or "inchi")

    Returns
    -------
    Callable[[str], Any]
        Function converting molecular identifier to function.

    Raises
    ------
    ValueError
        If `id_type` is not "smiles" or "inchi"
    """
    if id_type.lower() == "smiles":
        function = Chem.MolFromSmiles
    elif id_type.lower() == "inchi":
        function = Chem.MolFromInchi
    else:
        raise ValueError("`id_type` must be one of ('smiles', 'inchi').")
    return function


def get_atoms(identifier: str, id_type: str = "smiles") -> list[str]:
    """Get all atoms in a molecule represented by input SMILES.

    Parameters
    ----------
    identifier : str
        Value of molecular identifier (i.e. a SMILES or InChI)
    id_type : str
        Type of molecular identifier ("smiles" or "inchi")

    Returns
    -------
    list[str]
        List of all individual atoms contained in molecule
    """
    f = mol_function(id_type)
    mol = f(identifier)
    return [a.GetSymbol() for a in mol.GetAtoms()]


def has_bad_atoms(
    smiles: str,
    inchi: str,
    allowed_atoms: list[str] = ["C", "H", "N", "O", "P", "S"],
) -> bool:
    """Determine whether input SMILES contains disallowed atoms.

    Parameters
    ----------
    smiles : str
        SMILES string
    inchi : str
        InChI string
    allowed_atoms : list[str], optional
        List of allowed atoms, by default
            ["C", "H", "N", "O", "P", "S"]

    Returns
    -------
    bool
        Returns True if any atoms in SMILES are not in the allowed
            list.
    """
    try:
        atoms = get_atoms(smiles, id_type="smiles")
    except AttributeError:  # try inchi if smiles invalid
        try:
            atoms = get_atoms(inchi, id_type="inchi")
        except AttributeError:  # give up if nothing valid
            return True
    return any([a not in allowed_atoms for a in atoms])


def neutralize_atoms(mol):
    """_summary_

    From: https://www.rdkit.org/docs/Cookbook.html#neutralizing-molecules

    Parameters
    ----------
    mol : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol
