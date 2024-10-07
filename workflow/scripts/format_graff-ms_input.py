# imports
import pandas as pd
import utils
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.rdmolops import GetFormalCharge

# allowed atoms
GRAFFMS_ATOMS = ["C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I"]
GRAFFMS_ADDUCTS = ["[M+H]+", "[M-H]-"]

# explode data to include desired CEs
data = pd.read_csv(snakemake.input.data, sep="\t")
energies = snakemake.config["collision_energy"]
data["collision_energy"] = [energies] * len(data)

# filter out compounds with non-accepted atoms and adducts
data["bad_atoms"] = [
    utils.has_bad_atoms(s, i, GRAFFMS_ATOMS) for s, i in zip(data["smi"], data["inchi"])
]
data = data[~data["bad_atoms"]]
data = data[data["adduct"].isin(GRAFFMS_ADDUCTS)]
data = data.dropna(subset="smi")

# enforce formal charge of < 2 -- problem mols appear to be >= +2
data["formal_charge"] = [GetFormalCharge(MolFromSmiles(s)) for s in data["smi"]]
data = data[data["formal_charge"] <= 1]

# expand all CEs
data = data.explode("collision_energy").reset_index(drop=True)

# should be four tab-delimited fields: a unique identifier for the spectrum, a SMILES string, a precursor type ([M+H]+ or [M-H]-) and a normalized collision energy
data = data[["adduct_id", "smi", "adduct", "collision_energy"]].rename(
    columns={"smi": "smiles"}
)
data.to_csv(snakemake.output.tsv, sep="\t", index=False)
