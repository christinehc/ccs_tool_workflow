# imports
import pandas as pd
import utils

# cannot handle anything with >2 characters minus Co, lol
DEEPCCS_ATOMS = [
    "C",
    "H",
    "N",
    "O",
    "P",
    "S",
    "F",
    "I",
    # "B",  # I'm getting KeyError with "B"
    "Br",
    "Cl",
    "Co",
    "As",
    "Se",
]
# Maximum Smile length that can be used with DeepCCS.
MAX_SMILES_LENGTH = 250.0

# Accepted adducts by DeepCCS
DEEPCCS_ADDUCTS = ["[M+H]+", "[M+Na]+", "[M-H]-", "[M-2H]2-"]

# run code
data = pd.read_csv(snakemake.input[0], sep="\t")
data = data.rename(columns={"smi": "SMILES", "adduct": "Adducts"})

# filter out compounds with non-accepted atoms and adducts
data["bad_atoms"] = [
    utils.has_bad_atoms(s, i, DEEPCCS_ATOMS)
    for s, i in zip(data["SMILES"], data["inchi"])
]
data = data[~data["bad_atoms"]]
data = data[data["Adducts"].isin(DEEPCCS_ADDUCTS)]
data = data[
    data["SMILES"].apply(lambda x: len(str(x)) > 1)
]  # remove compounds with a SMILES length of 1 (throws error)
data = data.dropna(subset="SMILES")

# reduce cols to retip input format
data = data[["SMILES", "Adducts", "adduct_id"]]
data["Adducts"] = [a.split("[")[1].split("]")[0] for a in data["Adducts"]]
data.to_csv(snakemake.output[0], index=False)
