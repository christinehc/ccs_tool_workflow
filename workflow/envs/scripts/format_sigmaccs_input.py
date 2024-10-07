# imports
import pandas as pd
import utils

# define allowed atoms
SIGMACCS_ATOMS = ["C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I", "Co", "As", "Se"]
SIGMACCS_ADDUCTS = ["[M+H]+", "[M+Na]+", "[M-H]-"]

# format table to match sigmaccs reqs
data = pd.read_csv(snakemake.input[0], sep="\t")
data["True CCS"] = 100
data = data.rename(columns={"adduct": "Adduct", "smi": "SMILES"})

# filter out compounds with non-accepted atoms and adducts
data["bad_atoms"] = [
    utils.has_bad_atoms(s, i, SIGMACCS_ATOMS)
    for s, i in zip(data["SMILES"], data["inchi"])
]
data = data[~data["bad_atoms"]].reset_index(drop=True)
data = data[data["Adduct"].isin(SIGMACCS_ADDUCTS)]
print("before enforcing valid smi", len(data))
data = data.dropna(subset="SMILES")
print("after enforcing valid smi", len(data))

data.to_csv(snakemake.output[0], index=False)
