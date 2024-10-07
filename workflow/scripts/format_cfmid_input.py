# imports
from os.path import basename, splitext
import pandas as pd

# unclear whether limitations exist -- below are deepccs
# CFMID_ATOMS = [
#     "C",
#     "H",
#     "N",
#     "O",
#     "P",
#     "S",
#     "F",
#     "I",
#     "B",
#     "Br",
#     "Cl",
#     "Co",
#     "As",
#     "Se",
# ]

# web server has many available adducts; CLI/docker only 2
# CFMID_ADDUCTS = [
#     "[M+H]+",
#     # "[M]+",
#     # "[M+NH4]+",
#     # "[M+Na]+",
#     # "[M+K]+",
#     # "[M+Li]+",
#     "[M-H]-",
#     # "[M]-",
#     # "[M+Cl]-",
#     # "[M+HCOOH-H]-",
#     # "[M+CH3COOH-H]-",
#     # "[M-2H]2-",
# ]

# run code
data = pd.read_csv(snakemake.input.data, sep="\t")
# data = data[data["adduct"] == snakemake.wildcards.adduct]

# write data to txt file
# format: a .txt file containing a list of space-separated (id, smiles_or_inchi)
for a in snakemake.params.allowed_adducts:
    for filename in snakemake.output:
        if splitext(basename(filename))[0] == a:
            with open(filename, "w") as f:
                tmp = data[data["adduct"] == a]
                for aid, smiles in zip(tmp["adduct_id"], tmp["smi"]):
                    f.write(f"{aid} {smiles}\n")
