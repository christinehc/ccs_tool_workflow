import pandas as pd

RETIP_ADDUCTS = ["M+H]+", "[M-H]-"]

df = pd.read_csv(snakemake.input[0], sep="\t")[
    ["adduct_id", "aFormulas.form", "smi", "inchi_key"]
].rename(
    columns={
        "adduct_id": "name",
        "aFormulas.form": "formula",
        "smi": "smiles",
        "inchi_key": "inchikey",
    }
)
df.columns = [c.upper() for c in df.columns]
# df["INCHIKEY"] = [Chem.MolToInchiKey(Chem.MolFromSmiles(smi)) for smi in df["SMILES"]]
df.to_excel(snakemake.output[0], index=False)
