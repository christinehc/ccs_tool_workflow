# imports
from os.path import splitext
import pandas as pd

# concat data
# need: rt, adduct_id, src_id, 
data = list()
for f in snakemake.input:
    tmp = pd.read_excel(f)
    if "Unnamed: 0" in tmp.columns:
        tmp = tmp.drop(columns="Unnamed: 0")
    tmp.columns = [c.lower() for c in tmp.columns]
    tmp = tmp.rename(columns={"rtp": "rt", "name": "adduct_id"})
    tmp["adduct"] = splitext(f)[0].split("__")[1].split("_")[0]
    tmp["lc"] = splitext(f)[0].split("__")[1].split("_")[1].upper()
    data.append(tmp)

data = (
    pd.concat(data, ignore_index=True)
    .drop_duplicates()
    .sort_values(by=["adduct_id", "adduct", "lc"])
)
data = data[["rt", "adduct_id", "adduct", "lc", "formula", "smiles"]]
data["src"] = "retip"
data.to_csv(snakemake.output[0], index=False, sep="\t")
