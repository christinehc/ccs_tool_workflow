# imports
import pandas as pd

# annotate input idppdb format with ccs prediction output
data = list()
for f, d in zip(snakemake.input.files, snakemake.input.data):
    input_data = pd.read_csv(d, sep="\t")
    tmp = (
        pd.read_csv(f)
        .rename(columns={"Adduct": "adduct", "SMILES": "smi", "Predicted CCS": "ccs"})
        .drop(columns=["True CCS"])
    )
    tmp = tmp.merge(input_data, how="left", on=["smi", "adduct"])
    data.append(tmp)

data = pd.concat(data, ignore_index=True)
data["src"] = "sigmaccs"
data = data[["ccs", "adduct_id", "src"]]  # need ccs, adduct_id, src_id
data.to_csv(snakemake.output[0], sep="\t", index=False)
