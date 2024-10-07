# imports
import pandas as pd

# annotate input idppdb format with ccs prediction output
data = list()
for f, d in zip(snakemake.input.files, snakemake.input.data):
    input_data = pd.read_csv(d, sep="\t")
    tmp = pd.read_csv(f, sep="\t")
    if "Unnamed: 0" in tmp.columns:
        tmp = tmp.drop(columns="Unnamed: 0")
    tmp = tmp.rename(columns={"adductID": "adduct_id","rtp":"retention_time"})
    tmp = tmp.merge(input_data, how="left", on=["adduct_id"])
    data.append(tmp)

data = pd.concat(data, ignore_index=True)
data["src"] = "rtp"
data = data[["retention_time", "adduct_id", "src"]]  # need ccs, adduct_id, src_id
data.to_csv(snakemake.output[0], sep="\t", index=False)