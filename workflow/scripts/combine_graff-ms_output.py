# imports
import pandas as pd

# annotate input idppdb format with ccs prediction output
data = list()
for f, d in zip(snakemake.input.files, snakemake.input.data):
    input_data = pd.read_csv(d, sep="\t")
    tmp = pd.read_csv(f, sep="\t").rename(
        columns={
            "compound": "adduct_id",
            "collision_energy": "ms2_ce",
            "mz": "ms2_mz",
            "intensity": "ms2_i",
        }
    )
    tmp = tmp.merge(input_data, how="left", on=["adduct_id", "adduct"])
    data.append(tmp)

# DB format: ms2_mz, ms2_i, adduct_id, src_id, ms2_ce
data = pd.concat(data, ignore_index=True)
data["src"] = "graff-ms"
data = data[
    ["ms2_mz", "ms2_i", "adduct_id", "src", "ms2_ce"]
]  # need ccs, adduct_id, src_id
data.to_csv(snakemake.output[0], sep="\t", index=False)
