# =========================================================
# Imports
# =========================================================
import sys
from os import makedirs
from os.path import exists, join

sys.path.append("../repos/idpp/")
import pandas as pd

# =========================================================
# Functions
# =========================================================


# adapted from https://stackoverflow.com/questions/44729727/pandas-slice-large-dataframe-into-chunks/44729807#comment126483047_44729807
def chunk_df_by_compound_id(df: pd.DataFrame, chunk_size: int):
    cids = df["Adducts.cmpd_id"].unique()
    start = 0

    # if df is smaller than chunk, yield df
    if len(cids) <= chunk_size:
        yield df[:]
        return

    # yield valid chunks
    while start + chunk_size <= len(cids):
        subset = cids[start : start + chunk_size]
        yield df[df["Adducts.cmpd_id"].isin(subset)].reset_index(drop=True)
        start += chunk_size

    # yield final chunk
    if (start + chunk_size > len(cids)) & (start != len(cids)):
        subset = cids[start:]
        yield df[df["Adducts.cmpd_id"].isin(subset)].reset_index(drop=True)


# =========================================================
# Load DB, pull data, and format
# =========================================================
# if splitext(snakemake.input[0])[1] == "db":
#     db = IdPPdb(snakemake.input[0])

#     # Capture adduct table
#     adduct_data = {
#         "adduct_id": [],
#         "adduct": [],
#         "adduct_z": [],
#         "adduct_mz": [],
#         "Adducts.cmpd_id": [],
#         "cmpd_name": [],
#         "Adducts.form_id": [],
#         "aFormulas.form": [],
#     }
#     for i, adduct, z, mz, cid, cname, fid, formula in db.fetch_adduct_data():
#         adduct_data["adduct_id"] += [i]
#         adduct_data["adduct"] += [adduct]
#         adduct_data["adduct_z"] += [z]
#         adduct_data["adduct_mz"] += [mz]
#         adduct_data["Adducts.cmpd_id"] += [cid]
#         adduct_data["cmpd_name"] += [cname]
#         adduct_data["Adducts.form_id"] += [fid]
#         adduct_data["aFormulas.form"] += [formula]

#     adduct_data = pd.DataFrame(adduct_data)

#     # Capture compound table
#     compound_data = {
#         "Adducts.cmpd_id": [],
#         "cmpd_name": [],
#         "Adducts.form_id": [],
#         "formula": [],
#         "smi_id": [],
#         "smi": [],
#         "inchi_id": [],
#         "inchi": [],
#         "inchikey": [],
#     }
#     for (
#         cid,
#         cname,
#         fid,
#         formula,
#         sid,
#         smiles,
#         iid,
#         inchi,
#         inchikey,
#     ) in db.fetch_cmpd_data():
#         compound_data["Adducts.cmpd_id"] += [cid]
#         compound_data["cmpd_name"] += [cname]
#         compound_data["Adducts.form_id"] += [fid]
#         compound_data["formula"] += [formula]
#         compound_data["smi_id"] += [sid]
#         compound_data["smi"] += [smiles]
#         compound_data["inchi_id"] += [iid]
#         compound_data["inchi"] += [inchi]
#         compound_data["inchikey"] += [inchikey]

#     compound_data = pd.DataFrame(compound_data)

#     # Merge tables
#     data = adduct_data.merge(
#         compound_data, how="left", on="compound_id", suffixes=(None, "_compound")
#     )

# else:
data = pd.read_csv(snakemake.input[0], sep="\t")

# Remove invalid/unwanted entries for now
data = data[
    (data["smi"].notna()) & (data["adduct"].isin(snakemake.config["adducts"]))
].reset_index(drop=True)

# chunk and save
digits = len(str(len(data)))  # get num digits

if not exists(snakemake.output[0]):
    makedirs(snakemake.output[0])
for i, chunk in enumerate(
    chunk_df_by_compound_id(data, snakemake.config["db"]["chunk_size"])
):
    chunk.to_csv(
        join(snakemake.output[0], f"{i:0{digits}d}.tsv"), index=False, sep="\t"
    )
