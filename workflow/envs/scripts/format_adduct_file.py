import pandas as pd

data = pd.read_csv(snakemake.input[0])

if "exposome" in snakemake.config["match_file"]:
    data = data[data["UniqueID4DfileNames"] != "AgTune"].reset_index(drop=True)

# assign unique numerical ID per compound type
compounds = {
    cpd: i + 1
    for i, cpd in enumerate(
        data[snakemake.config["adducts"]["msac"]["id_col"]].unique()
    )
}
data["compound_id"] = [
    compounds[cpd] for cpd in data[snakemake.config["adducts"]["msac"]["id_col"]]
]

# get approximate chemical formula from smiles
data["formula"] = [
    utils.smiles_to_atoms(smi)
    for smi in data[snakemake.config["adducts"]["msac"]["smiles_col"]]
]

# separate adduct columns from non-adduct columns
# adduct_cols = [
#     col for col in data.columns if col.split("_")[-1][:3] == "ESI"
#     ]
# id_cols = [col for col in data.columns if col not in adduct_cols]
# data = pd.melt(data, id_vars=id_cols, value_vars=adduct_cols,
#                var_name='adduct', value_name='adduct_mass')

# clean adduct identifiers and make separate charge column
ref = pd.read_csv(params.ADDUCT_FILE)
ref = {add: z for add, z in zip(ref["adduct"], ref["charge"])}

# data['adduct'] = [adduct.split("_")[0] for adduct in data['adduct']]
data["adduct_mass"] = data["adduct mass"]
data["adduct_charge"] = [ref[adduct] for adduct in data["adduct"]]
data["esi_mode"] = [
    "neg" if int(z) < 0 else "pos"
    for adduct, z in zip(data["adduct"], data["adduct_charge"])
]

# check whether adducts can form based on formula
data["adduct_in_parent"] = [
    utils.adduct_in_parent(add, fmla, ignore_h=True)
    for add, fmla in zip(data["adduct"], data["formula"])
]
print(len(data))
data = data[data["adduct_in_parent"]].reset_index(drop=True)
print("after adduct_in_parent", len(data))

# remove checker column and drop duplicate or invalid rows
data = data.drop(columns=["formula", "adduct_in_parent", "adduct mass"])
data = data.drop_duplicates().reset_index(drop=True)
print("after drop_duplicates", len(data))
data = data[data["adduct_mass"].notna()].reset_index(drop=True)
print("after adduct_mass notna()", len(data))

data.to_csv(snakemake.output[0], index=False)
