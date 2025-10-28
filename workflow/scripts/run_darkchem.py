# print version
import sys
from os.path import join

import numpy as np
import pandas as pd

import darkchem

print(sys.version)

# note that these paths need to be RELATIVE to where the notebook file is
#   MODEL_PATH = "/Users/imal967/pnnl/idpp/darkchem_models/ex_name/"  # downloaded from Deception, will always be 1 directory (ex: N7b_[M+H])
MODEL_PATH = {
    "[M+H]+": join(snakemake.config["darkchem_path"], "N7b_[M+H]"),
    "[M-H]-": join(snakemake.config["darkchem_path"], "N7b_[M-H]"),
}
# DATA_PATH = "./data/test_batch_size-10.tsv"  # some tsv file
DATA_PATH = snakemake.input[0]

# load data
data = pd.read_csv(DATA_PATH, sep="\t", header=0)
data = data[data["adduct"].isin(snakemake.config["adducts"])]  # restrict adducts

predicted_ccs = list()
for adduct in snakemake.config["adducts"]:
    tmp = data[data["adduct"] == adduct]

    # load models
    model = darkchem.utils.load_model(MODEL_PATH[adduct])
    x = np.array(
        [darkchem.utils.struct2vec(smi) for smi in tmp["smi"]]
    )  # extract only SMILES notation
    # essentially what this does is that in a for loop, it will go over each SMILES notation and turn it into a numerical encoding of it and create a table where each row is a SMILES notation

    # generate latent space
    x_latent = model.encoder.predict(x)

    # generate property predictions
    y_pred = model.predictor.predict(x_latent)

    predict_df = pd.DataFrame({"adductID": tmp["adduct_id"], "darkchem": y_pred[:, 1]})
    predicted_ccs.append(predict_df)

predicted_ccs = pd.concat(predicted_ccs, ignore_index=True)
predicted_ccs.to_csv(
    snakemake.output[0], sep="\t"
)  # this is likely the snakemake output
