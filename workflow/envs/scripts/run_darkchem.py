import darkchem
import numpy as np
import pandas as pd
from os.path import join

# print version
import sys

print(sys.version)

# note that these paths need to be RELATIVE to where the notebook file is
#   MODEL_PATH = "/Users/imal967/pnnl/idpp/darkchem_models/ex_name/"  # downloaded from Deception, will always be 1 directory (ex: N7b_[M+H])
MODEL_PATH = join(
    snakemake.config["darkchem_path"], "N7b_[M+H]"
)  # probably something like this?
# DATA_PATH = "./data/test_batch_size-10.tsv"  # some tsv file
DATA_PATH = snakemake.input[0]  # something like this????

# load model
model = darkchem.utils.load_model(
    MODEL_PATH
)  # arguments.txt must be present in this folder, as well as respective network weights

# load data
# x = np.load(DATA_PATH)
data = pd.read_csv(DATA_PATH, sep="\t", header=0)  # read data in

x = np.array(
    [darkchem.utils.struct2vec(smi) for smi in data["smi"]]
)  # extract only SMILES notation
# essentially what this does is that in a for loop, it will go over each SMILES notation and turn it into a numerical encoding of it and create a table where each row is a SMILES notation

# generate latent space
x_latent = model.encoder.predict(x)

# generate property predictions
y_pred = model.predictor.predict(x_latent)

predicted_ccs = pd.DataFrame({"adductID": data["adduct_id"], "darkchem": y_pred[:, 1]})
# predicted_ccs.to_csv("./data/test_output_darkchem.tsv", sep="\t")
predicted_ccs.to_csv(
    snakemake.output[0], sep="\t"
)  # this is likely the snakemake output
