# imports
import pickle
from os.path import join
import pandas as pd
from rdkit.Chem import MolFromSmiles

# enable c3sdb import
import sys

sys.path.append(snakemake.config["c3sdb_path"])
from c3sdb.ml.data import data_for_inference


# generate input data (m/z + encoded adduct + MQNs, centered and scaled) for inference
# mzs, adducts, and smis are all numpy arrays with same length
# included is a boolean mask with same shape as input arrays, indicating which rows
# are included in the generated dataset (some might fail to compute MQNs from SMILES)
DATA_PATH = DATA_PATH = snakemake.input[0]  # some tsv file
MODEL_PATH = join(snakemake.config["c3sdb_path"], "pretrained")

OHEncoder = join(MODEL_PATH, "c3sdb_OHEncoder.pkl")
SScaler = join(MODEL_PATH, "c3sdb_SScaler.pkl")
kmcm_svr = join(MODEL_PATH, "c3sdb_kmcm_svr.pkl")

data = pd.read_csv(DATA_PATH, sep="\t", header=0)  # read data in

# Filter out compounds that do not pass RDKit conversion
data["valid"] = [MolFromSmiles(s) is not None for s in data["smi"]]
data = data[data["valid"]].reset_index(drop=True)

# Create arrays
mzs = data["adduct_mz"]
adducts = data["adduct"]
smis = data["smi"]

X, included = data_for_inference(mzs, adducts, smis, OHEncoder, SScaler)

# load the trained model
with open(kmcm_svr, "rb") as pf:
    kmcm_svr = pickle.load(pf)

# do inference
y_pred = kmcm_svr.predict(X)

# format output files
predicted_ccs = pd.DataFrame({"adductID": data["adduct_id"], "c3sdb": y_pred})
predicted_ccs.to_csv(snakemake.output[0], sep="\t")
