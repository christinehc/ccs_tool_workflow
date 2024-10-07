# import python packages
import joblib
import sys
import pandas as pd
from rdkit.Chem import MolFromSmiles

sys.path.append(snakemake.config["rtp_path"])  # adding rtp to the path

# importing functions from rtp repo
from rtp.data import FeatureSet, get_features
from rtp.sklearn_models import model_predictions

# paths to data, pre-trained model, and params
model_path = snakemake.config["rtp_model_path"]
data_path = snakemake.input[0]
fp_size = snakemake.config["fp_size"]

# read in data
data = pd.read_csv(data_path, sep="\t", header=0)

# Filter out compounds that do not pass RDKit conversion
data["valid"] = [MolFromSmiles(s) is not None for s in data["smi"]]
data = data[data["valid"]].reset_index(drop=True)

pretrained_model = joblib.load(model_path)
smiles_list = list(data["smi"])

# get test features
test_features, smiles_none = get_features(smiles_list, FeatureSet.FP, fp_size)

# get model predictions
svm_pred = model_predictions(pretrained_model, test_features)

# format and output data
predicted_rtp = pd.DataFrame({"adductID": data["adduct_id"], "rtp": svm_pred})
predicted_rtp.to_csv(snakemake.output[0], sep="\t")
