# imports
from os.path import basename, sep, splitext

import numpy as np
import pandas as pd


# functions
def cfmid2spectra(filename: str, output_format: str = "both") -> dict:
    """_summary_

    Parameters
    ----------
    filename : str
        _description_
    output : str, optional
        _description_, by default False

    Returns
    -------
    dict
        _description_
    """
    output = {"energy0": [], "energy1": [], "energy2": []}
    with open(filename, "r") as readf:
        get_spec = [False, False, False]
        for line in readf:
            if any(get_spec):
                try:
                    assert isinstance(
                        float(line.strip().split(" ")[0]), float
                    ), "Does not pass float conversion."
                    i = get_spec.index(True)
                    output[f"energy{i}"] += [
                        np.array(
                            [
                                float(line.strip().split(" ")[0]),
                                float(line.strip().split(" ")[1]),
                            ]
                        )
                    ]
                except ValueError:
                    pass
            if "energy" in line.strip():
                i = int(line.strip().replace("energy", ""))
                get_spec = [False, False, False]
                get_spec[i] = True
            if line.strip() == "":
                get_spec = [False, False, False]

    for energy, spec in output.items():
        if output_format == "mz":
            output[energy] = [f"{s[0]}" for s in spec]
            output[energy] = ",".join(output[energy])
        if output_format == "intensity":
            output[energy] = [f"{s[1]}" for s in spec]
            output[energy] = ",".join(output[energy])
        else:
            output[energy] = [f"{s[0]},{s[1]}" for s in spec]
            output[energy] = ";".join(output[energy])
    return output


# code
data = list()
# for filedir in snakemake.input.data:
#     output_files = glob.glob(join(glob.escape(filedir), "*.log"))
for of in snakemake.input.files:
    # get adduct id / adduct from filename
    adduct_id = splitext(basename(of.split(sep)[-1]))[0]
    adduct = of.split(sep)[-2]

    # get mz and intensity arrays
    mzs = cfmid2spectra(of, output_format="mz")
    mzs = pd.DataFrame({k: [v] for k, v in mzs.items()}).T.rename(columns={0: "ms2_mz"})
    iys = cfmid2spectra(of, output_format="intensity")
    iys = pd.DataFrame({k: [v] for k, v in iys.items()}).T.rename(columns={0: "ms2_i"})

    # reformat collision energy data
    mzs = (
        mzs.reset_index()
        .rename(columns={"index": "ms2_ce"})
        .replace({"energy0": 10.0, "energy1": 20.0, "energy2": 40.0})
    )
    iys = iys.reset_index(drop=True)

    tmp = mzs.join(iys, how="left").drop_duplicates()
    tmp["adduct_id"] = adduct_id
    tmp["adduct"] = adduct
    data.append(tmp)

data = pd.concat(data, ignore_index=True)
data["src"] = "cfm-id"
data = data[
    ["ms2_mz", "ms2_i", "adduct_id", "src", "ms2_ce"]
]  # need ccs, adduct_id, src_id
data.to_csv(snakemake.output.tsv, sep="\t", index=False)
