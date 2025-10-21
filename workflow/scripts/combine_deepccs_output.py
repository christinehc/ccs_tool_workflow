# imports
import re

import pandas as pd


# define function to detect adduct charge from adduct name
def get_adduct_charge(adduct):
    """Find adduct charge given an adduct name.

    Note: works on most adducts, may need tweaking to include more.

    Args:
    ----
    adduct : str
        Name of adduct as a string (e.g. 'M+H')

    Returns
    -------
    int
        Adduct charge (e.g. 1)

    """
    charge = 0
    neg_ions = ["Cl", "Br", "I", "OH", "HCO2", "H2"]
    h2_pattern = re.compile("-H-H(?![0-9|CHNOPSF])")
    h2_match = h2_pattern.search(str(adduct))
    if h2_match:
        adduct = adduct.replace("-H-H", "-H2")

    # only count -H toward negative charge
    for neg in adduct.split("-")[1:]:
        # if adduct is M-
        if not any(neg):
            charge -= 1
        else:
            if neg[0].isdigit() and neg[1:].split("+")[0] == "H":
                charge -= int(neg[0])
            elif neg.split("+")[0] == "H":
                charge -= 1
            elif neg.split("+")[0] in neg_ions:
                charge -= 1

    # count all +el toward positive charge
    for pos in adduct.split("+")[1:]:
        # if adduct is M+
        if not any(pos):
            charge += 1
        else:
            ions = pos[1:].split("-")[0]
            if pos[0].isdigit() and ions in neg_ions:
                charge -= int(pos[0])
            elif pos in neg_ions:
                charge -= 1
            elif pos[0].isdigit() and ions not in neg_ions:
                charge += int(pos[0])
            else:
                charge += 1

    # # account for coefficients
    # if any(adduct.split("M")[0]):
    #     charge = charge / int(adduct.split("M")[0])

    return charge


def format_adduct(adduct: str, z: int) -> str:
    """Apply [M+X]z format to adduct string.

    Parameters
    ----------
    adduct : str
        Adduct (e.g. "M-H")
    z : int
        Charge of adduct (e.g. -1)

    Returns
    -------
    str
        Formatted adduct (e.g. "[M-H]-")
    """
    # ensure not already bracketed
    if "[" in adduct and "]" in adduct:
        return adduct

    if z == 0:
        return f"[{adduct}]"

    # include charge information
    sign = "+" if z > 0 else "-"
    z = abs(z)

    if z == 1:
        return f"[{adduct}]{sign}"
    return f"[{adduct}]{z}{sign}"


# need ccs, adduct_id, src_id
data = list()
for f, d in zip(snakemake.input.files, snakemake.input.data):
    input_data = pd.read_csv(d, sep="\t")
    tmp = pd.read_csv(f).rename(
        columns={"Adducts": "adduct", "SMILES": "smi", "CCS_DeepCCS": "ccs"}
    )
    tmp["adduct_z"] = [get_adduct_charge(a) for a in tmp["adduct"]]
    tmp["adduct"] = [
        format_adduct(a, z) for a, z in zip(tmp["adduct"], tmp["adduct_z"])
    ]
    tmp = tmp.merge(
        input_data, how="left", on=["smi", "adduct_id", "adduct", "adduct_z"]
    )
    data.append(tmp)

data = pd.concat(data, ignore_index=True)
data["src"] = "deepccs"
data = data[["ccs", "adduct_id", "src"]]
data.to_csv(snakemake.output[0], sep="\t", index=False)
