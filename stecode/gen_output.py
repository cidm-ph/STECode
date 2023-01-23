"""
Subscript of STECode that converges the 3 input files from the parsers and generate the output file.
"""

import datetime
import logging
from .parsers import eaesubtype_parse as ep
from .parsers import recAstxeae_parse as rp
from .parsers import stecvir_parse as vp
import pandas as pd
import os


def merge_all_NNs(stfile, recAfile, virfile, reads, name, longread):
    NN1_df = ep.eaesubtype_input(stfile, name)
    if reads is True:
        NN2_df = rp.recA_input(recAfile)
    else:
        NN2_df = rp.run_skip(longread)
    NN3456_df = vp.stecfinder_input(virfile)
    merge_df = pd.concat([NN1_df, NN2_df, NN3456_df], axis=1, join="inner")
    if "tox1" not in merge_df:
        merge_df["tox1"] = "00"
    if "tox2" not in merge_df:
        merge_df["tox2"] = "00"
    if "tox3" not in merge_df:
        merge_df["tox3"] = "00"
    if "tox4" not in merge_df:
        merge_df["tox4"] = "00"
    cols = ["eae_sub", "iso_tox", "tox1", "tox2", "tox3", "tox4"]
    merge_df["Virulence_barcode"] = merge_df[cols].apply(
        lambda row: "-".join(row.values.astype(str)), axis=1
    )
    # Need to change all column names to below.
    keep_cols = [
        "#Sequence_ID",
        "eae_sub",
        "iso_tox",
        "tox1",
        "tox2",
        "tox3",
        "tox4",
        "Virulence_barcode",
    ]
    merge_df = merge_df.reindex(columns=keep_cols)
    return merge_df


def gen_output(name, output, go_df):
    now = datetime.datetime.now()
    date = now.strftime("%Y%m%d")
    outfile = os.path.join(output, name, name + "_virbarcode_" + date + ".tab")
    go_df.to_csv(outfile, sep="\t", index=False)
    logging.info("Here is your barcode:")
    print(go_df.to_string(index=False))
