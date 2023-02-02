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
import sys

cols = ["eae_sub", "iso_tox", "tox1", "tox2", "tox3", "tox4"]

def pre_merge_check(recAfile, virfile):
    mid_recA = rp.recA_input(recAfile)
    mid_recA.drop(columns=["iso_tox"], inplace=True)
    mid_stx = vp.stecfinder_input(virfile)
    merge_df = mid_recA.join(mid_stx, how="outer")
    discrepant_values = merge_df.loc[merge_df['virgene'].ne(merge_df['GENE'])].fillna('')
    if discrepant_values.empty is False:
        dv_list = discrepant_values.T
        dv_string = dv_list.to_string(index=False, header=False)
        x = ', '.join(dv_string.split())
        msg = f"Discrepant stx genes ({x})were detected between mapping and abricate, please check raw files"
        logging.error(msg)
        sys.exit(1)

def merge_all_NNs(stfile, recAfile, virfile, reads, name, longread):
    NN1_df = ep.eaesubtype_input(stfile, name)
    if reads is True:
        NN2_df = rp.recA_cont(recAfile)
    else:
        NN2_df = rp.run_skip(longread)
    NN3456_df = vp.stecfinder_cont(virfile)
    merge_df = pd.concat([NN1_df, NN2_df, NN3456_df], axis=1, join="inner")
    if "tox1" not in merge_df:
        merge_df["tox1"] = "00"
    if "tox2" not in merge_df:
        merge_df["tox2"] = "00"
    if "tox3" not in merge_df:
        merge_df["tox3"] = "00"
    if "tox4" not in merge_df:
        merge_df["tox4"] = "00"
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
    if go_df.empty is True:
        with open(outfile, "w") as emptyfile:
            emptyfile.write(
                "Error: An error has occurred during the generation of your barcode, please check raw files, or retry."
            )
            emptyfile.close()
    go_df.to_csv(outfile, sep="\t", index=False)
    go_string = go_df["Virulence_barcode"].to_string(index=False, header=False)
    logging.info(f"Here is your barcode: {go_string}")