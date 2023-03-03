"""
Getting the necessary vfdb data
"""

import pandas as pd
import logging
import sys

pd.set_option("display.max_rows", None)
pd.options.mode.chained_assignment = None


def stecfinder_input(test):
    vfdb_df = pd.read_csv(test, sep="\t", header=0)
    if vfdb_df.empty:
        msg = "This sample did not contain any stx genes, please check again if this sample is STEC or investigate the mapping results. Exiting"
        logging.error(msg)
        sys.exit(1)
    filt_vfdb_df = vfdb_df[vfdb_df["GENE"].str.contains("stx")]
    filt_vfdb_df = filt_vfdb_df[["GENE"]].drop_duplicates()
    return filt_vfdb_df


def stecfinder_cont(test):
    filt_vfdb_df = stecfinder_input(test)
    filt_vfdb_df["GENE"] = filt_vfdb_df["GENE"].str.upper()
    filt_vfdb_df["Stxtype"] = filt_vfdb_df["GENE"].str[-2:]
    temp_vfdb_df = (
        filt_vfdb_df.sort_values("Stxtype", ascending=True)
        .reset_index(drop=True)
        .drop(columns=["GENE"])
    )
    temp_vfdb_df.insert(0, "tox", range(1, 1 + len(temp_vfdb_df)))
    temp_vfdb_df["tox"] = "tox" + temp_vfdb_df["tox"].map(str)
    out_vfdb_df = temp_vfdb_df.set_index("tox").T
    return out_vfdb_df.reset_index(drop=True)
