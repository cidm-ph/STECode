"""
Getting the necessary samtools data
"""

import pandas as pd
import numpy as np
import logging
import sys

heirarchy = ["12", "02", "01", "T2", "T1", "00"]

def recA_input(file):
    stx_df = pd.read_csv(file, sep="\t", header=0)
    if stx_df.empty:
        msg = "This sample did not contain any stx genes, please check again if this sample is STEC. Exiting"
        logging.error(msg)
        sys.exit(1)
    if stx_df["#rname"].str.contains("vicoSt_normalisation_recA").any() == False:
        msg = "The sample in question did not successfully detect the recA gene in samtools, retry or check manually"
        logging.error(msg)
        sys.exit(1)
    re_stx_df = stx_df[stx_df["coverage"] >= 95]
    re_stx_df["virgene"] = (
        re_stx_df["#rname"].str.split("_").str.get(-1)
    )  # change this to last not second last.
    filt_stx_df = re_stx_df[re_stx_df["virgene"].str.contains("stx")]
    div = re_stx_df.loc[re_stx_df["#rname"] == "vicoSt_normalisation_recA"][
        "meandepth"
    ].values[0]
    filt_stx_df["Normalised"] = filt_stx_df["meandepth"] / div
    conditions = [
        (filt_stx_df["Normalised"] > 2.2)
        & (filt_stx_df["virgene"].str.contains("stx1" and "stx2").any()),
        (filt_stx_df["Normalised"] > 2.2)
        & (filt_stx_df["virgene"].str.contains("stx2").any()),
        (filt_stx_df["Normalised"] > 2.2)
        & (filt_stx_df["virgene"].str.contains("stx1").any()),
        (filt_stx_df["Normalised"] <= 2.2)
        & (filt_stx_df["Normalised"] >= 2)
        & (filt_stx_df["virgene"].str.contains("stx2").any()),
        (filt_stx_df["Normalised"] <= 2.2)
        & (filt_stx_df["Normalised"] >= 2)
        & (filt_stx_df["virgene"] .str.contains("stx1").any()),
        (filt_stx_df["Normalised"] < 2),
    ]
    filt_stx_df["iso_tox"] = np.select(conditions, heirarchy, default="XX")
    sub_stx_df = filt_stx_df.drop_duplicates(subset=["virgene"])
    keep_cols = ["virgene", "iso_tox"]
    sub_stx_df = sub_stx_df.reindex(columns=keep_cols)
    sub_stx_df["iso_tox"] = pd.Categorical(
        sub_stx_df["iso_tox"], ordered=True, categories=heirarchy
    )
    sub_stx_df.sort_values(by="iso_tox", inplace=True)
    sub_stx_df.reset_index(drop = True, inplace=True)
    result_df = sub_stx_df.truncate(after=0)
    return result_df


def run_skip(longread):
    if longread is True:
        stx_df = pd.DataFrame({"virgene": ["N/A"], "iso_tox": ["CG"]})
    else:
        stx_df = pd.DataFrame({"virgene": ["N/A"], "iso_tox": ["DG"]})
    return stx_df
