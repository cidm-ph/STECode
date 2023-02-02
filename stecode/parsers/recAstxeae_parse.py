"""
Getting the necessary samtools data
"""
import logging
import sys
import os
import pandas as pd
import numpy as np

choices = ["02", "01", "T2", "T1", "00"]
heirarchy = ["12", "02", "01", "T2", "T1", "00"]

def recA_input(file):
    stx_df = pd.read_csv(file, sep="\t")
    if stx_df.empty:
        msg = "This sample did not contain any stx genes, please check again if this sample is STEC. Exiting"
        logging.error(msg)
        sys.exit(1)
    if stx_df["#rname"].str.contains("STECode_normalisation_recA").any() == False:
        msg = "The sample in question did not successfully detect the recA gene in samtools, retry or check manually"
        logging.error(msg)
        sys.exit(1)
    re_stx_df = stx_df[stx_df["coverage"] >= 98]
    if re_stx_df.empty:
        msg = "stx genes not detected in sequencing reads. Exiting"
        logging.error(msg)
        sys.exit(1)
    re_stx_df["virgene"] = (
        re_stx_df["#rname"].str.split("_").str.get(-1)
    )  # change this to last not second last.
    filt_stx_df = re_stx_df[re_stx_df["virgene"].str.contains("stx")]
    div = re_stx_df.loc[re_stx_df["#rname"] == "STECode_normalisation_recA"][
        "meandepth"
    ].values[0]
    filt_stx_df["Normalised"] = filt_stx_df["meandepth"] / div
    conditions = [
        (filt_stx_df["Normalised"] > 2.1)
        & (filt_stx_df["virgene"].str.contains("stx2")),  # Triggers 02
        (filt_stx_df["Normalised"] > 2.1)
        & (filt_stx_df["virgene"].str.contains("stx1")),  # Triggers 01
        (filt_stx_df["Normalised"] <= 2.1)
        & (filt_stx_df["Normalised"] >= 2)
        & (filt_stx_df["virgene"].str.contains("stx2")),  # Triggers T2
        (filt_stx_df["Normalised"] <= 2.1)
        & (filt_stx_df["Normalised"] >= 2)
        & (filt_stx_df["virgene"].str.contains("stx1")),  # Triggers T1
        (filt_stx_df["Normalised"] < 2),  # Triggers 00
    ]
    filt_stx_df["iso_tox"] = np.select(conditions, choices, default="XX")
    sub_stx_df = filt_stx_df.drop_duplicates(subset=["virgene"])
    keep_cols = ["virgene", "iso_tox"]
    get_12 = sub_stx_df["iso_tox"].value_counts()
    if get_12.get("01", 0) > 0 and get_12.get("02", 0) > 0:
        sub_stx_df.loc[:, "iso_tox"] = "12"
    sub_stx_df.loc[:, "iso_tox"] = np.where(
        (sub_stx_df["iso_tox"] == "01") & (sub_stx_df["iso_tox"] == "02"),
        "12",
        sub_stx_df["iso_tox"],
    )
    sub_stx_df = sub_stx_df.reindex(columns=keep_cols)
    sub_stx_df["iso_tox"] = pd.Categorical(
        sub_stx_df["iso_tox"], ordered=True, categories=heirarchy
    )
    sub_stx_df.sort_values(by="iso_tox", inplace=True)
    sub_stx_df.reset_index(drop=True, inplace=True)
    return sub_stx_df

def recA_cont(file):
    sub_stx_df = recA_input(file)
    result_df = sub_stx_df.truncate(after=0)
    return result_df


def run_skip(longread):
    if longread is True:
        stx_df = pd.DataFrame({"virgene": ["N/A"], "iso_tox": ["CG"]})
    else:
        stx_df = pd.DataFrame({"virgene": ["N/A"], "iso_tox": ["DG"]})
    return stx_df