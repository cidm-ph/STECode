"""
Subscript of stecode that converges the 3 input files from the parsers and generate the output file.
"""

import datetime
import parsers.samtools_parse as sp
import parsers.recAstxeaehly_parse as rp
import parsers.stecvir_parse as vp
import pandas as pd
import os

def merge_all_NNs(stfile, recAfile, virfile, longread):
    NN1_df = sp.samtools_input(stfile)
    NN2_df = rp.recA_input(recAfile, longread)
    NN3456_df = vp.vfdb_input(virfile)
    merge_df = pd.concat([NN1_df, NN2_df, NN3456_df], axis=1, join="inner")
    if 'NN3' not in merge_df:
        merge_df['NN3'] = "00"
    if 'NN4' not in merge_df:
        merge_df['NN4'] = "00"
    if 'NN5' not in merge_df:
        merge_df['NN5'] = "00"
    if 'NN6' not in merge_df:
        merge_df['NN6'] = "00"
    cols = ['NN1', 'NN2', 'NN3', 'NN4', 'NN5', 'NN6']
    merge_df['Barcode'] = merge_df[cols].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)
    keep_cols = ['#FILE', 'NN1', 'NN2', 'NN3', 'NN4', 'NN5', 'NN6', 'Barcode']
    merge_df=merge_df.reindex(columns=keep_cols)
    return merge_df

def gen_output(stfile, recAfile, virfile, longread, name, output):
    now = datetime.datetime.now()
    date = now.strftime("%Y%m%d")
    outfile = os.path.join(output, name + "_virbarcode_" + date + ".tab")
    output_df = merge_all_NNs(stfile, recAfile, virfile, longread)
    output_df.to_csv(outfile, sep="\t", index=False)