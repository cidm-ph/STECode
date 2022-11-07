"""
Getting the necessary samtools data
"""

import pandas as pd
import numpy as np
import logging
import re

recAstxeaehly_test = "/Users/winx/Documents/Bioinformatics/stecode/tests/data/test1_2recAstxeaehly.txt"

def recA(file):
    stx_df = pd.read_csv(file, sep ='\t', header = 0)
#    if stx_df.empty:
#        msg = "the 2recAstxeaehly is empty"
#        logging.critical(msg,)
    re_stx_df = stx_df[stx_df['coverage'] >= 95]
    re_stx_df['rename'] = re_stx_df['#rname'].str.split('~~~', 2).str.get(-2)


    #else:
    #    print("This sample did not contain any Stx genes.")

    #keep_cols = ['#FILE', 'GENE', 'NN1']
    #samtools_df['NN1'] = samtools_df['GENE'].str[:2]
    #samtools_df=samtools_df.reindex(columns=keep_cols)
    return re_stx_df

print(recA(recAstxeaehly_test))