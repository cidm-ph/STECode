"""
Getting the necessary vfdb data
"""

import pandas as pd
import logging

stecvir_test = "/Users/winx/Documents/Bioinformatics/stecode/tests/data/test1_stecvir.tab"

pd.set_option('display.max_rows', None)

def vfdb_input(test):
    vfdb_df = pd.read_csv(test, sep ='\t', header = 0)
    #keep_cols = ['#FILE', 'GENE', 'NN1']
    #samtools_df['NN1'] = samtools_df['GENE'].str[:2]
    #samtools_df=samtools_df.reindex(columns=keep_cols)
    return vfdb_df