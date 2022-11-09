"""
Getting the necessary samtools data
"""

import pandas as pd
import logging


pd.set_option('display.max_rows', None)

def samtools_input(file):
    samtools_df = pd.read_csv(file, sep ='\t', header = 0)
    keep_cols = ['#FILE', 'GENE', 'NN1']
    samtools_df['NN1'] = samtools_df['GENE'].str[:2]
    samtools_df=samtools_df.reindex(columns=keep_cols)
    return samtools_df