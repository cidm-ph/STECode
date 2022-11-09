"""
Getting the necessary vfdb data
"""

import pandas as pd
import logging

pd.set_option('display.max_rows', None)
pd.options.mode.chained_assignment = None

def vfdb_input(test):
    vfdb_df = pd.read_csv(test, sep ='\t', header = 0)
    filt_vfdb_df = vfdb_df[vfdb_df['RESISTANCE'].str.contains('Stx')]
    filt_vfdb_df['RESISTANCE'] = filt_vfdb_df['RESISTANCE'].str.upper()
    filt_vfdb_df['Stxtype'], filt_vfdb_df['A or B'] = filt_vfdb_df['RESISTANCE'].str[3:5], filt_vfdb_df['RESISTANCE'].str[5]
    filt_vfdb_df = filt_vfdb_df[['Stxtype']].drop_duplicates()
    temp_vfdb_df = filt_vfdb_df.sort_values('Stxtype', ascending=True).reset_index(drop=True)
    temp_vfdb_df.insert(0, 'NN', range(3, 3 + len(temp_vfdb_df)))
    temp_vfdb_df['NN'] = 'NN' + temp_vfdb_df['NN'].map(str)
    out_vfdb_df = temp_vfdb_df.set_index('NN').T
    return out_vfdb_df.reset_index(drop=True)