"""
Getting the necessary samtools data
"""

import pandas as pd
import numpy as np
import logging

choices = ['12', '02', '01', 'T2', 'T1', '00']
heirachy = ['12', '02', '01', 'T2', 'T1', '00', 'CG']

def recA_input(file, longread):
    stx_df = pd.read_csv(file, sep ='\t', header = 0)
    re_stx_df = stx_df[stx_df['coverage'] >= 95]
    re_stx_df['virgene'] = re_stx_df['#rname'].str.split('~~~').str.get(-1)
    filt_stx_df = re_stx_df[re_stx_df['virgene'].str.contains('Stx')]
    filt_stx_df['Stxtype'], filt_stx_df['A or B'] = filt_stx_df['virgene'].str[0:4], filt_stx_df['virgene'].str[5]
    if longread == False:
        div = re_stx_df.loc[re_stx_df['#rname'] == 'recA_BA000007_3']['meandepth'].values[0]
        filt_stx_df['AvgReadDepth'] = filt_stx_df.groupby('Stxtype')['meandepth'].transform('mean')
        filt_stx_df['Normalised'] = filt_stx_df['AvgReadDepth']/div
        conditions = [
            (filt_stx_df['Normalised'] > 2.2) & (filt_stx_df['Stxtype'].str.contains("Stx1" and "Stx2").any()),
            (filt_stx_df['Normalised'] > 2.2) & (filt_stx_df['Stxtype'].str.contains("Stx2").any()),
            (filt_stx_df['Normalised'] > 2.2) & (filt_stx_df['Stxtype'].str.contains("Stx1").any()),
            (filt_stx_df['Normalised'] <= 2.2) & (filt_stx_df['Normalised'] >= 2) & (filt_stx_df['Stxtype'] == "Stx2"),
            (filt_stx_df['Normalised'] <= 2.2) & (filt_stx_df['Normalised'] >= 2) & (filt_stx_df['Stxtype'] == "Stx1"),
            (filt_stx_df['Normalised'] < 2)
        ]
        filt_stx_df['NN2'] = np.select(conditions, choices, default="-")
    else:
        filt_stx_df['NN2'] = "CG"
    sub_stx_df = filt_stx_df.drop_duplicates(subset = ['Stxtype'])
    keep_cols = ['Stxtype', 'NN2']
    sub_stx_df=sub_stx_df.reindex(columns=keep_cols)
    sub_stx_df['NN2'] = pd.Categorical(sub_stx_df['NN2'], ordered=True, categories=heirachy)
    sub_stx_df = sub_stx_df.sort_values(['Stxtype','NN2']).drop_duplicates('NN2')
    return sub_stx_df.reset_index(drop=True)
