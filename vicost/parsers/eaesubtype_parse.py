"""
Getting the necessary subtype from abricate data
"""

import pandas as pd
import logging


pd.set_option('display.max_rows', None)

def eaesubtype_input(file, name):
    eaesubtype_df = pd.read_csv(file, sep ='\t', header = 0)
    keep_cols = ['#FILE', 'GENE', 'NN1']
    eaesubtype_df=eaesubtype_df.reindex(columns=keep_cols)
    if eaesubtype_df.empty:
        msg = "This sample does not contain eae"
        logging.info(msg)
        eaesubtype_df = eaesubtype_df.append(pd.Series("00", index = eaesubtype_df.columns), ignore_index = True)
    eaesubtype_df['#FILE']= name
    eaesubtype_df['NN1'] = eaesubtype_df['GENE'].str[:2]
    return eaesubtype_df