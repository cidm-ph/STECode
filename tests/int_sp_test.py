import pytest
import os 
from vicost.parsers import stecvir_parse as sp
import pandas as pd
import numpy as np
import pandas.testing as pdt


data_dir = "tests/data/"

@pytest.mark.parametrize(
    "input, expected",
    [
        ("test1_stecvir.tab", {'NN':['NN3', 'NN4'], 'Stxtype': ['1A', '2C']}),
        ("test2_stecvir.tab", {'NN':['NN3'], 'Stxtype': ['1A']}),
    ]
)

def test_stecvir_parser(input, expected):
    temp = pd.DataFrame(expected)
    shouldbe = temp.set_index("NN").T
    shouldbe.reset_index(drop=True, inplace=True)
    result = sp.vfdb_input(os.path.join(data_dir, input))
    pdt.assert_frame_equal(result, shouldbe)