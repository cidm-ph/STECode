import pytest
import os 
from vicost.parsers import eaesubtype_parse as ep
import pandas as pd
import pandas.testing as pdt

data_dir = "tests/data/"

@pytest.mark.parametrize(
    "input, name, expected",
    [
        ("test1_eaesubtype.tab", "test1", {'#FILE': ['test1'], 'GENE': ['06_Gamma1'], 'NN1': ['06']}),
        ("test2_eaesubtype.tab", "test2", {'#FILE': ['test2'], 'GENE': ['03_Beta1'], 'NN1': ['03']}),
        ("null_eaesubtype.tab", "null", {'#FILE': ['null'], 'GENE': ['00'], 'NN1': ['00']}),
    ]
)

def test_eaesubtype_parser(input, name, expected):
    shouldbe = pd.DataFrame(expected)
    result = ep.eaesubtype_input(os.path.join(data_dir, input), name)
    pdt.assert_frame_equal(result, shouldbe)
