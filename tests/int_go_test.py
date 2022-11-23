import pytest
import os 
from vicost import gen_output as go
import pandas as pd
import pandas.testing as pdt

data_dir = "tests/data/"

@pytest.mark.parametrize(
    "stfile, recAfile, virfile, longread, name, expected",
    [
        ("test1_eaesubtype.tab",
            "test1_stecvir.tab",
            "test1_2recAstxeaehly.txt",
            False,
            "test1",
            {'#FILE': ['test1'], 'NN1': ['06'], 'NN2': ['00'], 'NN3': ['1A'], 'NN4': ['2C'], 'NN5': ['00'], 'NN6': ['00'], 'Barcode': ['06-00-1A-2C-00-00']}),
        ("test2_eaesubtype.tab",
            "test2_stecvir.tab",
            "test2_2recAstxeaehly.txt",
            False,
            "test2",
            {'#FILE': ['test1'], 'GENE': ['06_Gamma1'], 'NN1': ['06']})
    ]
)

def test_gen_output(stfile, recAfile, virfile, longread, name, expected):
    shouldbe = pd.DataFrame(expected)
    print(shouldbe)
    result = go.merge_all_NNs(
        os.path.join(data_dir, stfile),
        os.path.join(data_dir, recAfile),
        os.path.join(data_dir, virfile),
        longread,
        name
    )
    print(result)
    pdt.assert_frame_equal(result, shouldbe)
