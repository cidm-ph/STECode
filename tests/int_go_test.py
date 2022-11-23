import pytest
import os 
from vicost import gen_output as go
import pandas as pd
import pandas.testing as pdt

data_dir = "tests/data/"
heirarchy = ["12", "02", "01", "T2", "T1", "00", "CG"]

@pytest.mark.parametrize(
    "stfile, recAfile, virfile, longread, name, expected",
    [
        ("test1_eaesubtype.tab",
            "test1_2recAstxeaehly.txt",
            "test1_stecvir.tab",
            False,
            "test1",
            {'#FILE': ['test1'], 'NN1': ['06'], 'NN2': ['00'], 'NN3': ['1A'], 'NN4': ['2C'], 'NN5': ['00'], 'NN6': ['00'], 'Barcode': ['06-00-1A-2C-00-00']}),
        ("test2_eaesubtype.tab",
            "test2_2recAstxeaehly.txt",
            "test2_stecvir.tab",
            False,
            "test2",
            {'#FILE': ['test2'], 'NN1': ['03'], 'NN2': ['00'], 'NN3': ['1A'], 'NN4': ['00'], 'NN5': ['00'], 'NN6': ['00'], 'Barcode': ['03-00-1A-00-00-00']})
    ]
)

def test_gen_output(stfile, recAfile, virfile, longread, name, expected):
    shouldbe = pd.DataFrame(expected)
    shouldbe["NN2"] = pd.Categorical(
        shouldbe["NN2"], ordered=True, categories=heirarchy
    )
#    print(shouldbe)
    result = go.merge_all_NNs(
        os.path.join(data_dir, stfile),
        os.path.join(data_dir, recAfile),
        os.path.join(data_dir, virfile),
        longread,
        name
    )
    print(result)
    pdt.assert_frame_equal(result, shouldbe)
