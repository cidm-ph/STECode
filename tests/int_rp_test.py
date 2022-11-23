import pytest
import os 
from vicost.parsers import recAstxeaehly_parse as rp
import pandas as pd
import pandas.testing as pdt

data_dir = "tests/data/"
heirarchy = ["12", "02", "01", "T2", "T1", "00", "CG"]

@pytest.mark.parametrize(
    "input, longread, expected",
    [
        ("test1_2recAstxeaehly.txt", False, {'Stxtype': ['Stx1'], 'NN2': ['00']}),
        ("test2_2recAstxeaehly.txt", False, {'Stxtype': ['Stx1'], 'NN2': ['00']}),
        ("test1_2recAstxeaehly.txt", True, {'Stxtype': ['Stx1'], 'NN2': ['CG']}),
        ("test2_2recAstxeaehly.txt", True, {'Stxtype': ['Stx1'], 'NN2': ['CG']})
    ]
)

def test_recA_parser(input, longread, expected):
    shouldbe = pd.DataFrame(expected)
    shouldbe["NN2"] = pd.Categorical(
        shouldbe["NN2"], ordered=True, categories=heirarchy
    )
    result = rp.recA_input(os.path.join(data_dir, input), longread)
    pdt.assert_frame_equal(result, shouldbe)

def test_sys_exit():
    input = os.path.join(data_dir, "null_2recAstxeaehly.txt")
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        rp.recA_input(input, False)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1