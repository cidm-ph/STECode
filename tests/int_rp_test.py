import pytest
import os
from vicost.parsers import recAstxeae_parse as rp
import pandas as pd
import pandas.testing as pdt

data_dir = "tests/data/"
heirarchy = ["12", "02", "01", "T2", "T1", "00"]


@pytest.mark.parametrize(
    "input, expected",
    [
        ("test1_2recAstxeae.txt", {"virgene": ["stx1a"], "iso_tox": ["00"]}),
        ("test2_2recAstxeae.txt", {"virgene": ["stx1a"], "iso_tox": ["00"]}),
    ],
)
def test_recA_parser(input, expected):
    shouldbe = pd.DataFrame(expected)
    shouldbe["iso_tox"] = pd.Categorical(
        shouldbe["iso_tox"], ordered=True, categories=heirarchy
    )
    result = rp.recA_input(os.path.join(data_dir, input))
    pdt.assert_frame_equal(result, shouldbe)


def test_sys_exit():
    input = os.path.join(data_dir, "null_2recAstxeae.txt")
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        rp.recA_input(input)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1
