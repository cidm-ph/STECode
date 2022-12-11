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
        ("test1_sfindAbricate.tab", {"tox": ["tox1", "tox2"], "Stxtype": ["1A", "2C"]}),
        ("test2_sfindAbricate.tab", {"tox": ["tox1"], "Stxtype": ["1A"]}),
    ],
)
def test_stecvir_parser(input, expected):
    temp = pd.DataFrame(expected)
    shouldbe = temp.set_index("tox").T
    shouldbe.reset_index(drop=True, inplace=True)
    result = sp.stecfinder_input(os.path.join(data_dir, input))
    pdt.assert_frame_equal(result, shouldbe)


def test_sys_exit():
    input = os.path.join(data_dir, "null_sfindAbricate.tab")
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        sp.stecfinder_input(input)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1
