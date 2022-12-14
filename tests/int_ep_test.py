import pytest
import os
from stecode.parsers import eaesubtype_parse as ep
import pandas as pd
import pandas.testing as pdt

data_dir = "tests/data/"


@pytest.mark.parametrize(
    "input, name, expected",
    [
        (
            "test1_eaesubtype.tab",
            "test1",
            {"#Sequence_ID": ["test1"], "GENE": ["06_Gamma1"], "eae_sub": ["06"]},
        ),
        (
            "test2_eaesubtype.tab",
            "test2",
            {"#Sequence_ID": ["test2"], "GENE": ["03_Beta1"], "eae_sub": ["03"]},
        ),
        (
            "null_eaesubtype.tab",
            "null",
            {"#Sequence_ID": ["null"], "GENE": ["00"], "eae_sub": ["00"]},
        ),
    ],
)
def test_eaesubtype_parser(input, name, expected):
    shouldbe = pd.DataFrame(expected)
    result = ep.eaesubtype_input(os.path.join(data_dir, input), name)
    pdt.assert_frame_equal(result, shouldbe)
