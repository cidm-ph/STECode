import pytest
import os
from stecode import gen_output as go
import pandas as pd
import pandas.testing as pdt

data_dir = "tests/data/"
heirarchy = ["12", "02", "01", "T2", "T1", "00"]


@pytest.mark.parametrize(
    "stfile, recAfile, virfile, name, expected",
    [
        (
            "test1_eaesubtype.tab",
            "test1_targetstx.txt",
            "test1_sfindAbricate.tab",
            "test1",
            {
                "#Sequence_ID": ["test1"],
                "eae_sub": ["06"],
                "iso_tox": ["02"],
                "tox1": ["1A"],
                "tox2": ["2C"],
                "tox3": ["00"],
                "tox4": ["00"],
                "Virulence_barcode": ["06-02-1A-2C-00-00"],
            },
        ),
        (
            "test2_eaesubtype.tab",
            "test2_targetstx.txt",
            "test2_sfindAbricate.tab",
            "test2",
            {
                "#Sequence_ID": ["test2"],
                "eae_sub": ["03"],
                "iso_tox": ["00"],
                "tox1": ["1A"],
                "tox2": ["00"],
                "tox3": ["00"],
                "tox4": ["00"],
                "Virulence_barcode": ["03-00-1A-00-00-00"],
            },
        ),
    ],
)
def test_gen_output(stfile, recAfile, virfile, name, expected):
    shouldbe = pd.DataFrame(expected)
    shouldbe["iso_tox"] = pd.Categorical(
        shouldbe["iso_tox"], ordered=True, categories=heirarchy
    )
    result = go.merge_all_NNs(
        os.path.join(data_dir, stfile),
        os.path.join(data_dir, recAfile),
        os.path.join(data_dir, virfile),
        True,
        name,
        False,
    )
    print(shouldbe)
    print(result)
    pdt.assert_frame_equal(result, shouldbe)
