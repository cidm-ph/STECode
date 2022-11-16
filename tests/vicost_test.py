import pytest
import os
import datetime
from .utils import cmd_runner

vicost_cmd = cmd_runner(["python", "vicost"])

@pytest.mar.parametrize(
    "options,expected",
    [
        (["--name"], "argument --name/-n: expected one argument"),
        ([], "the following arguments are required: --name/-n, --input/-i"),
        (["--name", "test1"], "the following arguments are required: --input/-i"),
        (["--name", "test1", "--outdir", 'outdir'],
    "the following arguments are required: --input/-i"),
])

def test_cmd(out_dir):
    """
    test_cmd
    """
    now = datetime.datetime.now()
    date = now.strftime("%Y%m%d")
    data_dir = "tests/data/"

    out, err, code = vicost_cmd(["--R1", "--R1", "--name", "test1", "--outdir", out_dir, "--input", data_dir])

    assert code == 0

    expected_files = [
        "test1_2recAstxeaehly.txt",
        "test1_eaesubtype.tab",
        "test1_stecvir.tab",
        "test1_virbarcode %s" % (date),
    ]

    for file in expected_files:
        output_file_path = os.path.join(out_dir, file)
        expected_file_path = os.path.join(data_dir, file)

        # check files exists
        assert os.path.isfile(output_file_path)
        with open(expected_file_path) as expected_file:
            expected_line = expected_file.readlines()
        with open(output_file_path) as output_file:
            # check file contents match
            assert expected_line == output_file.readlines()
    