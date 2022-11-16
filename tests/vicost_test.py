import pytest
import os
import datetime
from .utils import cmd_runner

vicost_cmd = cmd_runner(["python", "vicost"])

@pytest.mark.parametrize(
    "options,expected",
    [
        (["--name"], "argument --name/-n: expected one argument"),
        ([], "the following arguments are required: --R1, --R2, --name/-n"),
        (["--name", "test1"], "the following arguments are required: --R1, --R2"),
        (["--name", "test1", "--outdir", 'outdir'],
        "the following arguments are required: --R1, --R2"),
        (["--R1", "test1", "--R2", "test"], "the following arguments are required: --name/-n")
    ])

def test_missing_args(options, expected):
    out, err, code = vicost_cmd(options)
    assert code == 2
    assert expected in err
    assert out == ''

def test_cmd(out_dir):
    """
    test_cmd
    """
    now = datetime.datetime.now()
    date = now.strftime("%Y%m%d")
    data_dir = "tests/data/"
    r1 = os.path.join(data_dir, "test1_trimmed.paired_R1.fq.gz")
    r2 = os.path.join(data_dir, "test1_trimmed.paired_R2.fq.gz")

    out, err, code = vicost_cmd(["--R1", r1, "--R2", r2, "--name", "test1", "--outdir", out_dir])

    assert code == 0

    expected_files = [
        "test1_2recAstxeaehly.txt",
        "test1_eaesubtype.tab",
        "test1_stecvir.tab",
        "test1_virbarcode_%s.tab" % (date),
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
    