import pytest
import os
import subprocess

data_dir = "tests/data/"
r1 = os.path.join(data_dir, "test1_trimmed.paired_R1.fq.gz")
r2 = os.path.join(data_dir, "test1_trimmed.paired_R2.fq.gz")
fasta = os.path.join(data_dir, "test1.contigs.fa")
name = "test1"


@pytest.mark.parametrize(
    "options,expected",
    [
        ([], b"stecode: error: the following arguments are required: --name/-n"),
        (
            ["--R1", r1],
            b"stecode: error: the following arguments are required: --name/-n",
        ),
    ],
)
@pytest.mark.skipif(
    not (os.path.exists(r1) and os.path.exists(r2)), reason="Test data not available"
)
def test_missing_args(options, expected):
    # testing the exit of not applying an R2
    result = subprocess.run(["python3", "-m", "stecode"] + options, capture_output=True)

    # assert that capture output is matching the expected
    lines = result.stderr.splitlines()
    last_line = lines[-1]
    assert result.returncode == 2
    assert last_line == expected
