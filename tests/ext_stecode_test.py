import pytest
import os
import datetime
import subprocess

data_dir = "tests/data/"


@pytest.mark.parametrize(
    "r1, r2, fasta, name",
    [
        (
            os.path.join(data_dir, "test1_trimmed.paired_R1.fq.gz"),
            os.path.join(data_dir, "test1_trimmed.paired_R2.fq.gz"),
            os.path.join(data_dir, "test1.contigs.fa"),
            "test1",
        ),
        (
            os.path.join(data_dir, "test2_trimmed.paired_R1.fq.gz"),
            os.path.join(data_dir, "test2_trimmed.paired_R2.fq.gz"),
            os.path.join(data_dir, "test2.contigs.fa"),
            "test2",
        ),
    ],
)
def test_cmd(out_dir, r1, r2, fasta, name):
    """
    test_cmd
    """
    now = datetime.datetime.now()
    date = now.strftime("%Y%m%d")

    result = subprocess.run(
        ["stecode", "--R1", r1, "--R2", r2, "-f", fasta, "-n", name, "-o", out_dir],
        capture_output=True,
    )

    output_file_path = os.path.join(out_dir, name + "_virbarcode_" + date + ".tab")
    expected_file_path = os.path.join(data_dir, name + "_virbarcode.tab")

    # check files exists
    assert os.path.isfile(output_file_path)

    # compare lines
    with open(expected_file_path) as expected_file:
        expected_line = expected_file.readlines()
    with open(output_file_path) as output_file:
        # check file contents match
        assert expected_line == output_file.readlines()

    # compare the final line with expected output
    lines = result.stderr.splitlines()
    last_line = lines[-1]
    assert result.returncode == 0
    assert (
        str(last_line)
        == "b'INFO:root:Complete :D please check "
        + str(out_dir)
        + " for the STEC barcode for your sample'"
    )
