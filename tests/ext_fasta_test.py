import pytest
import os
import datetime
import subprocess
import sys

data_dir = "tests/data/"


@pytest.mark.parametrize(
    "fasta, name, longread, expected",
    [
        (
            os.path.join(data_dir, "test1.contigs.fa"),
            "test1",
            True,
            "#Sequence_ID\teae_sub\tiso_tox\ttox1\ttox2\ttox3\ttox4\tVirulence_barcode\ntest1\t06\tCG\t1A\t2C\t00\t00\t06-CG-1A-2C-00-00\n",
        ),
        (
            os.path.join(data_dir, "test1.contigs.fa"),
            "test1",
            False,
            "#Sequence_ID\teae_sub\tiso_tox\ttox1\ttox2\ttox3\ttox4\tVirulence_barcode\ntest1\t06\tDG\t1A\t2C\t00\t00\t06-DG-1A-2C-00-00\n",
        ),
        (
            os.path.join(data_dir, "test2.contigs.fa"),
            "test2",
            True,
            "#Sequence_ID\teae_sub\tiso_tox\ttox1\ttox2\ttox3\ttox4\tVirulence_barcode\ntest2\t03\tCG\t1A\t00\t00\t00\t03-CG-1A-00-00-00\n",
        ),
        (
            os.path.join(data_dir, "test2.contigs.fa"),
            "test2",
            False,
            "#Sequence_ID\teae_sub\tiso_tox\ttox1\ttox2\ttox3\ttox4\tVirulence_barcode\ntest2\t03\tDG\t1A\t00\t00\t00\t03-DG-1A-00-00-00\n",
        ),
    ],
)
def test_cmd(out_dir, fasta, name, longread, expected):
    """
    test_cmd
    """
    now = datetime.datetime.now()
    date = now.strftime("%Y%m%d")
    if longread is False:
        result = subprocess.run(
            [
                "python3",
                "-m",
                "stecode",
                "-f",
                fasta,
                "-n",
                name,
                "-o",
                out_dir,
            ],
            capture_output=True,
        )
    else:
        result = subprocess.run(
            [
                "python3",
                "-m",
                "stecode",
                "-f",
                fasta,
                "-n",
                name,
                "-o",
                out_dir,
                "-l",
            ],
            capture_output=True,
        )

    assert result.returncode == 0

    output_file_path = os.path.join(
        out_dir, name, name + "_virbarcode_" + date + ".tab"
    )

    # check files exists
    assert os.path.isfile(output_file_path)

    # compare lines
    with open(output_file_path, "r") as outfile:
        out_results = outfile.read()
    assert out_results == expected
