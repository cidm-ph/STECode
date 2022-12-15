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
            "b'#Sequence_ID eae_sub iso_tox tox1 tox2 tox3 tox4 Virulence_barcode\\n       test1      06      CG   1A   2C   00   00 06-CG-1A-2C-00-00\\n'",
        ),
        (
            os.path.join(data_dir, "test1.contigs.fa"),
            "test1",
            False,
            "b'#Sequence_ID eae_sub iso_tox tox1 tox2 tox3 tox4 Virulence_barcode\\n       test1      06      DG   1A   2C   00   00 06-DG-1A-2C-00-00\\n'",
        ),
        (
            os.path.join(data_dir, "test2.contigs.fa"),
            "test2",
            True,
            "b'#Sequence_ID eae_sub iso_tox tox1 tox2 tox3 tox4 Virulence_barcode\\n       test2      03      CG   1A   00   00   00 03-CG-1A-00-00-00\\n'",
        ),
        (
            os.path.join(data_dir, "test2.contigs.fa"),
            "test2",
            False,
            "b'#Sequence_ID eae_sub iso_tox tox1 tox2 tox3 tox4 Virulence_barcode\\n       test2      03      DG   1A   00   00   00 03-DG-1A-00-00-00\\n'",
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

    output_file_path = os.path.join(out_dir, name + "_virbarcode_" + date + ".tab")

    # check files exists
    assert os.path.isfile(output_file_path)

    # compare lines
    assert str(result.stdout) == expected