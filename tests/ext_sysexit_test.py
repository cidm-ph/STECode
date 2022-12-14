import pytest
import os

def test_missing_R2():
    #importing the subprocess
    import subprocess

    # testing the exit of not applying an R2
    data_dir = "tests/data/"
    r1 = os.path.join(data_dir, "test1_trimmed.paired_R1.fq.gz")
    name = "test1"
    result = subprocess.run(["stecode", "--R1", r1, "-n", name], capture_output=True)

    #assert that capture output is matching the expected
    lines = result.stderr.splitlines()
    last_line = lines[-1]
    assert result.returncode == 1
    assert last_line == b'ERROR:root:R2 was not provided, please provide the paired reads'