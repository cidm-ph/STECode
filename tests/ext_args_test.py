import pytest
import os
import datetime
from .utils import cmd_runner

vicost_cmd = cmd_runner(["python", "-m", "vicost"])


@pytest.mark.parametrize(
    "options,expected",
    [
        ([], "the following arguments are required: --name/-n"),
    ],
)
def test_missing_args(options, expected):
    out, err, code = vicost_cmd(options)
    print(out, err, code)
    assert code == 2
    assert expected in err
    assert out == ""
