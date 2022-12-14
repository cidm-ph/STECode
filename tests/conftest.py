"""
Written by @Wytamma
Modified for stecode by @LilWinx
"""

import tempfile
import pytest


@pytest.fixture()
def out_dir():
    "conftest"
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir
