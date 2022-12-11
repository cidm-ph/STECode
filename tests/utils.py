"""
Written by @Wytamma
Modified for vicost by @LilWinx
"""

from dataclasses import dataclass
import subprocess
from typing import List, Tuple


@dataclass
class cmd_runner:
    """
    Commands for testing.
    """

    command: List[str]

    def __call__(self, args: List[str]) -> Tuple[str, str, int]:
        proc = subprocess.Popen(
            self.command + args,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        out, err = (output.decode("utf-8") for output in proc.communicate())
        return out, err, proc.returncode
