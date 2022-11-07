"""
Welcome to the secondary script of stecode
Written by Winkie Fong - winkie.fong@health.nsw.gov.au

A mini version of stecode that runs directly on the output files generated from samtools, and abricate
"""

import os
import sys
import argparse

__version__ = "0.0.1"

def stecode():
    parser = argparse.ArgumentParser(description="MiniStecode")
    parser.add_argument("--outdir", "-o", help="Output directory to write to")
    parser.add_argument("--input", "-i", help="Input directory or file", required=True)
    args = vars(parser.parse_args())

    print(
        "Launching MiniStecode v%s on %s and writing output files to directory %s"
        % (__version__, args["input"], args["outdir"])
    )


def check_files(args):
    