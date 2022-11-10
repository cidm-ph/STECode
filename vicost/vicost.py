"""
Welcome to the primary script of vicoSt
Written by Winkie Fong - winkie.fong@health.nsw.gov.au
"""

import os
import sys
import argparse
import logging
import assists
import cmd_runners

__version__ = "0.0.1"

def vicost():
    parser = argparse.ArgumentParser(description="vicoSt")
    parser.add_argument("--outdir", "-o", help="Output directory to write to")
    parser.add_argument("--R1", help="Path to R1 file", required=True)
    parser.add_argument("--R2", help="Path to R2 file", required=True)
    parser.add_argument("--name", "-n", help="Name of sample", required=True)
    parser.add_argument("--longread", "-l", help="Samples were sequenced with long read sequencing", action="store_true")
    args = vars(parser.parse_args())

    if not args["outdir"]:
        indir = os.path.dirname(args["R1"])
        if indir == "":
            args["outdir"] = "."
        else:
            args["outdir"] = indir

    logging.info(
        "Launching vicoSt v%s on %s and writing output files to directory %s"
        % (__version__, args["name"], args["outdir"])
    )
    ref = os.path.join(os.path.dirname(__file__), "database/stx_eae_hly_recA.fasta")

    assists.check_files(args['R1'])
    assists.check_files(args['R2'])
    assists.check_folders(args['outdir'])

    cmd_runners.recAstxeaehly_runner(
        args['R1'],
        args['R2'],
        ref,
        args['name'],
        args['outdir']
    )

