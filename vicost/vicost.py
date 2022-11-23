"""
Welcome to the primary script of vicoSt
Written by Winkie Fong - winkie.fong@health.nsw.gov.au
"""

import argparse
import logging
import os
from vicost import assists
from vicost import cmd_runners
from vicost import gen_output as go

__version__ = "0.0.1"
logging.getLogger().setLevel(logging.INFO)

dependency_list = ["abricate", "samtools", "bwa", "skesa"]


def vicost():
    """
    Running order of vicoSt
    """
    parser = argparse.ArgumentParser(description="vicoSt")
    parser.add_argument("--outdir", "-o", help="Output directory to write to")
    parser.add_argument("--R1", help="Path to R1 file", required=True)
    parser.add_argument("--R2", help="Path to R2 file", required=True)
    parser.add_argument("--name", "-n", help="Name of sample", required=True)
    parser.add_argument(
        "--longread",
        "-l",
        help="Samples were sequenced with long read sequencing",
        action="store_true",
    )
    args = vars(parser.parse_args())

    if not args["outdir"]:
        indir = os.path.dirname(args["R1"])
        if indir == "":
            args["outdir"] = "."
        else:
            args["outdir"] = indir

    logging.info(
        "Launching vicoSt v%s on %s and writing output files to directory %s",
        __version__,
        args["name"],
        args["outdir"],
    )
    logging.info("Checking installs of dependencies")
    for dependency in dependency_list:
        assists.check_dependencies(dependency)

    ref = os.path.join(os.path.dirname(__file__), "database/stx_eae_hly_recA.fasta")

    # checking file integrity and existence of output directory
    assists.check_files(args["R1"])
    assists.check_files(args["R2"])
    assists.check_folders(args["outdir"])

    # Run bwa, samtools, skesa and abricate
    cmd_runners.run_bwa(args["R1"], args["R2"], ref, args["name"], args["outdir"])
    cmd_runners.run_skesa(args["R1"], args["R2"], args["name"], args["outdir"])
    cmd_runners.run_abricate("eaesub", "stecvir", args["name"], args["outdir"])

    # vicost portion - file check
    file1 = os.path.join(args["outdir"], args["name"] + "_eaesubtype.tab")
    file2 = os.path.join(args["outdir"], args["name"] + "_2recAstxeaehly.txt")
    file3 = os.path.join(args["outdir"], args["name"] + "_stecvir.tab")

    assists.check_files(file1)
    assists.check_files(file2)
    assists.check_files(file3)

    # run vicost
    go.gen_output(file1, file2, file3, args["longread"], args["name"], args["outdir"])
    logging.info(
        "Complete :D please check %s for the STEC barcode for your sample",
        args["outdir"]
    )


if __name__ == "__main__":
    vicost()
