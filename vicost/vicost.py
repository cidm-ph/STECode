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

__version__ = "0.0.2"
logging.getLogger().setLevel(logging.INFO)

dependency_list = ["abricate", "samtools", "bwa", "skesa"]


def vicost():
    """
    Running order of vicoSt
    """
    parser = argparse.ArgumentParser(description="vicoSt")
    parser.add_argument("--outdir", "-o", help="Output directory to write to")
    parser.add_argument("--R1", help="Path to R1 file")
    parser.add_argument("--R2", help="Path to R2 file")
    parser.add_argument("--fasta", "-f", help="Path to Fasta file")
    parser.add_argument("--longread", "-")
    parser.add_argument("--name", "-n", help="Name of sample", required=True)
    parser.add_argument(
        "--version",
        "-v",
        action="version",
        help="get vicoSt version",
        version="vicoSt v%s" % __version__,
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
    if "abricate" in dependency_list:
        assists.check_abricate()

    ref = os.path.join(os.path.dirname(__file__), "database/stx_recA_eae.fasta")
    is_assembly = bool(args["fasta"] is not None)

    # checking file integrity and existence of output directory
    if is_assembly == True:
        assists.check_files(args["fasta"])

        # run only abricate
        cmd_runners.run_solo_abricate("eaesub", "stecfinder", args["name"], args['fasta'], args["outdir"])
    else:
        assists.check_files(args["R1"])
        assists.check_files(args["R2"])
        assists.check_folders(args["outdir"])

        # Run bwa, samtools, skesa and abricate
        cmd_runners.run_bwa(args["R1"], args["R2"], ref, args["name"], args["outdir"])
        cmd_runners.run_skesa(args["R1"], args["R2"], args["name"], args["outdir"])
        cmd_runners.run_abricate("eaesub", "stecfinder", args["name"], args["outdir"])

    # vicost portion - file check
    file1 = os.path.join(args["outdir"], args["name"] + "_eaesubtype.tab")
    file3 = os.path.join(args["outdir"], args["name"] + "_sfindAbricate.tab")

    assists.check_files(file1)
    assists.check_files(file3)

    if is_assembly is False:
        file2 = os.path.join(args["outdir"], args["name"] + "_2recAstxeae.txt")
        assists.check_files(file2)
    else:
        file2 = "skip"

    # run vicost
    go_df = go.merge_all_NNs(
        file1, file2, file3, is_assembly, args["name"], args["longread"]
    )
    go.gen_output(args["name"], args["outdir"], go_df)
    logging.info(
        "Complete :D please check %s for the STEC barcode for your sample",
        args["outdir"],
    )


if __name__ == "__main__":
    vicost()
