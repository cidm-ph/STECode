"""
Welcome to the secondary script of vicoSt
Written by Winkie Fong - winkie.fong@health.nsw.gov.au

A mini version of vicoSt that runs directly on the output files generated from samtools, and abricate
"""

import os
import sys
import argparse
import logging
import gen_output as go

__version__ = "0.0.1"
logging.getLogger().setLevel(logging.INFO)

def mini_vicost():
    parser = argparse.ArgumentParser(description="Mini_vicoSt")
    parser.add_argument("--outdir", "-o", help="Output directory to write to")
    parser.add_argument("--path", "-p", help="Path to input files", required=True)
    parser.add_argument("--name", "-n", help="Name of samples", required=True)
    parser.add_argument("--longread", "-l", help="Samples were sequenced with long read sequencing", action="store_true")
    args = vars(parser.parse_args())

    if not args["outdir"]:
        indir = os.path.dirname(args["path"])
        if indir == "":
            args["outdir"] = "."
        else:
            args["outdir"] = indir

    logging.info(
        "Launching Mini-vicoSt v%s on %s and writing output files to directory %s"
        % (__version__, args["name"], args["outdir"])
    )

    file1 = os.path.join(args['path'], args["name"] + "_eaesubtype.tab")
    file2 = os.path.join(args['path'], args["name"] + "_2recAstxeaehly.txt")
    file3 = os.path.join(args['path'], args["name"] + "_stecvir.tab")

    check_files(file1)
    check_files(file2)
    check_files(file3)
    check_folders(args['outdir'])

    go.gen_output(
        file1,
        file2,
        file3,
        args["longread"],
        args['name'],
        args['outdir']
    )

def check_files(file):
    """
    Check input files if they exist and have contents
    """
    if os.path.isfile(file) == True and os.stat(file).st_size != 0:
        truemsg = file + " exists and not empty, continuing..."
        logging.info(truemsg)
    else:
        msg = file + " either file is does not exist or is empty, please check files. Exiting."
        logging.critical(msg)
        sys.exit(1)

def check_folders(folder):
    """
    Check the output folder if it exists, if not make new directory.
    """
    if os.path.exists(folder) == True:
        truemsg = folder + " output folder exists"
        logging.info(truemsg)
    else:
        os.makedirs(folder)
        msg = folder + " does not exist, making output folder"
        logging.info(msg)

if __name__ == "__main__":
    mini_vicost()