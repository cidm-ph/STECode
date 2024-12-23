"""
Welcome to the primary script of STECode
Written by Winkie Fong - winkie.fong@health.nsw.gov.au
"""

import datetime
import logging
import os
import sys
import warnings
import glob
import concurrent.futures
from stecode import arguments
from stecode import assists
from stecode import cmd_runners
from stecode import gen_output as go


__version__ = "1.1.0"
logging.getLogger().setLevel(logging.INFO)
warnings.simplefilter(action="ignore", category=FutureWarning)
formatter = logging.Formatter(
    "STECode:%(levelname)s:%(asctime)s: %(message)s", datefmt="%y/%m/%d %I:%M:%S %p"
)

dependency_list = ["abricate", "samtools", "bwa", "skesa"]
ref_list = []


def stecode():
    """
    Running order of STECode
    """
    parser = arguments.create_parser()  # pylint: disable=E1101
    args = parser.parse_args()
    is_assembly = bool(args.fasta is not None)
    is_reads = bool(args.R1 is not None)
    now = datetime.datetime.now()
    date = now.strftime("%Y%m%d")
    logger = logging.getLogger()

    # set outdir defaults - if no outdir is set, it will default to either the fasta or R1 location
    if args.outdir is None and args.fasta is not None:
        default = os.path.dirname(args.fasta)
        outdir = default
    elif args.outdir is None and args.R1 is not None:
        default = os.path.dirname(args.R1)
        outdir = default
    else:
        outdir = args.outdir

    # force creation of new folder within set outdir
    maindir = outdir + "/" + args.name
    newdir = maindir + "/bams"
    folder_exists = os.path.exists(newdir)
    if not folder_exists:
        os.makedirs(newdir)

    # errorlog
    errorlog = os.path.join(maindir, args.name + "_stecode_" + date + ".log")

    stdout_handler = logging.StreamHandler(sys.stdout)
    file_handler = logging.FileHandler(errorlog, mode="w+")
    for handler in [stdout_handler, file_handler]:
        handler.setLevel(logging.INFO)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    # log of folder creation
    if not folder_exists:
        logging.info("%s does not exist, creating directory.", newdir)
    else:
        logging.info("%s exists, skipping directory creation", newdir)

    # cmd checks
    if is_reads is True:
        if args.R2 is None:
            logging.error("R2 was not provided, please provide the paired reads")
            sys.exit(1)

    # set threads defaults - if no threads are set, it will default to 4 threads
    if args.threads is None:
        default_threads = 4
    else:
        default_threads = args.threads

    # launch line
    logging.info(
        "Launching STECode v%s on %s and writing output files to directory %s using %s threads",
        __version__,
        args.name,
        outdir,
        default_threads,
    )

    # checking all the versions and installations of dependencies.
    logging.info("Checking installs of dependencies")
    for dependency in dependency_list:
        assists.check_dependencies(dependency)
    if "abricate" in dependency_list:
        assists.check_abricate()

    # checking file integrity and existence of output directory
    if all(item is not None for item in [args.fasta, args.R1, args.R2]):
        assists.check_files(args.R1)
        assists.check_files(args.R2)
        assists.check_files(args.fasta)
        logging.info("Found fasta, R1 and R2, skipping Skesa")

        # skip skesa
        if args.longread is True:
            logging.info(
                "Running only Abricate on the long read assembled genome, iso_tox will be CG"
            )
            cmd_runners.run_solo_abricateLR(
                "eaesub", "stecodeDB", args.name, args.fasta, outdir
            )

        else:
            ref_path = os.path.join(os.path.dirname(__file__), "database/stxrecaeae/")
            ref = "STECode_normalisation_stxrecAeae.fasta"
            cmd_runners.run_bwa(outdir, ref_path + ref, default_threads)
            subref_path = os.path.join(newdir, args.name + "_stxrecAeae.txt")
            subref_fa = cmd_runners.get_subref(
                ref_path, subref_path, os.path.join(outdir, args.name)
            )
            cmd_runners.run_bwa_index(subref_fa)
            cmd_runners.run_bwa(outdir, subref_fa, default_threads)
            cmd_runners.run_solo_abricate(
                "eaesub", "stecodeDB", args.name, args.fasta, outdir
            )

    elif is_assembly is True and is_reads is False:
        assists.check_files(args.fasta)
        if args.longread is True:
            logging.info(
                "Running only Abricate on long read assembled genomes, iso_tox will be CG"
            )
            cmd_runners.run_solo_abricateLR(
                "eaesub", "stecodeDB", args.name, args.fasta, outdir
            )
        else:
            logging.info(
                "Running only Abricate on short read draft genomes, iso_tox will be DG"
            )
            cmd_runners.run_solo_abricate(
                "eaesub", "stecodeDB", args.name, args.fasta, outdir
                )

    else:
        assists.check_files(args.R1)
        assists.check_files(args.R2)
        assists.check_folders(outdir)

        # Run bwa, samtools, skesa and abricate
        ref_path = os.path.join(os.path.dirname(__file__), "database/stxrecaeae/")
        ref = "STECode_normalisation_stxrecAeae.fasta"
        cmd_runners.run_bwa(outdir, ref_path + ref, default_threads)
        subref_path = os.path.join(newdir, args.name + "_stxrecAeae.txt")
        subref_fa = cmd_runners.get_subref(
            ref_path, subref_path, os.path.join(outdir, args.name)
        )

        """
        if args.parallel is True:
            with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
                future_to_bam = {
                    executor.submit(
                        cmd_runners.run_bwa,
                        outdir,
                        ref_path + subref + ".fasta",
                        default_threads,
                    ): subref
                    for subref in subref_list
                }
                for future in concurrent.futures.as_completed(future_to_bam):
                    bam = future_to_bam[future]
                    try:
                        data = future.result()
                    except Exception as exc:
                        logging.error("%s generated an exception: %s", bam, exc)
        else:
            for subref in subref_list:
                cmd_runners.run_bwa(
                    outdir, ref_path + subref + ".fasta", default_threads
                )
        cmd_runners.combine_stxrecaeae(args.name, outdir)
        """
        cmd_runners.run_bwa_index(subref_fa)
        cmd_runners.run_bwa(outdir, subref_fa, default_threads)
        cmd_runners.run_skesa(args.R1, args.R2, default_threads / 2, args.name, outdir)
        cmd_runners.run_abricate("eaesub", "stecodeDB", args.name, outdir)

    # stecode portion - file check
    file1 = os.path.join(maindir, args.name + "_eaesubtype.tab")
    file3 = os.path.join(maindir, args.name + "_sfindAbricate.tab")

    assists.check_files(file1)
    assists.check_files(file3)
    if args.longread is True:
        file2 = "skip"
    elif is_assembly is True and is_reads is False:
        file2 = "skip"
    else:
        file2 = os.path.join(newdir, args.name + "_targetstx.txt")
        assists.check_files(file2)
        go.pre_merge_check(file2, file3, args.longread)

    # clean-up

    for filepath in glob.glob(
        os.path.join(maindir, "STECode_normalisation_targetstx.fasta*")
    ):
        os.remove(filepath)

    # run stecode
    go_df = go.merge_all_NNs(file1, file2, file3, is_reads, args.name, args.longread)
    go.gen_output(args.name, outdir, go_df)
    logging.info(
        "Complete :D we have also made it into a file, please check %s for the STEC barcode for your sample",
        outdir,
    )


if __name__ == "__main__":
    stecode()
