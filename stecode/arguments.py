import argparse

__version__ = "0.0.5"


def create_parser():
    parser = argparse.ArgumentParser(description="STECode", prog="stecode")
    parser.add_argument("--outdir", "-o", help="Output directory to write to")
    parser.add_argument("--R1", help="Path to R1 file")
    parser.add_argument("--R2", help="Path to R2 file")
    parser.add_argument("--fasta", "-f", help="Path to Fasta file")
    parser.add_argument(
        "--longread",
        "-l",
        action="store_true",
        help="Genome was assembled by long reads",
    )
    parser.add_argument(
        "--threads",
        "-t",
        nargs="?",
        const=4,
        type=int,
        help="Specify number of threads used",
    )
    parser.add_argument("--name", "-n", help="Name of sample", required=True)
    parser.add_argument(
        "--version",
        "-v",
        action="version",
        help="get STECode version",
        version=f"STECode v{__version__}",
    )
    """
    parser.add_argument(
        "--parallel",
        "-p",
        action="store_true",
        help="Parallelise the mapping portion if computer allows",
    )
    """
    return parser
