import os
import subprocess
import sys
import logging
import shutil
import os.path

vicost_db_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "database")


def run_cmd(command):
    """
    Run commands with error outputs.
    """
    logging.info("Running command: %s", command)
    result = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
    )
    result.stdout = result.stdout.decode()
    result.stderr = result.stderr.decode()
    if result.returncode != 0:
        logging.critical("Failed to run command: %s", result.args)
        logging.critical("stdout: %s", result.stdout)
        logging.critical("stderr: %s", result.stderr)
        sys.exit(1)
    return result


def check_files(file):
    """
    Check input files if they exist and have contents
    """
    if os.path.isfile(file) is True and os.stat(file).st_size != 0:
        truemsg = file + " exists and not empty, continuing..."
        logging.info(truemsg)
    else:
        msg = (
            file
            + " either file is does not exist or is empty, please check files. Exiting."
        )
        logging.critical(msg)
        sys.exit(1)


def check_folders(folder):
    """
    Check the output folder if it exists, if not make new directory.
    """
    if os.path.exists(folder) is True:
        truemsg = folder + " output folder exists"
        logging.info(truemsg)
    else:
        os.makedirs(folder)
        msg = folder + " does not exist, making output folder"
        logging.info(msg)


def check_dependencies(cmd_exec):
    cmd_path = shutil.which(cmd_exec)
    if cmd_path is not None:
        msg = "Located " + cmd_exec + " in " + cmd_path
        logging.info(msg)
    else:
        msg = cmd_exec + " was not found, please check installation on your device"
        logging.critical(msg)
        sys.exit(1)


def check_abricate():
    result = subprocess.run(
        ["abricate", "--list", "--datadir", vicost_db_dir],
        capture_output=True,
        text=True,
    )
    if result.returncode > 0:
        logging.critical("Abricate database is not prepared")
        logging.critical(
            "correct by running:   abricate --setupdb --datadir " + vicost_db_dir
        )
        sys.exit(1)
    dbs = [x.split("\t")[0] for x in result.stdout.splitlines()[1:]]
    if any(x not in dbs for x in ["eaesub", "stecfinder"]):
        logging.critical("unable to find vicoSt databases")
        sys.exit(1)
