import os
import subprocess
import sys
import logging
import shutil
import os.path
import pkg_resources

stecode_db_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "database")


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
    if cmd_exec == "abricate" or cmd_exec == "samtools" or cmd_exec == "skesa":
        vcmd = subprocess.run([cmd_exec, "--version"], capture_output=True, text=True)
        result = vcmd.stdout.splitlines()
        if cmd_exec == "abricate":
            version = " ".join(result).replace("abricate ", "")
            if pkg_resources.parse_version(version) < pkg_resources.parse_version(
                "1.0.0"
            ):
                logging.critical("Abricate version too old, please upgrade to v1.0.0+")
                sys.exit(1)
        if cmd_exec == "samtools":
            version = result[0].replace("samtools ", "")
            if pkg_resources.parse_version(version) < pkg_resources.parse_version(
                "1.10"
            ):
                logging.critical("Samtools version too old, please upgrade to v1.10.0+")
                sys.exit(1)
        if cmd_exec == "skesa":
            version = " ".join(result).replace("SKESA ", "")
    if cmd_exec == "bwa":
        vcmd = subprocess.run(["bwa"], capture_output=True, text=True)
        result = vcmd.stderr.splitlines()
        int_res = result[2].replace("Version: ", "")
        v_parts = int_res.split("-")
        version = "".join(v_parts[:1])

    if cmd_path is not None:
        msg = "Located " + cmd_exec + " " + version + " in " + cmd_path
        logging.info(msg)
    else:
        msg = cmd_exec + " was not found, please check installation on your device"
        logging.critical(msg)
        sys.exit(1)


def check_abricate():
    result = subprocess.run(
        ["abricate", "--list", "--datadir", stecode_db_dir],
        capture_output=True,
        text=True,
    )
    if result.returncode > 0:
        logging.critical("Abricate database is not prepared")
        logging.critical(
            "correct by running:   abricate --setupdb --datadir " + stecode_db_dir
        )
        sys.exit(1)
    dbs = [x.split("\t")[0] for x in result.stdout.splitlines()[1:]]
    if any(x not in dbs for x in ["eaesub", "stecodeDB"]):
        logging.critical("unable to find STECode databases")
        sys.exit(1)
