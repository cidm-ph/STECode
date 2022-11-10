import os
import subprocess
import sys
import logging

def execute_command(command, check=True):
    logging.debug('command: %s', command)
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    result.stdout = result.stdout.decode()
    result.stderr = result.stderr.decode()
    if check and result.returncode != 0:
        logging.critical('Failed to run command: %s', result.args)
        logging.critical('stdout: %s', result.stdout)
        logging.critical('stderr: %s', result.stderr)
        sys.exit(1)
    return 

def check_files(file):
    if os.path.isfile(file) == True and os.stat(file).st_size != 0:
        truemsg = file + " exists and not empty, continuing..."
        logging.info(truemsg)
    else:
        msg = file + " either file is does not exist or is empty, please check files. Exiting."
        logging.critical(msg)
        sys.exit(1)

def check_folders(folder):
    if os.path.exists(folder) == True:
        truemsg = folder + " output folder exists"
        logging.info(truemsg)
    else:
        os.makedirs(folder)
        msg = folder + " does not exist, making output folder"
        logging.info(msg)

