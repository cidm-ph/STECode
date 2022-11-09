from distutils.version import LooseVersion
import argparse
import logging
import shutil
import re
import sys
import subprocess

__version__ = "0.0.1"

class vicostargs:
    def arguments():
        # argparse
        parser = argparse.ArgumentParser(description="vicoSt")
        parser.add_argument(
            "--outdir", "-o", 
            help="Output directory to write to"
        )
        parser.add_argument(
            "--input", "-i", 
            help="Input directory or file",
            required=True
        )
        parser.add_argument(
            "--version",
            action="version",
            help="get SABRes version",
            version="SABRes v%s" % __version__,
        )
        parser.add_argument(
            '-h', '--help',
            action='help',
            help='Show this help message and exit'
        )
        vars(parser.parse_args())

def check_arguments(args):

def execute_command(command, check=True):
    logging.debug('command: %s', command)
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    # Manually decoding as subprocess.run decoding replaces \r with \n
    result.stdout = result.stdout.decode()
    result.stderr = result.stderr.decode()
    if check and result.returncode != 0:
        logging.critical('Failed to run command: %s', result.args)
        logging.critical('stdout: %s', result.stdout)
        logging.critical('stderr: %s', result.stderr)
        sys.exit(1)
    return result

def check_dependencies():
    logging.info('Checking dependencies')
    dependencies = {
        'bwa': {
            'vcommand': 'bwa -version',
            'vregex': re.compile(r'^blastn: (.+)\n'),
            'vrequired': '2.2.28'},
        'makeblastdb': {
            'vcommand': 'makeblastdb -version',
            'vregex': re.compile(r'^makeblastdb: (.+)\n'),
            'vrequired': '2.2.28'},
        'prodigal': {
            'vcommand': 'prodigal -v 2>&1',
            'vregex': re.compile(r'.*Prodigal V(.+?):'),
            'vrequired': '2.6.1'}
        }
    for dependency, version_data in dependencies.items():
        if not shutil.which(dependency):
            logging.critical('Could not find dependency %s' % dependency)
            sys.exit(1)
        result = execute_command(version_data['vcommand'], check=False)
        try:
            version = version_data['vregex'].search(result.stdout).group(1)
        except AttributeError:
            logging.critical('Unable to determine version for %s' % dependency)
            sys.exit(1)
        if LooseVersion(version) < LooseVersion(version_data['vrequired']):
            msg = '%s version %s or high is required'
            logging.critical(msg % (dependency, version_data['vrequired']))
            sys.exit(1)
        else:
            logging.debug('Found %s version %s' % (dependency, version))