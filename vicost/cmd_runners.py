"""
All the runner scripts and commands 
"""

import assists

def recAstxeaehly_runner(fq1, fq2, ref, name, outdir):
    """
    ITS NOT EMPTY
    """
    command = 'bwa mem -t 30 %s %s %s > %s/%s.sam' % (ref, fq1, fq2, outdir, name)
    assists.execute_command(command)
