"""
All the runner scripts and commands 
"""

from . import assists


def run_bwa(fq1, fq2, ref, name, outdir):
    """
    Run BWA and Samtools to generate the 2recAstxeae text file.
    """
    command1 = (
        "bwa mem -t 8 %s %s %s | samtools view --threads 8 -b -S | samtools sort --threads 8 -o %s/%s_2recAstxeae.sorted.bam"
        % (ref, fq1, fq2, outdir, name)
    )
    command2 = "samtools index -b %s/%s_2recAstxeae.sorted.bam" % (outdir, name)
    command3 = (
        "samtools coverage %s/%s_2recAstxeae.sorted.bam > %s/%s_2recAstxeae.txt"
        % (outdir, name, outdir, name)
    )
    assists.run_cmd(command1)
    assists.run_cmd(command2)
    assists.run_cmd(command3)


def run_skesa(fq1, fq2, name, outdir):
    """
    Run Skesa for input into Abricate.
    """
    command4 = "skesa --fastq %s,%s --contigs_out %s/%s.contigs.fa" % (
        fq1,
        fq2,
        outdir,
        name,
    )
    assists.run_cmd(command4)


def run_abricate(eaesub_db, stecvir_db, name, outdir):
    """
    Run Abricate to get the eaesubtype and stecvirulence
    """

    command5 = (
        "abricate --datadir %s --db %s --minid 90.0 --mincov 90.0 %s/%s.contigs.fa > %s/%s_eaesubtype.tab"
        % (assists.vicost_db_dir, eaesub_db, outdir, name, outdir, name)
    )
<<<<<<< HEAD
    command6 = "abricate -db %s --mincov 21 %s/%s.contigs.fa > %s/%s_sfindAbricate.tab" % (
=======
    command6 = "abricate --datadir %s --db %s %s/%s.contigs.fa > %s/%s_stecvir.tab " % (
        assists.vicost_db_dir,
>>>>>>> fa5e87d713c3b21012fada79fe2d7285fb0b5aa4
        stecvir_db,
        outdir,
        name,
        outdir,
        name,
    )
    assists.run_cmd(command5)
    assists.run_cmd(command6)
