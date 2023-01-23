"""
All the runner scripts and commands 
"""

from . import assists


def run_bwa(fq1, fq2, ref, name, outdir):
    """
    Run BWA and Samtools to generate the 2recAstxeae text file.
    """
    refname = ref.split("_")[2][:-6]
    bamdir = outdir + "/" + name + "/bams"
    stuffdir = outdir + "/" + name
    command1 = f"bwa mem -t 8 {ref} {fq1} {fq2} | samtools view --threads 8 -b -S | samtools sort --threads 8 -o {bamdir}/{name}_{refname}.sorted.bam"
    command2 = f"samtools index -b {bamdir}/{name}_{refname}.sorted.bam"
    command3 = f"samtools coverage {bamdir}/{name}_{refname}.sorted.bam > {bamdir}/{name}_{refname}.txt"

    assists.run_cmd(command1)
    assists.run_cmd(command2)
    assists.run_cmd(command3)


def run_skesa(fq1, fq2, name, outdir):
    """
    Run Skesa for input into Abricate.
    """
    stuffdir = outdir + "/" + name
    command4 = f"skesa --fastq {fq1},{fq2} --contigs_out {stuffdir}/{name}.contigs.fa"
    assists.run_cmd(command4)


def run_abricate(eaesub_db, stecvir_db, name, outdir):
    """
    Run Abricate to get the eaesubtype and stecvirulence
    """
    stuffdir = outdir + "/" + name
    command5 = f"abricate --datadir {assists.stecode_db_dir} --db {eaesub_db} --minid 90.0 --mincov 90.0 {outdir}/{name}.contigs.fa > {stuffdir}/{name}_eaesubtype.tab"
    command6 = f"abricate --datadir {assists.stecode_db_dir} --mincov 21 --db {stecvir_db} {stuffdir}/{name}.contigs.fa > {stuffdir}/{name}_sfindAbricate.tab"
    assists.run_cmd(command5)
    assists.run_cmd(command6)


def run_solo_abricate(eaesub_db, stecvir_db, name, infile, outdir):
    """
    Run Abricate to get the eaesubtype and stecvirulence for FASTAs only
    """
    stuffdir = outdir + "/" + name
    command5 = f"abricate --datadir {assists.stecode_db_dir} --db {eaesub_db} --minid 90.0 --mincov 90.0 {infile} > {stuffdir}/{name}_eaesubtype.tab"
    command6 = f"abricate --datadir {assists.stecode_db_dir} --mincov 21 --db {stecvir_db} {infile} > {stuffdir}/{name}_sfindAbricate.tab"
    assists.run_cmd(command5)
    assists.run_cmd(command6)
