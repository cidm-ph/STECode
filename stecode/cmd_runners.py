"""
All the runner scripts and commands 
"""

from . import assists


def run_bwa(fq1, fq2, ref, name, outdir):
    """
    Run BWA and Samtools to generate the 2recAstxeae text file.
    """
    refname = ref.split("_")[2][:-6]
    command1 = f"bwa mem -t 8 {ref} {fq1} {fq2} | samtools view --threads 8 -b -S | samtools sort --threads 8 -o {outdir}/{name}_{refname}.sorted.bam"
    command2 = f"samtools index -b {outdir}/{name}_{refname}.sorted.bam"
    command3 = f"samtools coverage {outdir}/{name}_{refname}.sorted.bam > {outdir}/{name}_{refname}.txt && sed -n '2p' {outdir}/{name}_{refname}.txt >> {outdir}/{name}_2recAstxeae.txt"

    assists.run_cmd(command1)
    assists.run_cmd(command2)
    assists.run_cmd(command3)


def run_skesa(fq1, fq2, name, outdir):
    """
    Run Skesa for input into Abricate.
    """
    command4 = f"skesa --fastq {fq1},{fq2} --contigs_out {outdir}/{name}.contigs.fa"
    assists.run_cmd(command4)


def run_abricate(eaesub_db, stecvir_db, name, outdir):
    """
    Run Abricate to get the eaesubtype and stecvirulence
    """

    command5 = f"abricate --datadir {assists.stecode_db_dir} --db {eaesub_db} --minid 90.0 --mincov 90.0 {outdir}/{name}.contigs.fa > {outdir}/{name}_eaesubtype.tab"
    command6 = f"abricate --datadir {assists.stecode_db_dir} --mincov 21 --db {stecvir_db} {outdir}/{name}.contigs.fa > {outdir}/{name}_sfindAbricate.tab"
    assists.run_cmd(command5)
    assists.run_cmd(command6)


def run_solo_abricate(eaesub_db, stecvir_db, name, infile, outdir):
    """
    Run Abricate to get the eaesubtype and stecvirulence for FASTAs only
    """

    command5 = f"abricate --datadir {assists.stecode_db_dir} --db {eaesub_db} --minid 90.0 --mincov 90.0 {infile} > {outdir}/{name}_eaesubtype.tab"
    command6 = f"abricate --datadir {assists.stecode_db_dir} --mincov 21 --db {stecvir_db} {infile} > {outdir}/{name}_sfindAbricate.tab"
    assists.run_cmd(command5)
    assists.run_cmd(command6)
