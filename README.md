![STECode_Brand_large_crop](https://user-images.githubusercontent.com/93765714/224181287-feba57a1-c336-48f1-a304-9a1ee7eb7464.png)
# STECode

A pipeline for deployment of a Shiga Toxin-Producing Escherichia coli Virulence Barcode.

The STEC virulence barcode made up of 12 digits in the format of “XX-XX-XX-XX-XX-XX” reflecting their STEC virulence makeup. Where each set of two digits represented a particular virulence factor.

## Installation

PIP (Recommended)
```
pip install git+https://github.com/cidm-ph/STECode
```
GITHUB
```
git clone https://github.com/cidm-ph/STECode.git
```

Then initialise the abricate database with
```
abricate --setupdb --datadir $STECode/stecode/database
```

## Usage
Example usage

```
stecode --R1 $PATH/$R1.fq.gz --R2 $PATH/$R2.fq.gz --name $FILENAME --outdir $OUTDIR
```

FLAGS

```
--outdir, -o [PATH]             optional folder to write output files to
--threads, -t [INT]             specify number of threads used (default = 4)
--R1 [PATH]                     R1 fastq of sample (can be gzipped files)
--R2 [PATH]                     R2 fastq of sample (can be gzipped files)
--fasta, -f [PATH]              optional fasta file which will skip SKESA, can be used in conjunction with --longread.              
--longread, -l                  turns on long read mode.
--name, -n [STR]                name of the file you wish it to be [REQUIRED!]
--version, -v                   print version
```

SKESA Genome Assembly is the longest portion of this pipeline, so if you already have a genome assembly you can bypass SKESA by supplying a FASTA file. A FASTA only input can also be performed however, the second 'XX' will not show isogenic stx genes.

## Output
A few files are coalesced from mapping and abricate into a virulence barcode. 

The first set of two digits represented the presence (to the subtype level) or absence of the eae gene. The next set of two digits represented inference of possible multiple, isogenic stx genes not assembled via short read sequencing. The last four sets of 2-mers each reflected the presence (to the subtype level) or absence of stx. This representation allowed up to four different stx operons to be captured, which is currently the maximum number observed both in vitro and in isolates.

The result of your barcode will appear on the console, log file and its own file (--name_virbarcode_YYYYMMDD.tab).

If a discrepancy between the mapping and abricate is found the program will stop and tell you to look at the raw output files. Raw output files that are most useful include:
- sfindAbricate 
- eaesubtype
- targetstx

## Dependencies
- Python == 3.8.13
- [BWA](https://sourceforge.net/projects/bio-bwa/) >= 0.7.17
- [Samtools](http://www.htslib.org/) >= 1.16.1
- [Abricate](https://github.com/tseemann/abricate) >= 1.0.1
- [Skesa](https://github.com/ncbi/SKESA) >= 2.4.0


## Associated Citations
Sim, E. M., Kim, R., Gall, M., Arnott, A., Howard, P., Valcanis, M., . . . Sintchenko, V. (2021). Added value of genomic surveillance of virulence factors in shiga toxin-producing escherichia coli in New South Wales, Australia. Frontiers in Microbiology, 12. [doi:10.3389/fmicb.2021.713724](https://www.frontiersin.org/articles/10.3389/fmicb.2021.713724/full)

## Licence
Copyright (C) 2022 Western Sydney Local Health District, NSW Health

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the [GNU General Public License](./LICENSE) along with this program. If not, see <https://www.gnu.org/licenses/>. 
