# vicoSt

A pipeline for deployment of a Shiga Toxin-Producing Escherichia coli Virulence Barcode.

The STEC virulence barcode made up of 12 digits in the format of “XX-XX-XX-XX-XX-XX” reflecting their STEC virulence makeup. Where each set of two digits represented a particular virulence factor.

## Installation

PIP (Recommended)
```
pip install git+https://github.com/cidm-ph/vicoSt
```
GITHUB
```
git clone https://github.com/cidm-ph/vicoSt.git
```

## Usage
Example usage (mini version)

```
vicost --R1 $PATH/$R1.fq.gz --R2 $PATH/$R2.fq.gz --name $FILENAME --outdir $OUTDIR
```

FLAGS

```
--outdir, -o [Folder] optional folder to write output files to
--longread, -l turns on long read mode.
--R1, R1 fastq of sample (can be gzipped files)
--R2, R2 fastq of sample (can be gzipped files)
--name, -n [Name] name of the file you wish it to be
```

## Output

The first set of two digits represented the presence (to the subtype level) or absence of the eae gene. The next set of two digits represented inference of possible multiple, isogenic stx genes not assembled via short read sequencing. The last four sets of 2-mers each reflected the presence (to the subtype level) or absence of stx. This representation allowed up to four different stx operons to be captured, which is currently the maximum number observed both in vitro and in isolates.

## Dependencies
- Python == 3.8.13
- [BWA](https://sourceforge.net/projects/bio-bwa/) >= 0.7.27
- [Samtools](http://www.htslib.org/) >= 1.16.1
- [Abricate](https://github.com/tseemann/abricate) >= 1.0.1


## Nomenclature
(VI)rulence Bar(CO)ding for (ST)EC

## Associated Citations
Sim, E. M., Kim, R., Gall, M., Arnott, A., Howard, P., Valcanis, M., . . . Sintchenko, V. (2021). Added value of genomic surveillance of virulence factors in shiga toxin-producing escherichia coli in New South Wales, Australia. Frontiers in Microbiology, 12. [doi:10.3389/fmicb.2021.713724](https://www.frontiersin.org/articles/10.3389/fmicb.2021.713724/full)

## Licence
Copyright (C) 2022 Western Sydney Local Health District, NSW Health

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the [GNU General Public License](./LICENSE) along with this program. If not, see <https://www.gnu.org/licenses/>. 
