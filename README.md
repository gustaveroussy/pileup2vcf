# pileup2vcf

## Synopsys

p2v is a mpileup-based variant caller that intends to allow you get everything and to filter variants according to your criterias.

## Dependencies

1) You should use samtools-0.1.18 in order to generate your pileup files. 

cd Dependencies
tar xjvf *.bz2
cd samtools*
make
( sudo make install )

2) You should install the following python libraries:
* scipy
* numpy
* fisher
* argparse
* matplotlib (this might require 
* biopython
* pyfaidx

For simple installation, use pip install <dependency name> 
If not root: pip install --user <dependency name>

3) Generate a samtools file

samtools mpileup -A -s -O -B -o <output> -f <reference genome> <your bam file>

4) Use p2v

For preset clinical constitutionnal filters:

python p2v -i <output from samtools mpileup> --outputPrefix <output prefix> --regions <bed file> --refGenome <reference genome> --clinics


## Author

Yannick Boursin

## Contributors

Birama Ndiaye
Alec Guyomard
Jeremie Pagnac

## Licence

Yet to be chosen

