# pileup2vcf

## Synopsys

p2v is a mpileup-based variant caller that intends to allow you get everything and to filter variants according to your criterias.

## Dependencies

1) You should use samtools-0.1.18 in order to generate your pileup files. 

```
cd Dependencies
tar xjvf samtools-0.1.18.tar.bz2
cd samtools-0.1.18
make
( sudo make install )
```

2) You should install the following python libraries:
* scipy
* numpy
* fisher
* argparse
* matplotlib (this might require the installation of libpng-devel and freetype-devel)
* biopython
* pyfaidx

For simple installation, use `pip install <dependency name>` 
If not root: `pip install --user <dependency name>`

3) Generate a samtools file

```samtools mpileup -A -s -O -B -f <reference genome> <your bam file> > <output.mpileup>```

4) Use p2v

For preset clinical constitutionnal filters:

```python p2v -i <output from samtools mpileup> --outputPrefix <output prefix> --regions <bed file> --refGenome <reference genome> --clinics```


## Author

Yannick Boursin - 2015 & 2016 - Gustave Roussy

Until we choose a licence, please be advised that all rights are reserved.

## Contributors

* Birama Ndiaye - 2015 & 2016 - Gustave Roussy
* Alec Guyomard - 2016 - Gustave Roussy
* Jeremie Pagnac - 2015 - Gustave Roussy

## Licence

Yet to be chosen.
For now, please do not use this program without explicit authorization.
Please report any bug to yannick.boursin@gustaveroussy.fr

