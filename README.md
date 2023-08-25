cirCodAn
=======
**cirCodAn** (**cir**cular **Cod**ing sequence **An**notator) is a computational tool designed to predict CDS in circRNAs.

 - **cirCodAn** is an extension of [**CodAn**](https://github.com/pedronachtigall/CodAn), which was designed to predict CDS in linear RNAs.

 - The developmental repository of cirCodAn is available [here](https://github.com/denilsonfbar/cirCodAn-dev-v1).



## Quick Installation

Clone the cirCodAn repository:
```
git clone https://github.com/denilsonfbar/cirCodAn.git
```

Add the CirCodAn and CodAn/bin folder into your PATH:
```
export PATH=$PATH:$PWD/cicCodAn/
export PATH=$PATH:$PWD/cicCodAn/CodAn/bin/
```

Apply "execution permission" to executables if needed:
```
chmod +x $PWD/cicCodAn/*.py
chmod +x $PWD/cicCodAn/CodAn/bin/*
```

## Requirements

- [Python3](https://www.python.org/) and [Biopython](https://biopython.org/wiki/Download)
    - ```apt-get install python3-biopython```
- [Perl](https://www.perl.org/), [Bioperl](https://bioperl.org/) and [MCE](https://metacpan.org/release/MCE) (libmce-perl)
    - ```apt-get install bioperl libmce-perl```
- [CodAn]([https://www.perl.org/](https://github.com/pedronachtigall/CodAn))
- [ToPS](https://tops.sourceforge.net/)

:warning: **Conda environment installation**

If you takes advantage of [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) (i.e., one of the famous environment management system), you can follow the steps below:


Create Conda environment from ```environment.yml``` configuration file:
```
#create the enviroment
conda env create -f environment.yml
#add cirCodAn into your PATH and give execution permission
export PATH=$PATH:$PWD/cicCodAn/
chmod +x $PWD/cicCodAn/*.py
```

Activate the created environment:
```
conda activate circcodan
```

## Usage

Run the prediction example using Python interpreter:
```
python3 circodan.py -f example/circRNA_seqs.fa
```

OR

Give execution privileges to the cirCodAn script file:
```
chmod +x circodan.py
```

And run the example directly from the script file:
```
./circodan.py -f example/circRNA_seqs.fa
```

OR

If you need to run cirCodAn from anywhere on the file system, update the PATH variable:
```
export PATH=$PATH:path/to/cirCodAn
```

And run the example from anywhere on the file system:
```
circodan.py -f path/to/example/circRNA_seqs.fa
```

### Expected example output
```
2023-03-15 17:07:25 -> started cirCodAn v1.0
2023-03-15 17:07:26 -> prediction finished
Number of input sequences -> 138
Number of predicted CDSs  -> 104
GTF file with prediction annotation -> example/cirCodAn_output/CDS_predicted.gtf
Predicted CDS seqs FASTA file -> example/cirCodAn_output/CDS_predicted_seqs.fa
Predicted peptides FASTA file -> example/cirCodAn_output/CDS_predicted_seqs_aa.fa
```

The files with the outputs of cirCodAn execution are recorded at the addresses given.


## Manual installation

Install the following requirements:

- [Python3](https://www.python.org/), [Pandas](https://pandas.pydata.org/docs/getting_started/install.html) and [Biopython](https://biopython.org/wiki/Download)
- [Perl](https://www.perl.org/), [Bioperl](https://bioperl.org/) and [MCE](https://metacpan.org/release/MCE)
- [NCBI-BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279671/) (v2.9.0 or above)

Carry out of the same steps described in **Installation** section, except for creating and activating the Conda environment.


## Usage

```
Usage: circodan.py [options]

Options:
  -h, --help            show this help message and exit
  -f file, --file=file  Mandatory - input circRNAs file (FASTA format),
                        /path/to/circRNA_seqs.fa
  -o folder, --output=folder
                        Optional - path to output folder,
                        /path/to/output/folder/ if not declared, it will be
                        created at the circRNAs input folder
                        [default="cirCodAn_output"]
```

Basic example to find CDS in circRNA sequences:
```
python3 circodan.py -f example/circRNA_seqs.fa
```

## Reference

If you use or discuss circCodAn, please cite:

(circCodAn book chapter under peer review)

Nachtigall, P. G., Kashiwabara, A. Y., & Durham, A. M. (2021). CodAn: Predictive models for precise identification of coding regions in eukaryotic transcripts. Briefings in Bioinformatics, 22(3), 1–11. https://doi.org/10.1093/bib/bbaa045

## License

[GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)


## Contact

To report bugs, to ask for help and to give any feedback, please contact denilsonfbar@gmail.com
