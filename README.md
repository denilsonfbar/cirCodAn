circCodAn
=======
**circCodAn** (**circ**ular **Cod**ing sequence **An**notator) is a computational tool designed to CDS prediction in circRNAs.

**circCodAn** is an extension of the **CodAn** tool, available in [https://github.com/pedronachtigall/CodAn](https://github.com/pedronachtigall/CodAn).


## Installation

Clone the circCodAn repository:
```
git clone https://github.com/denilsonfbar/circCodAn.git
```

Enter the circCodAn folder:
```
cd circCodAn
```

Change the permissions of the ```CodAn/bin/``` folder, to allow the creation of temporary files necessary for running the tool:
```
chmod -R 777 CodAn/bin/
```

Create Conda environment from ```environment.yml``` configuration file:
```
conda env create -f environment.yml
```

Activate the created environment:
```
conda activate circcodan_env
```


## Running

Run the prediction example using Python interpreter:
```
python3 circ-codan.py -f example/circRNA_seqs.fa
```

OR

Give execution privileges to the circCodAn script file:
```
chmod +x circ-codan.py
```

And run the example directly from the script file:
```
./circ-codan.py -f example/circRNA_seqs.fa
```

OR

If you need to run circCodAn from anywhere on the file system, update the PATH variable:
```
export PATH=$PATH:path/to/circCodAn
```

And run the example from anywhere on the file system:
```
circ-codan.py -f path/to/example/circRNA_seqs.fa
```

### Expected example output
```
2023-03-15 17:07:25 -> started circCodAn v1.0
2023-03-15 17:07:26 -> prediction finished
Number of input sequences -> 138
Number of predicted CDSs  -> 86
GTF file with prediction annotation -> example/circCodAn_output/CDS_prediction.gtf
Predicted CDS seqs FASTA file -> example/circCodAn_output/CDS_predicted_seqs.fa
Predicted peptides FASTA file -> example/circCodAn_output/CDS_predicted_seqs_aa.fa
```

The files with the outputs of circCodAn execution are recorded at the addresses given.


## Manual installation

Install the following requirements:

- [Python3](https://www.python.org/), [Pandas](https://pandas.pydata.org/docs/getting_started/install.html) and [Biopython](https://biopython.org/wiki/Download)
- [Perl](https://www.perl.org/), [Bioperl](https://bioperl.org/) and [MCE](https://metacpan.org/release/MCE)
- [NCBI-BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279671/) (v2.9.0 or above)

Carry out of the same steps described in **Installation** section, except for creating and activating the Conda environment.


## Usage

```
Usage: circ-codan.py [options]

Options:
  -h, --help            show this help message and exit
  -f file, --file=file  Mandatory - input circRNAs file (FASTA format),
                        /path/to/circRNA_seqs.fa
  -o folder, --output=folder
                        Optional - path to output folder,
                        /path/to/output/folder/ if not declared, it will be
                        created at the circRNAs input folder
                        [default="circCodAn_output"]
```

Basic example to find CDS in circRNA sequences:
```
python3 circ-codan.py -f example/circRNA_seqs.fa
```

## Reference

If you use or discuss circCodAn, please cite:

(not published yet)


## License

[GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)


## Contact

To report bugs, to ask for help and to give any feedback, please contact denilsonfbar@gmail.com
