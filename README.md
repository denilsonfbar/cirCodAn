<div align="center">
<center>

 ![cirCodAn_logo](/circodan_logo.png)

</center>
</div>

cirCodAn
=======
**cirCodAn** (**cir**cular **Cod**ing sequence **An**notator) is a computational tool designed to predict CDS in circRNAs.

 - **cirCodAn** is an extension of [**CodAn**](https://github.com/pedronachtigall/CodAn), which was designed to predict CDS in linear RNAs.


## Quick Installation

- Clone the cirCodAn repository:
   - ```git clone https://github.com/denilsonfbar/cirCodAn.git```

- Add the cirCodAn and CodAn/bin folder into your PATH:
    - ```export PATH=$PATH:$PWD/cirCodAn/ && export PATH=$PATH:$PWD/cirCodAn/CodAn/bin/```
    - to add cirCodAn permanently to your PATH, add the previous "export" commands into your ```~/.bashrc``` or ```~/.bashprofile```

- Apply "execution permission" to executables if needed:
    - ```chmod +x $PWD/cirCodAn/*.py && chmod +x $PWD/cirCodAn/CodAn/bin/*```

## Requirements

- [Python3](https://www.python.org/), [Biopython](https://biopython.org/wiki/Download), and [Pandas](https://pypi.org/project/pandas/)
    - ```apt-get install python3-biopython python3-pandas```
- [Perl](https://www.perl.org/), [Bioperl](https://bioperl.org/) and [MCE](https://metacpan.org/release/MCE)
    - ```apt-get install bioperl libmce-perl```
- [CodAn](https://github.com/pedronachtigall/CodAn)
- [ToPS](https://tops.sourceforge.net/)

Ensure that all requirements are working properly.

:warning: **Conda environment installation**

If you are not the root user and takes advantage of [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) (i.e., one of the famous environment management system), you can follow the steps below:


- Create the conda environment with all dependencies:
   - ```conda create -n circodan_env -c bioconda codan pandas```
   - ```git clone https://github.com/denilsonfbar/cirCodAn.git```
   - ```export PATH=$PATH:$PWD/cirCodAn/```
       - to add it permanently to your PATH, add the command ```export PATH=$PATH:$PWD/cirCodAn/``` into your ```~/.bashrc``` or ```~/.bashprofile```
   - ```chmod +x $PWD/cirCodAn/*.py```

- Activate the environment before use:
   - ```conda activate circodan_env```

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
  -m folder, --model=folder
                        Optional - path to model folder
                        [default="models/VERT_circ"]
```

Basic usage:
```
circodan.py -f circRNA_seqs.fa
```

### Testing the installation

To test if cirCodAn is properly working, just run the testing set:
```
circodan.py -f example/circRNA_seqs.fa
```

**Expected example output**
```
2023-03-15 17:07:25 -> started cirCodAn v1.0
2023-03-15 17:07:26 -> prediction finished
Number of input sequences -> 138
Number of predicted CDSs  -> 104
GTF file with prediction annotation -> example/cirCodAn_output/CDS_predicted.gtf
Predicted CDS seqs FASTA file -> example/cirCodAn_output/CDS_predicted_seqs.fa
Predicted peptides FASTA file -> example/cirCodAn_output/CDS_predicted_seqs_aa.fa
```

### Outputs
```
cirCodAn_output/
├── CDS_predicted_seqs.fa
├── CDS_predicted_seqs_aa.fa
└── CDS_predicted.gtf
```

## This repository

- CodAn/: original CodAn files
- datasets/: datasets used in cirCodAn development and evaluation
- development/: raw data and scripts of cirCodAn development and evaluation
- example/: circRNA example file
- literature_search/: bibliometric analysis files
- models/: ToPS GHMM models for circRNA


## Citation

If you use or discuss cirCodAn, please cite the following:

(cirCodAn book chapter under peer review)

Nachtigall, P. G., Kashiwabara, A. Y., & Durham, A. M. (2021). CodAn: Predictive models for precise identification of coding regions in eukaryotic transcripts. Briefings in Bioinformatics, 22(3), 1–11. https://doi.org/10.1093/bib/bbaa045

## License

[GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)


## Contact

To report bugs, to ask for help and to give any feedback, please contact denilsonfbar@gmail.com
