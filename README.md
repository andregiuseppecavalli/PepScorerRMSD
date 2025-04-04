# PepScorer::RMSD

In this section are included the information to install and use **PepScorer::RMSD**.

## Installation

Creating a new python environment, with `python >= 3.12`, is recommended before running the requirements installation.
After downloading and extracting the folder, install all the requirements:

```bash
pip install -r requirements.txt
```

## Testing

Once installed, run the following command to test the installation:

```bash
python predict.py -p test/1BE9_protein.mol2 -l test/poses/ -r test/rescore.xlsx
```

## Usage

**PepScorer::RMSD** works with MOL2 files, each pose should be represented by a single MOL2 file.
Run `python predict.py -h` to see a list of all the available tags:

```bash
usage: predict.py [-h] -p PROTEIN -l LIGAND -r RESCORE [-f FEATURES]

Predicts the RMSD of a peptide binding pose with the model PepScorer::RMSD

options:
  -h, --help            show this help message and exit
  -p PROTEIN, --protein PROTEIN
                        Path of the protein mol2 file.
  -l LIGAND, --ligand LIGAND
                        Path of the directory with mol2 file/s of the ligand/s.
  -r RESCORE, --rescore RESCORE
                        Path of CSV or XLSX file with Rescore+ features.
  -f FEATURES, --features FEATURES
                        Save or not calculated features CSV file. Default is "False", write "True" to save.
```

The `-r` tag is used to pass Rescore+ features calculated with **VegaZZ Rescore+**.
The scores to calculate are: `CHARMM, APBS, ElectDD, RPSCORE, MLPInS3, ChemPlp, XScore`. The order of the columns in your feature CSV/XLSX file should be the same as in the `test/rescore.xlsx` file.

**VegaZZ** natively does not include `XScore` and `ChemPlp`. `XScore` can be downloaded freely at the following link: http://www.sioc-ccbg.ac.cn/?p=42&software=xscore. 
PLANTS can be obtained from the following GitHub repository: https://github.com/purnawanpp/plants to calculate `ChemPlp` scores.