# PepScorer::RMSD

This repository includes the information to install and use **PepScorer::RMSD**, a machine-learning model for protein-peptide pose selection. 
Read more in the original article https://doi.org/10.3390/ijms27020870.
![abstract_graphics](https://github.com/user-attachments/assets/9b47586a-e212-4b0c-b78e-e4e7b3b1b71e)


## Installation

Creating a new python environment, with `python v. 3.12`, is recommended before running the requirements installation.
First, download and extract the folder in your project directory.

## Setting Up a Python Virtual Environment
### 1. Create a virtual environment
Run the following command in your project directory:

```bash
python3 -m venv .venv
```
This will create a folder called `.venv` containing the virtual environment.

### 2. Activate the virtual environment
On Linux: `source .venv/bin/activate`
On Windows: `.venv\Scripts\activate`

## Install all the requirements:

```bash
pip install -r requirements.txt
```
## Testing

Once installed, run the following command to test the installation:

```bash
python predict.py -l test/poses/ -r test/rescore.xlsx
```

## Usage

**PepScorer::RMSD** works with MOL2 files, each pose should be represented by a single MOL2 file.
Run `python predict.py -h` to see a list of all the available options:

```bash
usage: predict.py [-h] -l LIGAND -r RESCORE [-f FEATURES]

Predicts the RMSD of a peptide binding pose with the model PepScorer::RMSD

options:
  -h, --help            show this help message and exit
  -l LIGAND, --ligand LIGAND
                        Path of the directory with mol2 file/s of the ligand/s.
  -r RESCORE, --rescore RESCORE
                        Path of CSV or XLSX file with Rescore+ features.
  -f FEATURES, --features FEATURES
                        Save or not calculated features CSV file. Default is "False", write "True" to save.
```

The `-r` tag is used to pass Rescore+ features calculated with **VegaZZ Rescore+**.
The scores to calculate are: `CHARMM, APBS, ELECT, ELECTDD, MLPINS, MLPINS2, MLPINS3, MLPINSF, RPSCORE, CHEMPLP, XSCORE`. The order of the columns in your feature CSV/XLSX file should be the same as in the `test/rescore.xlsx` file. Check that the order of the poses in your Rescore+ input file is the same as the one in the poses' directory, and that the calculated scores are 30.

**VegaZZ** natively does not include `XScore` and `ChemPlp`. `XScore` can be downloaded freely at the following link: http://www.sioc-ccbg.ac.cn/?p=42&software=xscore. 
PLANTS can be obtained from the following GitHub repository: https://github.com/purnawanpp/plants to calculate `ChemPlp` scores.

## Other files

The file `train.py` can be used to train a model with your own poses.

The file `train_and_evaluate.ipynb` is a jupiter notebook that can be used to replicate the results obtained on evaluation and external test sets.
