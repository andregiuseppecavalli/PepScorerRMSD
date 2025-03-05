# PepScorer::RMSD

**PepScorer::RMSD** is a peptide-specific machine-learning scoring function for pose selection after docking. The model predicts the root-mean-squared deviation (RMSD) of a given binding pose, therefore, the lower is the RMSD the better is the pose.

## Installation

Creating a new python environment is recomended before running the requirements installation.
Clone the repository and install dependencies:

```bash
git clone https://github.com/andregiuseppecavalli/PepScorerRMSD
cd PepScorerRMSD
pip install -r requirements.txt
```

## Testing

Once installed all the requirements, run the following command to test the installation:

```bash
python predict.py -p test/1BE9_fixedlig_nolig.mol2 -l test/poses/ -r test/rescore.xlsx
```

## Usage

**PepScorer::RMSD** works with MOL2 files, each pose should be represented by a single MOL2 file.
Run `python predict.py -h` to see a list of all the availabel tags:

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
The scores to calculate are: `CHARMM, APBS, ElectDD, RPSCORE, MLPInS3, ChemPlp, XScore`.