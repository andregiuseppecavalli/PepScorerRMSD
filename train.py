import argparse
from rdkit import Chem
import os
import pandas as pd
from tqdm import tqdm
import numpy as np
import joblib
from mordred import Calculator, descriptors
import mordred
from Bio.PDB import PDBParser
import Bio.PDB
import math
import openbabel.pybel as pybel
import warnings
from sklearn.ensemble import HistGradientBoostingRegressor
warnings.filterwarnings("ignore")

def convert_mol2_to_pdb(input_mol2, output_pdb):
    mol = next(pybel.readfile("mol2", input_mol2))
    mol.write("pdb", output_pdb, overwrite=True)

def desc_calc(mols):
    print('Calculating 3D descriptors ...')
    # descriptors that will be calculated

    calc = Calculator([mordred.CPSA.PNSA, mordred.CPSA.DPSA, mordred.CPSA.PPSA, mordred.CPSA.FNSA, mordred.CPSA.FPSA, mordred.CPSA.WNSA, mordred.CPSA.WPSA, mordred.CPSA.RNCS, 
                   mordred.CPSA.RPCS, mordred.CPSA.TASA, mordred.CPSA.TPSA, mordred.CPSA.RASA, mordred.CPSA.RPSA, mordred.GeometricalIndex.Diameter3D, mordred.GeometricalIndex.Radius3D, 
                   mordred.GeometricalIndex.GeometricalShapeIndex, mordred.GeometricalIndex.PetitjeanIndex3D, mordred.MomentOfInertia.MomentOfInertia, mordred.PBF.PBF])
    df_desc = calc.pandas(mols)
    # calculate number of heavy atoms and number of rotatable bonds for the descriptors' normalization
    rdk_desc = {
        'NumHeavyAtoms': [],
        'NumRotatableBonds': []
    }
    for mol in mols:
        nha = Chem.rdMolDescriptors.CalcNumHeavyAtoms(mol)
        nrb = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
        rdk_desc['NumHeavyAtoms'].append(nha)
        rdk_desc['NumRotatableBonds'].append(nrb)

    df_desc = pd.concat((df_desc, pd.DataFrame(rdk_desc)), axis=1)

    for col in df_desc.columns: # descriptors normalization
        df_desc[f'{col}:nha'] = df_desc[col] / df_desc['NumHeavyAtoms']
        df_desc[f'{col}:nrb'] = df_desc[col] / df_desc['NumRotatableBonds']

    df_desc.drop(['NumHeavyAtoms:nha', 'NumHeavyAtoms:nrb', 'NumRotatableBonds:nha', 'NumRotatableBonds:nrb', 'NumHeavyAtoms', 'NumRotatableBonds'], axis=1, inplace=True)
    desc_todrop = np.load(os.path.join(os.path.join(os.getcwd(), 'objects'), '3ddesc_todrop.npy'), allow_pickle=True)

    return df_desc.drop(desc_todrop, axis=1)

def ram_calc(ligand_path, phi_kde, psi_kde):
    print('Calculating Ramahandran index ...')
    df_ram = pd.DataFrame(None, columns=['Region 1','Region 2','Region 3', 'Region 4', 'phi_mean','psi_mean','phi_prob', 'psi_prob'])
    lengths = []
    for i, pose in tqdm(enumerate(os.listdir(ligand_path)), total=len(os.listdir(ligand_path))):
        pdb = pose.split('.')[0] + '.pdb'
        pose_path = os.path.join(ligand_path, pose)
        pdb_path = os.path.join(ligand_path, pdb)
        convert_mol2_to_pdb(pose_path, pdb_path)
        reg1 = 0
        reg2 = 0
        reg3 = 0
        reg4 = 0
        residues = 0
        phi_list = []
        psi_list = []
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_path)

        for model in structure:
            for chain in model:
                poly = Bio.PDB.Polypeptide.Polypeptide(chain)
                phi_psi_list = poly.get_phi_psi_list()
                for residue in phi_psi_list:
                    if None not in residue:
                        
                        phi = math.degrees(residue[0])
                        psi = math.degrees(residue[1])
                        phi_list.append(phi)
                        psi_list.append(psi)
                        phi_prob = np.mean(phi_kde.evaluate(phi_list))
                        psi_prob = np.mean(psi_kde.evaluate(psi_list))
                        
                        if ((phi > -130 and phi < -50) and (psi > 120 and psi < 180)) or ((phi > -75 and phi < -60) and (psi > -50 and psi < -25)):
                            reg1 += 1
                        elif ((phi > -150 and phi < -45) and (psi > 100 and psi < 180)) or ((phi > -90 and phi < -45) and (psi > -65 and psi < 0)):
                            reg2 += 1
                        elif ((phi > -180 and phi < -30) and (psi > -180 and psi < 180)) or ((phi > 30 and phi < 105) and (psi > -30 and psi < 90)):
                            reg3 += 1
                        else:
                            reg4 += 1

                        residues += 1

        os.remove(pdb_path)
        df_ram.loc[i, 'Region 1'] = reg1/residues
        df_ram.loc[i, 'Region 2'] = reg2/residues
        df_ram.loc[i, 'Region 3'] = reg3/residues
        df_ram.loc[i, 'Region 4'] = reg4/residues
        df_ram.loc[i, 'phi_mean'] = np.array(phi_list).mean()
        df_ram.loc[i, 'psi_mean'] = np.array(psi_list).mean()
        df_ram.loc[i, 'phi_prob'] = phi_prob
        df_ram.loc[i, 'psi_prob'] = psi_prob
        lengths.append(len(phi_psi_list))

    if len(df_ram) != len(lengths):
        raise Exception('Number of ligands in ramachandran index dataframe is different from number of calculated residues.')
    
    return df_ram, lengths

def train_model(ligands_path, rescore, y_path, out_features = False):

    # read ligands with rdkit
    ligands = [Chem.MolFromMol2File(os.path.join(ligands_path, mol), sanitize=False, removeHs=False) for mol in os.listdir(ligands_path)]
    for mol in ligands:
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol,Chem.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|Chem.SanitizeFlags.SANITIZE_SYMMRINGS,catchErrors=True)
        
    # read RMSD file
    if y_path.split('.')[-1] == 'csv':
        y = pd.read_csv(y_path).to_numpy()
    elif y_path.split('.')[-1] == 'xlsx':
        y = pd.read_excel(y_path).to_numpy()
    else:
        raise Exception("RMSD file must be a CSV or XLSX file.")
    
    # Calculate 3D descriptors
    df_desc = desc_calc(ligands)

    # calculate Ramachandran index
    phi_kde = joblib.load(os.path.join(os.path.join(os.getcwd(), 'objects'), 'phi_kde.joblib'))
    psi_kde = joblib.load(os.path.join(os.path.join(os.getcwd(), 'objects'), 'psi_kde.joblib'))
    df_ram, pep_lengths = ram_calc(ligands_path, phi_kde, psi_kde)
   
    # Read Rescore+ CSV file
    if rescore.split('.')[-1] == 'csv':
        df_resc = pd.read_csv(rescore)
    elif rescore.split('.')[-1] == 'xlsx':
        df_resc = pd.read_excel(rescore)
    else:
        raise Exception("Rescore+ file must be a CSV or XLSX file.")
        
    # Get pep length dummy
    pep_dummies = []
    for length in pep_lengths:
        du_encode = [0] * 8
        du_encode[length-3] = 1
        pep_dummies.append(du_encode)


    df_dummy = pd.DataFrame(pep_dummies, columns=['pep length_3', 'pep length_4', 'pep length_5',
                                                  'pep length_6', 'pep length_7', 'pep length_8', 'pep length_9', 'pep length_10'])

    # concatenate all features
    df = pd.concat((df_dummy, df_resc.drop(columns='Name'), df_desc, df_ram), axis=1)
    if df.shape[1] != 92:
        raise Exception(f'Number of features before SFS is {df.shape[1]}, it should be 92!')
    # select features with Sequential Feature Selector
    sfs = joblib.load(os.path.join(os.path.join(os.getcwd(), 'objects'), 'sfs.joblib'))
    X = sfs.transform(df)
    if out_features: # save features file
        pd.concat((df_resc['Name'], df.iloc[:,list(sfs.k_feature_idx_)]), axis=1).to_csv(r'features.csv', index=False)

    print('Predicting RMSD ...')
    # Train
    model = HistGradientBoostingRegressor()
    model.fit(X, y)

    # Save the model
    joblib.dump(model, "my_model.joblib")
    print('Model trained!')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train an Histogram Gradient Boosting regressor with your own poses. The output is 'my_model.joblib' file, saved in the current directory.")
    
    parser.add_argument("-l", "--ligand", type=str, required=True, help="Path of the directory with mol2 file/s of the ligand/s.")
    parser.add_argument("-r", "--rescore", type=str, required=True, help="Path of CSV or XLSX file with Rescore+ features.")
    parser.add_argument("-y", "--rmsd", type=str, required=True, help="Path of CSV or XLSX file with RMSD of the poses. The file should contain one column only with RMSD values and with the header.")
    parser.add_argument("-f", "--features", type=bool, help='Save or not calculated features CSV file. Default is "False", write "True" to save.', default=False)
    args = parser.parse_args()
    
    if args.features:
        train_model(args.ligand, args.rescore, args.rmsd, out_features = True)
    else:
        train_model(args.ligand, args.rescore, args.rmsd)