import argparse
from rdkit import Chem
import os
import pandas as pd
from tqdm import tqdm
import oddt
import numpy as np
import joblib
from oddt.fingerprints import PLEC
from mordred import Calculator, descriptors
import mordred
from Bio.PDB import PDBParser
import Bio.PDB
import math
import openbabel.pybel as pybel
import warnings
warnings.filterwarnings("ignore")

def convert_mol2_to_pdb(input_mol2, output_pdb):
    mol = next(pybel.readfile("mol2", input_mol2))
    mol.write("pdb", output_pdb, overwrite=True)

def plec_calc(protein_path, ligands):
    print('Calculating PLEC ...')
    plec_list = []
    for i in range(0,16384):
        plec_list.append(i)
    df_plec = pd.DataFrame(None, columns=plec_list)
    protein = next(oddt.toolkit.readfile('mol2', protein_path))
    protein.protein = True
    for mol in tqdm(ligands):
        ligand = oddt.toolkits.rdk.Molecule(mol)
        plec = PLEC(ligand=ligand, protein=protein, sparse=False)
        df_plec.loc[len(df_plec)] = pd.Series(plec)
    return df_plec

def desc_calc(mols):
    print('Calculating 3D descriptors ...')
    calc = Calculator([mordred.CPSA.PNSA, mordred.CPSA.DPSA, mordred.CPSA.FPSA, 
                   mordred.CPSA.RPCS, mordred.GeometricalIndex.Diameter3D,
                   mordred.GeometricalIndex.GeometricalShapeIndex, mordred.MomentOfInertia.MomentOfInertia])

    df_desc = calc.pandas(mols)

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

    for col in df_desc.columns:
        df_desc[f'{col}:nha'] = df_desc[col] / df_desc['NumHeavyAtoms']
        df_desc[f'{col}:nrb'] = df_desc[col] / df_desc['NumRotatableBonds']

    df_desc.drop(['NumHeavyAtoms:nha', 'NumHeavyAtoms:nrb', 'NumRotatableBonds:nha', 'NumRotatableBonds:nrb', 'NumHeavyAtoms', 'NumRotatableBonds'], axis=1, inplace=True)

    df_desc = df_desc[['RPCS', 'MOMI-Z', 'PNSA2:nrb', 'DPSA4:nha', 'DPSA4:nrb', 'FPSA1:nrb', 'RPCS:nha', 'GeomDiameter:nha', 'GeomShapeIndex:nha']]

    return df_desc

def ram_calc(ligand_path, phi_kde, psi_kde):
    print('Calculating Ramahandran index ...')
    df_ram = pd.DataFrame(None, columns=['Region 4', 'phi_prob', 'psi_prob'])
    lengths = []
    for i, pose in tqdm(enumerate(os.listdir(ligand_path))):
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

        df_ram.loc[i, 'Region 4'] = reg4/residues
        df_ram.loc[i, 'phi_prob'] = phi_prob
        df_ram.loc[i, 'psi_prob'] = psi_prob
        lengths.append(len(phi_psi_list))

    if len(df_ram) != len(lengths):
        raise Exception('Number of ligands in ramachandran index dataframe is different from number of calculated residues.')
    
    return df_ram, lengths

def predict_rmsd(protein, ligands_path, rescore, out_features = False):

    # read ligands with rdkit
    ligands = [Chem.MolFromMol2File(os.path.join(ligands_path, mol), sanitize=False, removeHs=False) for mol in os.listdir(ligands_path)]
    for mol in ligands:
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol,Chem.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|Chem.SanitizeFlags.SANITIZE_SYMMRINGS,catchErrors=True)
    # Calculate the PLEC fingerprints
    df_plec = plec_calc(protein, ligands)
    # Select PLEC based on Variance Threshold
    plec_todrop = np.load(os.path.join(os.path.join(os.getcwd(), 'objects'), 'plec_vt_todrop.npy'), allow_pickle=True).astype(int)
    X_plec_01 = df_plec.drop(plec_todrop, axis=1).to_numpy()
    # Select PLEC based on Random Forest selector
    plec_rf_selector = joblib.load(os.path.join(os.path.join(os.getcwd(), 'objects'), 'plec_selector.joblib'))
    X_plec = plec_rf_selector.transform(X_plec_01)
    # Final selection
    X_plec = X_plec[:, [6,10,15,16,17,21,30,34,35,39,40,57,65,67,71,78,79,86,87,95,100,107,110,114,115,117,128,131,132]]

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
    
    df_resc.sort_values('Name', inplace=True)
    df_resc_sel = df_resc[['APBS_Ligand', 'APBS_Binding', 'CHARMM', 'ElectDD', 'MlpInS3', 'PLANTS_CHEMPLP_NORM_WEIGHT', 'PLANTS_CHEMPLP_NORM_CRT_WEIGHT', 'RPScore_LigContacts', 'XS_HPScore', 'XS_Average', 'XS_Binding']]

    # Get pep length dummy
    pep_dummies = []
    for length in pep_lengths:
        du_encode = [0] * 10
        du_encode[length-1] = 1
        pep_dummies.append(du_encode)


    df_dummy = pd.DataFrame(pep_dummies, columns=['pep length_1', 'pep length_2', 'pep length_3', 'pep length_4', 'pep length_5',
                                                  'pep length_6', 'pep length_7', 'pep length_8', 'pep length_9', 'pep length_10'])
    df_dummy = df_dummy[['pep length_3', 'pep length_5', 'pep length_8', 'pep length_9', 'pep length_10']]

    # concatenate all features
    df = pd.concat((df_dummy, df_ram, df_desc, pd.DataFrame(X_plec), df_resc_sel), axis=1)

    if out_features:
        pd.concat((df_resc['Name'], df), axis=1).to_csv(r'PepScorerRMSD_features.csv', index=False)

    print('Predicting RMSD ...')
    # Predict
    model = joblib.load(os.path.join(os.path.join(os.getcwd(), 'objects'), 'PepScorerRMSD.joblib'))
    preds = model.predict(df.to_numpy())

    # Output the results
    pd.concat((df_resc['Name'], pd.Series(preds).T), axis=1).to_csv(r'PepScorerRMSD_output.csv', index=False)

    print('Done!')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Predicts the RMSD of a peptide binding pose with the model PepScorer::RMSD")
    
    parser.add_argument("-p", "--protein", type=str, required=True, help="Path of the protein mol2 file.")
    parser.add_argument("-l", "--ligand", type=str, required=True, help="Path of the directory with mol2 file/s of the ligand/s.")
    parser.add_argument("-r", "--rescore", type=str, required=True, help="Path of CSV or XLSX file with Rescore+ features.")
    parser.add_argument("-f", "--features", type=bool, help='Save or not calculated features CSV file. Default is "False", write "True" to save.', default=False)
    args = parser.parse_args()
    
    if args.features:
        predict_rmsd(args.protein, args.ligand, args.rescore, out_features = True)
    else:
        predict_rmsd(args.protein, args.ligand, args.rescore)