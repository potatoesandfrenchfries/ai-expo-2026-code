import pandas as pd
import numpy as np
import torch
import json 

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, FilterCatalog, rdFingerprintGenerator
from rdkit.Chem.FilterCatalog import FilterCatalogParams


df= pd.read_csv('src/Chemistry Formalisms/Data/250k_rndm_zinc_drugs_clean_3.csv')

print(f"Loaded {len(df)} molecules from the dataset.")

# print(df.head())

gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

#Sampling
df=df.sample(n=10000, random_state=42).reset_index(drop=True)

valid_rows = []
for i, row in df.iterrows():
    mol = Chem.MolFromSmiles(row['smiles'])
    if mol is not None:
        valid_rows.append(i)
df = df.iloc[valid_rows].reset_index(drop=True)

print(f"After filtering invalid SMILES, {len(df)} molecules remain.")

# features

fingerprints = []

for smi in df['smiles']:
    mol = Chem.MolFromSmiles(smi)
    fp = gen.GetFingerprintAsNumPy(mol)
    
    fingerprints.append(fp)

X = torch.FloatTensor(np.array(fingerprints))

print(f"Feature matrix shape: {X.shape}")

# labelling

y=torch.FloatTensor(df['qed'].values)

'''
print(f"Labels shape: {y.shape}")
print(f"min: {y.min()}, max: {y.max()}, mean: {y.mean()}, std: {y.std()}")
'''

# substructure flags

# PAINS filter : pan-assay interference compounds

pains_params = FilterCatalogParams()
pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)

pains_catalog = FilterCatalog.FilterCatalog(pains_params)

# Brenk filter: alerts for reactive/toxic groups (implementing only acid halides, peroxides, and anhydrides for simplicity) (maybe michael acceptors if rocq cooperates)

brenk_params = FilterCatalogParams()
brenk_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
brenk_catalog = FilterCatalog.FilterCatalog(brenk_params)

toxic_flags = []
for smi in df['smiles']:
    mol = Chem.MolFromSmiles(smi)
    
    pains_flag = pains_catalog.HasMatch(mol)
    brenk_flag = brenk_catalog.HasMatch(mol)
    
    toxic_flags.append(1 if (pains_flag or brenk_flag) else 0)
    
toxic_flags = torch.IntTensor(toxic_flags)

print(f"Toxic flags: {toxic_flags.sum().item()} out of {len(toxic_flags)} ({toxic_flags.float().mean()*100:.2f}%) ")


# ── Atom-level data for Rocq record generation ────────────────────────────────

# For each molecule, extract the per-atom information
# that verifier.py needs to auto-generate the Rocq source file

atom_data = []

for smi in df['smiles']:
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol) # Add explicit hydrogens to get accurate atom counts and properties

    atoms = []
    for atom in mol.GetAtoms():
        atoms.append({
            'symbol':      atom.GetSymbol(),
            'bonds':       atom.GetTotalDegree(),
            'is_aromatic': atom.GetIsAromatic()
        })

    mw = Descriptors.ExactMolWt(mol)

    atom_data.append({
        'smiles':    smi,
        'atoms':     atoms,
        'mw':        round(mw, 4),
        'toxic_flag': toxic_flags[len(atom_data)]
    })

print(f"Atom data collected for {len(atom_data)} molecules")
print(f"Example entry:\n  {atom_data[0]['smiles']}")
print(f"  MW: {atom_data[0]['mw']}")
print(f"  Atoms: {len(atom_data[0]['atoms'])}")
print(f"  Toxic: {atom_data[0]['toxic_flag']}")

# Files
torch.save(X,     'data/X.pt')
torch.save(y,     'data/y.pt')
torch.save(toxic_flags, 'data/toxic.pt')

with open('data/atom_data.json', 'w') as f:
    json.dump(atom_data, f, indent=2)

df[['smiles', 'qed']].to_csv('data/molecules_clean.csv', index=False)

print("\n✓ Done. Files saved:")
print(f"  X.pt             {tuple(X.shape)}  — BNN input features")
print(f"  y.pt             {tuple(y.shape)}   — QED labels")
print(f"  toxic.pt         {tuple(toxic_flags.shape)}   — toxic flags")
print(f"  atom_data.json   {len(atom_data)} molecules — Rocq bridge input")
print(f"  molecules_clean.csv")