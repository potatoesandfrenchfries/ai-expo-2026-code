import pandas as pd
import numpy as np
import torch

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams
from rdkit.Chem import rdFingerprintGenerator


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
