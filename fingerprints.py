## Python script
## Calculate morgan fingerprints from SMILES using RDKit

import sys
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from rdkit.Chem import PandasTools


## fingerprint params
radius=3
nBits=1024
ECFP6 = [AllChem.GetMorganFingerprintAsBitVect(x,radius=radius, nBits=nBits) for x in ontox_toxcast_no_na['ROMol']]
ECFP6[0]

## convert to fingerprints
ecfp6_name = [f'Bit_{i}' for i in range(nBits)]
ecfp6_bits = [list(l) for l in ECFP6]
df_morgan = pd.DataFrame(ecfp6_bits, index = ontox_toxcast_no_na.canonical_smiles, columns=ecfp6_name)
df_morgan.index.name = 'canonical_smiles'
df_morgan.reset_index(inplace=True)
df_morgan.to_csv('./inst/data_out/compound_morgan_fp_ontox_toxcast.csv', index=False, sep=',')
df_morgan.head(10)
