# Definition of file name
filename = 'aminoacids.txt'

# Libraries needed
import sys
import pubchempy as pcp
import pandas as pd

# We use pandas to read the file we defined previously, which does not have a header.
# We're adding a header called 'Compound' to the DataFrame for clarity.
df = pd.read_csv(filename,header=None,names=['Compound'])

# We create a list to add all SMILES
SMILES_list = []

# Iterate through each identifier in the 'Compound' column to obtain compound names from PubChem
for ids in df['Compound']:
    compound_name = pcp.get_compounds(ids,'name')
    
    # Once we have compound names, iterate through every name to obtain canonical SMILES from PubChem 
    for name in compound_name:
        smiles = name.canonical_smiles
        # Then, we add canonical SMILES to a list called SMILES_list
        SMILES_list.append(smiles)

# We turn the list into a dictionary to remove all repeated canonical SMILES while maintaining order        
canonical_SMILES = list(dict.fromkeys(SMILES_list))

# We create a DataFrame with canonical SMILES
smiles_df = pd.DataFrame(data=canonical_SMILES)
# Rename the column
smiles_df.columns = ['SMILES']
# Insert the DataFrame with compound names into the DataFrame with SMILES
smiles_df.insert(0,'Compound',df['Compound'],True)

# We create a new list called data to add all properties
data = []

# We iterate through every SMILES in smiles_df to get properties of every molecule
for molecule in smiles_df['SMILES']:
    # We add all properties to a list
    prop = pcp.get_properties(['MolecularFormula', 'MolecularWeight','InChI', 'InChIKey', 'IUPACName', 
                                'XLogP', 'ExactMass', 'MonoisotopicMass', 'TPSA', 'Complexity', 'Charge', 
                                'HBondDonorCount', 'HBondAcceptorCount', 'RotatableBondCount', 
                                'HeavyAtomCount', 'IsotopeAtomCount', 'AtomStereoCount', 
                                'DefinedAtomStereoCount', 'UndefinedAtomStereoCount', 'BondStereoCount', 
                                'DefinedBondStereoCount', 'UndefinedBondStereoCount', 'CovalentUnitCount', 
                                'Volume3D', 'XStericQuadrupole3D', 'YStericQuadrupole3D', 
                                'ZStericQuadrupole3D', 'FeatureCount3D', 'FeatureAcceptorCount3D', 
                                'FeatureDonorCount3D', 'FeatureAnionCount3D', 'FeatureCationCount3D', 
                                'FeatureRingCount3D', 'FeatureHydrophobeCount3D', 'ConformerModelRMSD3D', 
                                'EffectiveRotorCount3D', 'ConformerCount3D'],molecule,'smiles')
    # Once we have the properties, we add them to the data
    data.append(prop)

# We define a new empty list called rows
# and a list called columns with data keys
rows = []
columns = data[0][0].keys()

# We iterate through all data list
for i in range(len(data)):
    # And we add the properties of each molecule into the rows list
    rows.append(data[i][0].values())
# We define a new DataFrame called props_df with all properties we obtained before
props_df = pd.DataFrame(data=rows,columns=columns)

# We merge the DataFrame containing SMILES with the DataFrame containing all properties
names_props = pd.concat([smiles_df,props_df],axis=1)

# We create a new document with all information
names_props.to_csv('molecules_with_properties.csv',index=False)