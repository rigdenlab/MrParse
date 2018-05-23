#!/usr/bin/env ccp4-python

import pandas as pd

import sys

# http://download.cathdb.info/cath/releases/latest-release/cath-classification-data/README-cath-list-file-format.txt

df = pd.read_csv('tabularResults.csv')
df['PDB ID'] = df['PDB ID'].str.lower()
pdbs = df['PDB ID'].unique().tolist()

tuples = list(zip(df['PDB ID'].values, df['Chain ID'].values))
index = pd.MultiIndex.from_tuples(tuples, names=['pdb', 'chain'])
df.set_index(index, drop=True, inplace=True)

res = {}
with open('cath-domain-boundaries.txt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split()
        ndoms = int(fields[1][1:3])
        if ndoms < 3:
            continue
        pid = fields[0]
        pdb = pid[:4]
        if pdb in pdbs:
            chain = pid[4]
            idx = (pdb, chain)
            if idx in df.index:
                print pdb, chain, df.loc[idx, 'Chain Length']
            # look up chain Length



# Chain Length

# 2CWG - v. old - => use 2UVO

