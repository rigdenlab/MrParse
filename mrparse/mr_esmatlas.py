"""
Created on 04 Nov 2022

@author: rmk65 & hlasimpk
"""

# Simple function to download a predicted model from the ESMAtlas 

import os,sys

def get_esm_prediction(seqin, debug=False):
    """ Get a prediction from the ESM Atlas """

    sq=open(seqin, "r")
    sqlines=sq.readlines()
    sq.close()
    
    sequence=""
    for line in sqlines:
        if ">" not in line:
            sequence+=line.strip()

    if len(sequence) > 350:
        print("Sequence must be less than 350 residues")
        sys.exit()
    
    outpdb=os.path.join(os.getcwd(), "esm_prediction.pdb")
    
    import requests
    url = 'https://api.esmatlas.com/foldSequence/v1/pdb/'
    esmpdb=requests.post(url, data=sequence, timeout=100.0)

    # Recalcuate the pLDDT column (*100)
    plddt_pdb=""
    for line in esmpdb.text:
        if "ATOM" in line[0:4]:
            plddt_pdb+=("%s %.2lf %s\n" % (line.strip()[:60], float(line.strip()[61:66])*100.0, line.strip()[67:]))
        else:
            plddt_pdb+=line
    
    # Write out the corrected PDB
    op=open(outpdb,"w")
    op.write(plddt_pdb)
    op.close()

    return outpdb

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--seqin', help='input sequence file (fasta format, single sequence)')
parser.add_argument('-d', '--debug', action='store_true', help="debug output")
args=parser.parse_args()

esm_model=get_esm_prediction(args.seqin, debug=args.debug)
