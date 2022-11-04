"""
Created on 04 Nov 2022

@author: rmk65 & hlasimpk
"""

# Simple function to download a predicted model from the ESMAtlas 

import os,sys

def get_esm_prediction(sequence=None, fasta=None, pdbout=None, debug=False):
    """ Get a prediction from the ESM Atlas """

    # Read the input sequence. Sequence string takes priority if both are provided
    if sequence is not None:
        esm_sequence=sequence.strip()
    elif fasta is not None:
        sq=open(fasta, "r")
        sqlines=sq.readlines()
        sq.close()
    
        esm_sequence=""
        for line in sqlines:
            if ">" not in line:
                esm_sequence+=line.strip()
    else:
        print("A sequence file or a sequence string needs to be provided")
        return None

    if len(esm_sequence) > 350:
        print("Sequence must be less than 350 residues")
        return None
    
    import requests
    url = 'https://api.esmatlas.com/foldSequence/v1/pdb/'
    esmpdb=requests.post(url, data=esm_sequence, timeout=100.0)

    # Recalcuate the pLDDT column (*100)
    plddt_pdb=""
    for line in esmpdb.iter_lines():
        if "ATOM" in line.decode()[0:4]:
            plddt_pdb+=("%s %.2lf %s\n" % (line.decode().strip()[:60], float(line.decode().strip()[61:66])*100.0, line.decode().strip()[67:]))
        else:
            plddt_pdb+=line.decode() + "\n"
    
    # Write out the corrected PDB
    if pdbout is None:
        outpdb=os.path.join(os.getcwd(), "esm_prediction.pdb")
    else:
        outpdb=pdbout
    
    op=open(outpdb,"w")
    op.write(plddt_pdb)
    op.close()

    return outpdb

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--sequence', help="input sequence", default=None)
parser.add_argument('-d', '--debug', action='store_true', help="debug output")
parser.add_argument('-f', '--fasta', help='input fasta format sequence file (single sequence)', default=None)
parser.add_argument('-p', '--pdbout', default=None)
args=parser.parse_args()

esm_model=get_esm_prediction(sequence=args.sequence, fasta=args.fasta, pdbout=args.pdbout, debug=args.debug)
