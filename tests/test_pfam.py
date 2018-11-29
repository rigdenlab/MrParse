'''
Created on 26 Oct 2018

@author: jmht

https://omictools.com/protein-coiled-coil-domain-prediction-category

paircoil2 MIT
deepcoil **


'''
import set_mrparse_path
from mrparse.mr_classify import get_annotation
from mrparse.mr_annotation import get_annotation_chunks
from mrparse.mr_pfam import pfam_classification_dict



def test_basic():
    seqin = "/opt/MrParse/data/Q13586.fasta"
    topcons_dir = "/opt/MrParse/data/Q13586/topcons"
    annotation = get_annotation(seqin, topcons_dir=topcons_dir)
    
    # jpred_dir = "/opt/MrParse/data/Q13586/jpred"
    # pfam_data = secondary_structure_prediction(jpred_dir)
    
    chunks = get_annotation_chunks(annotation)
    print("CHUNKS ",chunks)
    
    print(pfam_classification_dict(chunks, len(annotation)))

