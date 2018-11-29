'''
Created on 26 Oct 2018

@author: jmht
'''
import set_mrparse_path
from mrparse.mr_classify import get_annotation
from mrparse.mr_annotation import get_annotation_chunks
from mrparse.mr_jpred import secondary_structure_prediction
from mrparse.mr_pfam import pfam_classification_dict



def test_basic():
    seqin = "/opt/MrParse/data/Q13586.fasta"
    topcons_dir = "/opt/MrParse/data/Q13586/topcons"
    annotation = get_annotation(seqin, topcons_dir=topcons_dir)
    annotation_chunks = get_annotation_chunks(annotation)
    class_d = pfam_classification_dict(annotation_chunks, len(annotation))
    assert class_d['length'] == 685    
    assert len(class_d['regions']) == 7
    
    jpred_dir = "/opt/MrParse/data/Q13586/jpred"
    sspred = secondary_structure_prediction(jpred_dir)
    sspred_chunks = get_annotation_chunks(sspred)
    sspred_d = pfam_classification_dict(sspred_chunks, len(sspred))
     
    assert sspred_d['length'] == 685    
    assert len(sspred_d['regions']) == 23
    

