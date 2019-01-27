#!/usr/bin/env ccp4-python
'''
Created on 26 Oct 2018

@author: jmht
'''
import set_mrparse_path
# from mrparse.mr_classify import get_annotation
# from mrparse.mr_annotation import get_annotation_chunks
# from mrparse.mr_jpred import JPred
# from mrparse.mr_pfam import pfam_classification_dict

from mrparse.mr_hit import find_hits
from mrparse.mr_homolog import homologs_from_hits
from mrparse.mr_region import RegionFinder
from mrparse.mr_pfam import add_pfam_dict_to_homologs
from mrparse.mr_sequence import read_fasta



def test_homologs(test_data):
    seqin = test_data.x2uvoA_fasta
    hits = find_hits(seqin)
    regions = RegionFinder().find_regions_from_hits(hits)
    homologs = homologs_from_hits(hits)
    seqlen = len(read_fasta(seqin))
    add_pfam_dict_to_homologs(regions, seqlen)
    
    print(homologs.values()[0])
    print(homologs.values()[0]._pfam_json)
    

# def test_local():
#     seqin = "/opt/MrParse/data/Q13586.fasta"
#     topcons_dir = "/opt/MrParse/data/Q13586/topcons"
#     annotation = get_annotation(seqin, topcons_dir=topcons_dir)
#     annotation_chunks = get_annotation_chunks(annotation)
#     class_d = pfam_classification_dict(annotation_chunks, len(annotation))
#     assert class_d['length'] == 685    
#     assert len(class_d['regions']) == 7
#     
#     jpred_rundir = "/opt/MrParse/data/Q13586/jpred"
#     sspred = JPred().secondary_structure_prediction(jpred_rundir=jpred_rundir)
#     sspred_chunks = get_annotation_chunks(sspred)
#     sspred_d = pfam_classification_dict(sspred_chunks, len(sspred))
#      
#     assert sspred_d['length'] == 685    
#     assert len(sspred_d['regions']) == 23
# 
# 
# def test_remote():
#     seqin = "/opt/MrParse/data/Q13586.fasta"
#     annotation = get_annotation(seqin)
#     annotation_chunks = get_annotation_chunks(annotation)
#     class_d = pfam_classification_dict(annotation_chunks, len(annotation))
#     assert class_d['length'] == 685    
#     assert len(class_d['regions']) == 7
#      
#     sspred = JPred().secondary_structure_prediction(seqin=seqin)
#     sspred_chunks = get_annotation_chunks(sspred)
#     sspred_d = pfam_classification_dict(sspred_chunks, len(sspred))
#       
#     assert sspred_d['length'] == 685    
#     assert len(sspred_d['regions']) == 23
#     
#     
if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])
