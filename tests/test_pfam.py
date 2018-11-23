'''
Created on 26 Oct 2018

@author: jmht

https://omictools.com/protein-coiled-coil-domain-prediction-category

paircoil2 MIT
deepcoil **


'''
import json
import os
import set_mrparse_path
from mrparse.mr_classify import MrClassifier
from mrparse.mr_jpred import secondary_structure_prediction
from mrparse.mr_annotation import get_annotation_chunks

# seqin = "/opt/MrParse/data/Q13586.fasta"
# topcons_dir = "/opt/MrParse/data/Q13586/topcons"
# classifier = MrClassifier()
# class_data = classifier.get_prediction(seqin, topcons_dir=topcons_dir)

jpred_dir = "/opt/MrParse/data/Q13586/jpred"
pfam_data = secondary_structure_prediction(jpred_dir)

print(get_annotation_chunks(pfam_data))

