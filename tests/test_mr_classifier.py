'''
Created on 17 Nov 2018

@author: jmht
'''
import set_mrparse_path
import conftest

from mrparse.mr_classify import MrClassifier

seqin = '/opt/MrParse/data/Q13586.fasta'
classifier = MrClassifier()
classifier.get_predictions(seqin, topcons_dir="/opt/MrParse/data/Q13586/topcons", jpred_dir="/opt/MrParse/data/Q13586/jpred")
print(classifier.classification)
print(classifier.sspred)
print(classifier.pfam_data())
#class_data = classifier.pfam_data()