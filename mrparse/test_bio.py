#!/usr/bin/env ccp4-python
#
from Bio import AlignIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

align_file = "phmmerAlignment.log"

alignment = AlignIO.read(align_file, "stockholm")


# Use .letter_annotations for (e.g. resSeq) information
# Use .annotations for random information
# Use .features for more structure information - domains?

rcrd = alignment[0]
rcrd.letter_annotations['resseq'] = list(range(len(rcrd)))
print rcrd.letter_annotations


example_feature1 = SeqFeature(FeatureLocation(5, 18), type="domain", qualifiers={'foo' : 'bar1'})
rcrd.features.append(example_feature1)
example_feature2 = SeqFeature(FeatureLocation(40, 45), type="domain", qualifiers={'foo' : 'bar2'})
rcrd.features.append(example_feature2)

for f in rcrd.features:
    if 10 in f:
        print("GOT FOR FEATURE ",f)

slice1 = rcrd[0:30]
print(slice1)
# Only contains one feature
