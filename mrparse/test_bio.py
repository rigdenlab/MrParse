#!/usr/bin/env ccp4-python
import sys
from Bio import AlignIO
from Bio.SeqFeature import SeqFeature, FeatureLocation


from Bio import SearchIO
align_file = "/opt/mrbump/tests/phmmer.log"
io = SearchIO.read(align_file, 'hmmer3-text')

h1 = io['3f85_A']
len(h1) == 3 # 3 domains

#h1.hsps or h1.fragments
hs = h1.hsps[0] # first domain



sys.exit()
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
