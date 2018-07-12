#!/usr/bin/env ccp4-python
import sys
from Bio import AlignIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SearchIO

sys.path.insert(0, '/opt/mrbump/include/seq_align')

from MRBUMP_phmmer import PHHit, Domains
import simpleSeqID

# Compare 1bik
align_file = "/opt/mrbump/tests/phmmer.log"
io = SearchIO.read(align_file, 'hmmer3-text')


assert len(io) == 31

included = io.hit_filter(lambda x: x.is_included)

assert len(included) == 30

resultsDict = {}
targetSequence = 'FOO' # This needs to be passed in

rank = 0
for i, hit in enumerate(included):
    for j, hsp in enumerate(hit.hsps):
        rank += 1
        ph = PHHit()
        ph.rank = rank
        ph.chainName = hsp.hit_id
        name, chain = hsp.hit_id.split('_')
        ph.pdbName = name
        ph.chainID = chain
        ph.score = hsp.bitscore
        ph.evalue = hsp.evalue
        ph.ndomains = len(hit)
        
        
        start, stop = hsp.hit_range
        ph.alnRange = "{}-{}".format(start + 1, stop)
        start, stop = hsp.query_range
        ph.tarRange = "{}-{}".format(start + 1, stop)
        ph.tarExtent = hsp.query_span - 1
        targetAlignment = "".join(hsp.aln[0].upper()) # assume the first Sequence is always the target
        ph.targetAlignment = targetAlignment
        alignment = "".join(hsp.aln[1].upper()) # assume the first Sequence is always the target
        ph.alignment = alignment

        local, overall = simpleSeqID.simpleSeqID().getPercent(alignment, targetAlignment, targetSequence)
        ph.localSEQID = local
        ph.overallSEQID = overall
        
        key = hit.id + str(hsp.domain_index)
        resultsDict[key] = ph


for k, v in resultsDict.items():
    print k, v
        


# Biopython includes 1toc_R1, which is not in mrbump as below threshold
# Filter on is_included


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
