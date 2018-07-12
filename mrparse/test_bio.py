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
targetSequence = 'AVLPQEEEGSGGGQLVTEVTKKEDSCQLGYSAGPCMGMTSRYFYNGTSMACETFQYGGCMGNGNNFVTEKECLQTCRTVAACNLPIVRGPCRAFIQLWAFDAVKGKCVLFPYGGCQGNGNKFYSEKECREYCGVPGDGDEELLRFSN'

rank = 0
for i, hit in enumerate(included):
    rank = i + 1
    for j, hsp in enumerate(hit.hsps):
        ph = PHHit()
        ph.rank = rank
        ph.chainName = hsp.hit_id
        name, chain = hsp.hit_id.split('_')
        ph.pdbName = name
        ph.chainID = chain
        #ph.score = hsp.bitscore
        ph.score = hit.bitscore
        ph.evalue = hsp.evalue # is i-Evalue - possibly evalue_cond in later BioPython
        ph.ndomains = len(hit)
        
        start, stop = hsp.hit_range
        ph.alnRange = "{}-{}".format(start + 1, stop)
        start, stop = hsp.query_range
        start_p1 = start + 1
        ph.tarRange = "{}-{}".format(start_p1, stop)
        ph.tarExtent = hsp.query_span - 1
        ph.tarMidpoint = ((float(stop) - float(start_p1)) / 2.0) + float(start_p1)
        
        targetAlignment = "".join(hsp.aln[0].upper()) # assume the first Sequence is always the target
        ph.targetAlignment = targetAlignment
        alignment = "".join(hsp.aln[1].upper()) # assume the first Sequence is always the target
        ph.alignment = alignment

        local, overall = simpleSeqID.simpleSeqID().getPercent(alignment, targetAlignment, targetSequence)
        ph.localSEQID = local
        ph.overallSEQID = overall
        
        key = hit.id + str(hsp.domain_index)
        resultsDict[key] = ph


# for k in sorted(resultsDict.keys()):
#     print k, resultsDict[k]
        

k = '3f85_A2'
assert(k in resultsDict.keys())
assert(resultsDict[k].ndomains == 3)
assert(resultsDict[k].rank == 29)
assert(resultsDict[k].chainID == 'A')
assert(resultsDict[k].chainName =='3f85_A' )
#assert(resultsDict[k].evalue == 0.00016)
assert(resultsDict[k].pdbName == '3f85')
assert(resultsDict[k].score == 22.7)
assert(resultsDict[k].alnRange == '221-252')
assert(resultsDict[k].tarExtent == 31)
assert(resultsDict[k].tarRange == '10-41')
assert(abs(resultsDict[k].tarMidpoint - 25.5) <= 0.0001)
assert(resultsDict[k].alignment == 'SWGGQVLSTTAKEFEAAALGYSASGVNGVSSS')
assert(resultsDict[k].targetAlignment == 'SGGGQLVTEVTKKEDSCQLGYSAGPCMGMTSR')
assert(resultsDict[k].localSEQID == 37)
assert(resultsDict[k].overallSEQID == 8)


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
