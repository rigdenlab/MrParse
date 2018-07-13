#!/usr/bin/env ccp4-python
import sys
from collections import OrderedDict
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

targetSequence = 'AVLPQEEEGSGGGQLVTEVTKKEDSCQLGYSAGPCMGMTSRYFYNGTSMACETFQYGGCMGNGNNFVTEKECLQTCRTVAACNLPIVRGPCRAFIQLWAFDAVKGKCVLFPYGGCQGNGNKFYSEKECREYCGVPGDGDEELLRFSN'


resultsDict = OrderedDict()
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
        
        name = hit.id + str(hsp.domain_index)
        # jmht - need to add this
        ph.name = name
        resultsDict[name] = ph


def calculate_domains(resultsDict):
    """Figure out the domains for the target that have been matched"""
    
    def update_or_create_domain(hit, targetDomainDict):
        for domain in targetDomainDict.values():
            if hit_within_domain(hit, domain):
                return update_domain(hit, domain)
        # Domain not found so add new
        add_new_domain(hit, targetDomainDict)
        return

    def hit_within_domain(hit, domain, extentTolerance=50, midpointTolerance=20):
        if hit.tarExtent >= domain.extent - extentTolerance and \
            hit.tarExtent <= domain.extent + extentTolerance and \
            hit.tarMidpoint >= domain.midpoint - midpointTolerance and \
            hit.tarMidpoint <= domain.midpoint + midpointTolerance:
            return False
        return True
    
    def add_new_domain(hit, targetDomainDict):
        if not targetDomainDict:
            domCount = 1
        else:
            domCount = max(targetDomainDict.keys()) + 1
        targetDomainDict[domCount] = Domains()
        targetDomainDict[domCount].ID = domCount
        targetDomainDict[domCount].midpoint = hit.tarMidpoint
        targetDomainDict[domCount].extent = hit.tarExtent
        targetDomainDict[domCount].matches.append(hit.name)
        targetDomainDict[domCount].ranges.append(hit.tarRange)
        return targetDomainDict
    
    def update_domain(hit, domain):
        domain.matches.append(hit.name)
        domain.ranges.append(hit.tarRange)
    
    # Set up first domain
    targetDomainDict = {}
    for hit in resultsDict.values():
        update_or_create_domain(hit, targetDomainDict)
    return targetDomainDict


targetDomainDict = calculate_domains(resultsDict)

print "GOT ",len(targetDomainDict)
# for k in sorted(resultsDict.keys()):
#     print k, resultsDict[k]
for k in sorted(targetDomainDict.keys()):
    print k, targetDomainDict[k]
sys.exit()

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


# Check domains
assert(len(targetDomainDict) == 4)
key = 3
dd = targetDomainDict[key]
assert(dd.ID == key)
assert(dd.extent == 52)
assert(dd.matches == ['4bd9_B2', '2ody_F2', '4isl_B1', '1tfx_D1', '3uir_D1', '4u32_X1', '1bz5_D1', '1zr0_B1', '2m01_A1', '2m99_A1', '4bqd_A1', 
                      '1jc6_A1', '1irh_A1', '2kcr_A1', '1dtk_A1', '1den_A1', '3wny_I1', '1y62_C1', '4nty_B1'])
assert(dd.midpoint == 50.0)
assert(dd.ranges == ['24-76', '24-78', '22-79', '22-77', '23-77', '24-76', '24-76', '25-76', '25-77', '25-76', '26-76', '25-77', '29-77', '24-78', '26-76', '26-76', '31-76', '25-76', '35-76'])


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
