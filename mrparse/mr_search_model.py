'''
Created on 18 Oct 2018

@author: jmht
'''
import os
from collections import OrderedDict


# from Bio import AlignIO
# from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SearchIO

from ample.util.sequence_util import Sequence

from pyjob import cexec
from pyjob.script import EXE_EXT    

from mrbump.seq_align.simpleSeqID import simpleSeqID
from mrbump.seq_align.MRBUMP_phmmer import PHHit, Domains


class DomainFinder(object):
    def __init__(self):
        pass
     
    def find_hits(self, seqin):
        
        self.run_phmmer(seqin)
        targetSequence = Sequence(fasta=seqin).sequence()
        
        io = SearchIO.read(self.logfile, 'hmmer3-text')
        included = io.hit_filter(lambda x: x.is_included)
        
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
        
                local, overall = simpleSeqID().getPercent(alignment, targetAlignment, targetSequence)
                ph.localSEQID = local
                ph.overallSEQID = overall
                
                name = hit.id + str(hsp.domain_index)
                # jmht - need to add this
                ph.name = name
                resultsDict[name] = ph

        return resultsDict

    def run_phmmer(self, seqin):
        self.logfile = "phmmer.log"
        self.alnfile = "phmmerAlignment.log"
        self.phmmerTblout = "phmmerTblout.log"
        self.phmmerDomTblout = "phmmerDomTblout.log"
        self.phmmerEXE = os.path.join(os.environ["CCP4"], "libexec", "phmmer")
        self.seqdb = os.path.join(os.environ["CCP4"], "share", "mrbump", "data", "pdb70.txt")
        # Set the command line
        cmd = [self.phmmerEXE + EXE_EXT,
               '--notextw',
               '--tblout', self.phmmerTblout,
               '--domtblout', self.phmmerDomTblout,
               '-A', self.alnfile,
               seqin, self.seqdb]
        
        stdout = cexec(cmd)

        with open(self.logfile, 'w') as f_out:
            f_out.write(stdout)

    def find_domains(self, seqin):
        hitDict = self.find_hits(seqin)
        targetDomainDict = self._find_domains(hitDict)
        return targetDomainDict

    def _find_domains(self, resultsDict):
        """Figure out the domains for the target that have been matched"""
        # Set up first domain
        targetDomainDict = {}
        for hit in resultsDict.values():
            self.update_or_create_domain(hit, targetDomainDict)
        return targetDomainDict
    
    def update_or_create_domain(self, hit, targetDomainDict):
        for domain in targetDomainDict.values():
            if self.hit_within_domain(hit, domain):
                return self.update_domain(hit, domain)
        # Domain not found so add new
        self.add_new_domain(hit, targetDomainDict)
        return

    def hit_within_domain(self, hit, domain, extentTolerance=50, midpointTolerance=20):
        if hit.tarExtent >= domain.extent - extentTolerance and \
            hit.tarExtent <= domain.extent + extentTolerance and \
            hit.tarMidpoint >= domain.midpoint - midpointTolerance and \
            hit.tarMidpoint <= domain.midpoint + midpointTolerance:
            return True
        return False
    
    def add_new_domain(self, hit, targetDomainDict):
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
    
    def update_domain(self, hit, domain):
        domain.matches.append(hit.name)
        domain.ranges.append(hit.tarRange)
        return
    

class SearchModelFinder(object):
    def __init__(self, seqin):
        self.seqin = seqin

