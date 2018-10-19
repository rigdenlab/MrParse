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

from simbad.util.pdb_util import PdbStructure

from mrbump.seq_align.simpleSeqID import simpleSeqID
from mrbump.seq_align.MRBUMP_phmmer import PHHit, Domains


def find_hits(seqin):
    phmmer_logfile = run_phmmer(seqin)
    targetSequence = Sequence(fasta=seqin).sequence()
    
    io = SearchIO.read(phmmer_logfile, 'hmmer3-text')
    included = io.hit_filter(lambda x: x.is_included)
    
    hitDict = OrderedDict()
    rank = 0
    for i, hit in enumerate(included):
        rank = i + 1
        for hsp in hit.hsps:
#             sequence_identity = float(sum([a==b for a, b in zip(*[str(s.upper().seq) for s in hsp.aln.get_all_seqs()])])) / float(hsp.aln_span) 
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
            
            name = hit.id + "_" + str(hsp.domain_index)
            # jmht - need to add this
            ph.name = name
            hitDict[name] = ph

    return hitDict

def run_phmmer(seqin):
    logfile = "phmmer.log"
    alnfile = "phmmerAlignment.log"
    phmmerTblout = "phmmerTblout.log"
    phmmerDomTblout = "phmmerDomTblout.log"
    phmmerEXE = os.path.join(os.environ["CCP4"], "libexec", "phmmer")
    seqdb = os.path.join(os.environ["CCP4"], "share", "mrbump", "data", "pdb70.txt")
    cmd = [phmmerEXE + EXE_EXT,
           '--notextw',
           '--tblout', phmmerTblout,
           '--domtblout', phmmerDomTblout,
           '-A', alnfile,
           seqin, seqdb]
    stdout = cexec(cmd)
    with open(logfile, 'w') as f_out:
        f_out.write(stdout)
    return logfile 


class DomainFinder(object):
    def __init__(self):
        pass
     
    def find_domains(self, seqin):
        hits = find_hits(seqin)
        targetDomains = self.find_domains_from_hits(hits)
        return targetDomains

    def find_domains_from_hits(self, hits):
        """Figure out the domains for the target that have been matched"""
        # Set up first domain
        targetDomains = []
        for hit in hits.values():
            self.update_or_create_domain(hit, targetDomains)
        return targetDomains
    
    def update_or_create_domain(self, hit, targetDomains):
        for domain in targetDomains:
            if self.hit_within_domain(hit, domain):
                return self.update_domain(hit, domain)
        # Domain not found so add new
        self.add_new_domain(hit, targetDomains)
        return

    def hit_within_domain(self, hit, domain, extentTolerance=50, midpointTolerance=20):
        if hit.tarExtent >= domain.extent - extentTolerance and \
            hit.tarExtent <= domain.extent + extentTolerance and \
            hit.tarMidpoint >= domain.midpoint - midpointTolerance and \
            hit.tarMidpoint <= domain.midpoint + midpointTolerance:
            return True
        return False
    
    def add_new_domain(self, hit, targetDomains):
        domain = Domains()
        domain.ID = len(targetDomains) + 1
        domain.midpoint = hit.tarMidpoint
        domain.extent = hit.tarExtent
        domain.matches.append(hit.name)
        domain.ranges.append(hit.tarRange)
        targetDomains.append(domain)
        return targetDomains
    
    def update_domain(self, hit, domain):
        # Should we update the midpoint and extent of the domain?
        domain.matches.append(hit.name)
        domain.ranges.append(hit.tarRange)
        return
    

def get_homologs(hits, domains):
    pdb_dir = 'pdb_downloads'
    if not os.path.isdir(pdb_dir):
        os.mkdir(pdb_dir)
    homologs = {}
    for domain in domains:
        for match in domain.matches:
            hit = hits[match]
            #pdb_code, chain_id, _ = match.split("_")
            pdb_struct = PdbStructure()
            pdb_struct.from_pdb_code(hit.pdbName)
            pdb_struct.select_chain_by_id(hit.chainID)
            pdb_struct.standardize()
            pdb_name = hit.pdbName + '_' + hit.chainID + '.pdb'
            fpath = os.path.join(pdb_dir, pdb_name)
            homologs[hit.name] = {}
            homologs[hit.name]['pdb'] = fpath
            homologs[hit.name]['mw'] = pdb_struct.molecular_weight
            homologs[hit.name]['seqid'] = hits[hit.name].localSEQID / 100.0
            pdb_struct.save(fpath)
    return homologs

    

class SearchModelFinder(object):
    def __init__(self, seqin):
        self.seqin = seqin

