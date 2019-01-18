'''
Created on 18 Oct 2018

@author: jmht
'''
from collections import OrderedDict
import os

from mrbump.seq_align.simpleSeqID import simpleSeqID
from pyjob import cexec
from pyjob.script import EXE_EXT
from mr_sequence import read_fasta

from Bio import SearchIO


class SequenceHit:
    def __init__(self):
        self.name = None
        self.chainName = None
        self.pdbName = None
        self.chainID = None
        self.domainID = None
        self.rank = None
        self.prob = 0.0
        self.evalue = 0.0
        self.pvalue = 0.0
        self.score = 0.0
        self.ndomains = 0
        self.alignment = ""
        self.alnRange = ""
        self.tarRange = ""
        self.tarExtent = 0
        self.tarMidpoint = 0.0
        self.cols = 0
        self.localSEQID = 0.0
        self.overallSEQID = 0.0
        self.targetAlignment = None
        self.alignments = dict([])

    def __str__(self):
        attrs = [k for k in self.__dict__.keys() if not k.startswith('_')]
        INDENT = "  "
        out_str = "Class: {}\nData:\n".format(self.__class__)
        for a in sorted(attrs):
            out_str += INDENT + "{} : {}\n".format(a, self.__dict__[a])
        return out_str

def find_hits(seqin):
    assert os.path.isfile(seqin), "Cannot find input file: %s" % seqin
    phmmer_logfile = run_phmmer(seqin)
    targetSequence = read_fasta(seqin)
    
    io = SearchIO.read(phmmer_logfile, 'hmmer3-text')
    included = io.hit_filter(lambda x: x.is_included)
    
    hitDict = OrderedDict()
    rank = 0
    for i, hit in enumerate(included):
        rank = i + 1
        for hsp in hit.hsps:
#             sequence_identity = float(sum([a==b for a, b in zip(*[str(s.upper().seq) for s in hsp.aln.get_all_seqs()])])) / float(hsp.aln_span) 
            ph = SequenceHit()
            ph.rank = rank
            ph.chainName = hsp.hit_id
            name, chain = hsp.hit_id.split('_')
            ph.pdbName = name
            ph.chainID = chain
            
            #ph.score = hsp.bitscore
            ph.score = hit.bitscore
            ph.evalue = hsp.evalue # is i-Evalue - possibly evalue_cond in later BioPython
            ph.ndomains = len(hit)
            
#             hstart, hstop = hsp.hit_range
            hstart = hsp.hit_start
            hstop = hsp.hit_end
            ph.alnRange = "{}-{}".format(hstart + 1, hstop)
            qstart, qstop = hsp.query_range
            qstart_p1 = qstart + 1
            ph.tarRange = "{}-{}".format(qstart_p1, qstop)
            ph.tarExtent = hsp.query_span - 1
            ph.tarMidpoint = ((float(qstop) - float(qstart_p1)) / 2.0) + float(qstart_p1)
            
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

def sort_hits_by_size(hits, ascending=False):
    reverse = not(ascending)
    return OrderedDict(sorted(hits.items(), key=lambda x: x[1].tarExtent, reverse=reverse))

def run_phmmer(seqin, dblvl=100):
    logfile = "phmmer.log"
    alnfile = "phmmerAlignment.log"
    phmmerTblout = "phmmerTblout.log"
    phmmerDomTblout = "phmmerDomTblout.log"
    phmmerEXE = os.path.join(os.environ["CCP4"], "libexec", "phmmer")
    if dblvl == 70:
        pfile = "pdb70.txt"
    else:
        pfile = "pdb100.txt"
    seqdb = os.path.join(os.environ["CCP4"], "share", "mrbump", "data", pfile)
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
