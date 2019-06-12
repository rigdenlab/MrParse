'''
Created on 18 Oct 2018

@author: jmht
'''
from collections import OrderedDict
import logging
import os

from mrbump.seq_align.simpleSeqID import simpleSeqID
from pyjob import cexec
from pyjob.script import EXE_EXT
from Bio import SearchIO

PHMMER = 'phmmer'
HHSEARCH = 'hhsearch'

logger = logging.getLogger(__name__)

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
        self.alnStart = None # Count from 1
        self.alnStop = None
        self.tarStart = None # Count from 1
        self.tarStop = None
        self.tarExtent = 0
        self.tarMidpoint = 0.0
        self.cols = 0
        self.local_sequence_identity = 0.0
        self.overall_sequence_identity = 0.0
        self.targetAlignment = None
        self.alignments = dict([])
        self.region = None
        self._homolog = None
    
    @property
    def alnRange(self):
        return "{}-{}".format(self.alnStart, self.alnStop)
    
    @property
    def length(self):
        return self.alnStop - self.alnStart
    
    @property
    def regionId(self):
        if self.region is not None and hasattr(self.region, 'ID'):
            return self.region.ID
        return None

    @property
    def tarRange(self):
        return "{}-{}".format(self.tarStart, self.tarStop)

    def __str__(self):
        attrs = [k for k in self.__dict__.keys() if not k.startswith('_')]
        INDENT = "  "
        out_str = "Class: {}\nData:\n".format(self.__class__)
        for a in sorted(attrs):
            out_str += INDENT + "{} : {}\n".format(a, self.__dict__[a])
        return out_str


def find_hits(seq_info, search_engine=PHMMER):
    if search_engine == PHMMER:
        logfile = run_phmmer(seq_info)
        searchio_type = 'hmmer3-text'
    elif search_engine == HHSEARCH:
        searchio_type = 'hhsuite2-text'
        logfile = run_hhsearch(seq_info)
    else:
        raise RuntimeError("Unrecognised search_engine: {}".format(search_engine))
    target_sequence = seq_info.sequence
    return _find_hits(logfile=logfile, searchio_type=searchio_type, target_sequence=target_sequence)
    

def _find_hits(logfile=None, searchio_type=None, target_sequence=None):
    assert logfile and searchio_type and target_sequence
    try:
        io = SearchIO.read(logfile, searchio_type)
    except AttributeError:
        # If the error is: "'NoneType' object has no attribute 'group'" then this is an error in Biopython that
        # has been fixed in later versions
        logger.exception("AttributeError while running Biopython")
        raise RuntimeError('Problem running Biopython SearchIO - you may need to update your version of Biopython.')
    except:
        logger.exception("Unexpected error while running Biopython")
        raise
        
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
            hstart = hsp.hit_start
            hstop = hsp.hit_end
            qstart, qstop = hsp.query_range
            qstart_p1 = qstart + 1
            ph.alnStart = qstart_p1
            ph.alnStop = qstop
            ph.tarStart = hstart
            ph.tarStop = hstop
            ph.tarExtent = hstop - hstart        
            ph.tarMidpoint = ((float(qstop) - float(qstart_p1)) / 2.0) + float(qstart_p1)
            target_alignment = "".join(hsp.aln[0].upper()) # assume the first Sequence is always the target
            ph.targetAlignment = target_alignment
            alignment = "".join(hsp.aln[1].upper()) # assume the first Sequence is always the target
            ph.alignment = alignment
            local, overall = simpleSeqID().getPercent(alignment, target_alignment, target_sequence)
            ph.local_sequence_identity = local
            ph.overall_sequence_identity = overall
            name = hit.id + "_" + str(hsp.domain_index)
            ph.name = name
            hitDict[name] = ph
    return hitDict

def sort_hits_by_size(hits, ascending=False):
    reverse = not(ascending)
    return OrderedDict(sorted(hits.items(), key=lambda x: x[1].tarExtent, reverse=reverse))

def run_phmmer(seq_info, dblvl=95):
    logfile = "phmmer.log"
    alnfile = "phmmerAlignment.log"
    phmmerTblout = "phmmerTblout.log"
    phmmerDomTblout = "phmmerDomTblout.log"
    phmmerEXE = os.path.join(os.environ["CCP4"], "libexec", "phmmer")
    if dblvl == 95:
        pfile = "pdb95.txt"
    else:
        pfile = "pdb100.txt"
    seqdb = os.path.join(os.environ["CCP4"], "share", "mrbump", "data", pfile)
    cmd = [phmmerEXE + EXE_EXT,
           '--notextw',
           '--tblout', phmmerTblout,
           '--domtblout', phmmerDomTblout,
           '-A', alnfile,
           seq_info.sequence_file, seqdb]
    stdout = cexec(cmd)
    with open(logfile, 'w') as f_out:
        f_out.write(stdout)
    return logfile

def run_hhsearch(seq_info):
    raise NotImplementedError("{} search needs to be implemented".format(HHSEARCH))

