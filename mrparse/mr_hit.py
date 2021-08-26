"""
Created on 18 Oct 2018

@author: jmht
"""
from collections import OrderedDict
import logging
import os

from mrparse.mr_util import run_cmd, EXE_EXT
from mrbump.seq_align.simpleSeqID import simpleSeqID
from Bio import SearchIO

PHMMER = 'phmmer'
HHSEARCH = 'hhsearch'

logger = logging.getLogger(__name__)


class SequenceHit:
    def __init__(self):
        self.name = None
        self.pdb_id = None
        self.chain_id = None
        self.rank = None
        self.prob = 0.0
        self.evalue = 0.0
        self.pvalue = 0.0
        self.score = 0.0
        self.seq_ali = None
        self.alignment = ""
        self.query_start = None 
        self.query_stop = None
        self.hit_start = None
        self.hit_stop = None
        self.overall_sequence_identity = 0.0
        self.target_alignment = None
        self.alignments = dict([])
        # pointers to objects
        self.region = None
        self.homolog = None
    
    @property
    def length(self):
        return self.query_extent
    
    @property
    def hit_extent(self):
        """python list indexing so need to subtract 1"""
        return (self.hit_stop - self.hit_start) - 1
    
    @property
    def hit_range(self):
        return (self.hit_start, self.hit_stop)
    
    @property
    def query_extent(self):
        """python list indexing so need to add 1"""
        return (self.query_stop - self.query_start) - 1
    
    @property
    def query_midpoint(self):
        return ((float(self.query_stop) - float(self.query_start + 1)) / 2.0) + float(self.query_start + 1)
    
    @property
    def query_range(self):
        return (self.query_start, self.query_stop)
    
    @property
    def region_id(self):
        return self._get_child_attr('region', 'id')
    
    @property
    def region_index(self):
        return self._get_child_attr('region', 'index')

    # Non-property methods
    def _get_child_attr(self, child, attr):
        if hasattr(self, child):
            child = getattr(self, child)
            if hasattr(child, attr):
                return getattr(child, attr)
        return None

    def __str__(self):
        attrs = [k for k in self.__dict__.keys() if not k.startswith('_')]
        INDENT = "  "
        out_str = "Class: {}\nData:\n".format(self.__class__)
        for a in sorted(attrs):
            out_str += INDENT + "{} : {}\n".format(a, self.__dict__[a])
        return out_str


def find_hits(seq_info, search_engine=PHMMER, dblvl=95):
    if search_engine == PHMMER:
        logfile = run_phmmer(seq_info, dblvl=dblvl)
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
            sh = SequenceHit()
            sh.rank = rank
            if "_" in hsp.hit_id:
                name, chain = hsp.hit_id.split('_')
                sh.pdb_id = name
                sh.chain_id = chain
            else:
                sh.pdb_id = hsp.hit_id
                sh.chain_id = "A"
    #sh.score = hsp.bitscore
            sh.score = hit.bitscore
            sh.evalue = hsp.evalue # is i-Evalue - possibly evalue_cond in later BioPython
            hstart = hsp.hit_start
            hstop = hsp.hit_end
            qstart, qstop = hsp.query_range
            seq_ali = zip(range(qstart, qstop), hsp.hit.seq)
            sh.seq_ali = [x[0] for x in seq_ali if x[1] != '-']
            sh.query_start = qstart
            sh.query_stop = qstop
            sh.hit_start = hstart
            sh.hit_stop = hstop
            target_alignment = "".join(hsp.aln[0].upper()) # assume the first Sequence is always the target
            sh.target_alignment = target_alignment
            alignment = "".join(hsp.aln[1].upper()) # assume the first Sequence is always the target
            sh.alignment = alignment
            local, overall = simpleSeqID().getPercent(alignment, target_alignment, target_sequence)
            sh.local_sequence_identity = local
            sh.overall_sequence_identity = overall
            hit_name = hit.id + "_" + str(hsp.domain_index)
            sh.name = hit_name
            hitDict[hit_name] = sh
    return hitDict


def sort_hits_by_size(hits, ascending=False):
    reverse = not(ascending)
    return OrderedDict(sorted(hits.items(), key=lambda x: x[1].length, reverse=reverse))


def run_phmmer(seq_info, dblvl=95):
    logfile = "phmmer_{}.log".format(dblvl)
    alnfile = "phmmerAlignment_{}.log".format(dblvl)
    phmmerTblout = "phmmerTblout_{}.log".format(dblvl)
    phmmerDomTblout = "phmmerDomTblout_{}.log".format(dblvl)
    phmmerEXE = os.path.join(os.environ["CCP4"], "libexec", "phmmer")
    if dblvl == 95:
        seqdb = os.path.join(os.environ["CCP4"], "share", "mrbump", "data", "pdb95.txt")
    elif dblvl == "af2":
        seqdb = os.path.join(os.environ["CCP4"], "share", "mrparse", "data", "af2_sequences.fasta")
    else:
        seqdb = os.path.join(os.environ["CCP4"], "share", "mrbump", "data", "pdb100.txt")
    cmd = [phmmerEXE + EXE_EXT,
           '--notextw',
           '--tblout', phmmerTblout,
           '--domtblout', phmmerDomTblout,
           '-A', alnfile,
           seq_info.sequence_file, seqdb]
    stdout = run_cmd(cmd)
    if os.name == 'nt':
        lines = stdout.split('\r\n')
        lines[0] = "# phmmer :: search a protein sequence against a protein database"
        stdout = "\n".join(lines)
    with open(logfile, 'w') as f_out:
        f_out.write(stdout)
    return logfile


def run_hhsearch(seq_info):
    raise NotImplementedError("{} search needs to be implemented".format(HHSEARCH))

