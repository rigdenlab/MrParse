"""
Created on 18 Oct 2018

@author: jmht & hlasimpk
"""
from Bio import SearchIO
from collections import OrderedDict
import json
import logging
import os
from pathlib import Path
from pyjob.script import EXE_EXT
import requests

from mrparse.mr_util import run_cmd
from mrbump.seq_align.simpleSeqID import simpleSeqID
from mrbump.tools import makeSeqDB

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
        self.search_engine = None
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
        return self.hit_start, self.hit_stop

    @property
    def query_extent(self):
        """python list indexing so need to add 1"""
        return (self.query_stop - self.query_start) - 1

    @property
    def query_midpoint(self):
        return ((float(self.query_stop) - float(self.query_start + 1)) / 2.0) + float(self.query_start + 1)

    @property
    def query_range(self):
        return self.query_start, self.query_stop

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
        indent = "  "
        out_str = f"Class: {self.__class__}\nData:\n"
        for a in sorted(attrs):
            out_str += indent + f"{a} : {self.__dict__[a]}\n"
        return out_str


def find_hits(seq_info, search_engine=PHMMER, hhsearch_exe=None, hhsearch_db=None, phmmer_dblvl=95):
    target_sequence = seq_info.sequence
    if search_engine == PHMMER:
        if phmmer_dblvl == "af2":
            json_file = run_phmmer_alphafold_api(seq_info)
            return _find_json_hits(json_file, target_sequence=seq_info)
        else:
            logfile = run_phmmer(seq_info, dblvl=phmmer_dblvl)
        searchio_type = 'hmmer3-text'
    elif search_engine == HHSEARCH:
        searchio_type = 'hhsuite2-text'
        logfile = run_hhsearch(seq_info, hhsearch_exe, hhsearch_db)
    else:
        raise RuntimeError(f"Unrecognised search_engine: {search_engine}")
    return _find_hits(logfile=logfile, searchio_type=searchio_type, target_sequence=target_sequence)


def _find_hits(logfile=None, searchio_type=None, target_sequence=None):
    assert logfile and searchio_type and target_sequence
    try:
        io = SearchIO.read(logfile, searchio_type)
    except ValueError:
        logger.exception("ValueError while running Biopython")
        raise RuntimeError('Problem running Biopython SearchIO - you may need to update your version of Biopython.')
    except:
        logger.exception("Unexpected error while running Biopython")
        raise

    included = io.hit_filter(lambda x: x.is_included)
    hitDict = OrderedDict()
    for i, hit in enumerate(included):
        rank = i + 1
        for hsp in hit.hsps:
            sh = SequenceHit()
            sh.rank = rank
            sh.pdb_id, sh.chain_id = hsp.hit_id.split('_')
            sh.evalue = hsp.evalue  # is i-Evalue - possibly evalue_cond in later BioPython
            hstart = hsp.hit_start
            hstop = hsp.hit_end
            qstart, qstop = hsp.query_range
            seq_ali = zip(range(qstart, qstop), hsp.hit.seq)
            sh.seq_ali = [x[0] for x in seq_ali if x[1] != '-']
            sh.query_start = qstart
            sh.query_stop = qstop
            sh.hit_start = hstart
            sh.hit_stop = hstop
            target_alignment = "".join(hsp.aln[0].upper())  # assume the first Sequence is always the target
            sh.target_alignment = target_alignment
            alignment = "".join(hsp.aln[1].upper())  # assume the first Sequence is always the target
            sh.alignment = alignment
            local, overall = simpleSeqID().getPercent(alignment, target_alignment, target_sequence)
            sh.local_sequence_identity = local
            sh.overall_sequence_identity = overall

            if searchio_type == "hmmer3-text":
                sh.score = hit.bitscore
                hit_name = hit.id + "_" + str(hsp.domain_index)
                sh.search_engine = "phmmer"
            else:
                sh.score = hit.score
                hit_name = hit.id + "_" + str(hsp.output_index)
                sh.search_engine = "hhsearch"
            sh.name = hit_name
            hitDict[hit_name] = sh
    return hitDict


def _find_json_hits(json_file, target_sequence):
    hitDict = OrderedDict()
    with open(json_file, 'r') as f_in:
        data = json.load(f_in)
        for i, hit in enumerate(data['results']['hits']):
            try:
                sh = SequenceHit()
                sh.rank = i + 1
                sh.pdb_id = hit['name']
                sh.evalue = hit['evalue']

                alignment_info = hit['domains'][0]
                qstart = alignment_info['alihmmfrom']
                qstop = alignment_info['alihmmto']
                sh.query_start = qstart
                sh.query_stop = qstop
                sh.hit_start = alignment_info['alisqfrom']
                sh.hit_stop = alignment_info['alisqto']
                seq_ali = zip(range(qstart, qstop), alignment_info['aliaseq'])
                sh.seq_ali = [x[0] for x in seq_ali if x[1] != '-']
                alignment = alignment_info['aliaseq'].upper()
                target_alignment = alignment_info['alimodel'].upper()
                sh.target_alignment = alignment
                sh.alignment = target_alignment
                local, overall = simpleSeqID().getPercent(alignment, target_alignment, target_sequence)
                sh.local_sequence_identity = local
                sh.overall_sequence_identity = overall

                sh.score = hit['score']
                hit_name = hit['name'] + "_" + str(hit['ndom'])
                sh.name = hit_name
                sh.search_engine = "phmmer"
                hitDict[hit_name] = sh
            except Exception:
                logger.debug(f"Issue with target {hit['name']}")
    return hitDict


def sort_hits_by_size(hits, ascending=False):
    reverse = not (ascending)
    return OrderedDict(sorted(hits.items(), key=lambda x: x[1].length, reverse=reverse))


def run_phmmer(seq_info, dblvl=95):
    logfile = f"phmmer_{dblvl}.log"
    alnfile = f"phmmerAlignment_{dblvl}.log"
    phmmerTblout = f"phmmerTblout_{dblvl}.log"
    phmmerDomTblout = f"phmmerDomTblout_{dblvl}.log"
    phmmerEXE = Path(os.environ["CCP4"], "libexec", "phmmer")
    delete_db = False
    if dblvl == "af2":
        seqdb = Path(os.environ["CCP4"], "share", "mrbump", "data", "afdb.fasta")
    else:
        sb = makeSeqDB.sequenceDatabase()
        seqdb = Path(sb.makePhmmerFasta(RLEVEL=dblvl))
        delete_db = True

    cmd = [str(phmmerEXE) + EXE_EXT,
           '--notextw',
           '--tblout', phmmerTblout,
           '--domtblout', phmmerDomTblout,
           '-A', alnfile,
           seq_info.sequence_file, str(seqdb)]
    stdout = run_cmd(cmd)
    if os.name == 'nt':
        lines = stdout.split('\r\n')
        lines[0] = "# phmmer :: search a protein sequence against a protein database"
        stdout = "\n".join(lines)
    with open(logfile, 'w') as f_out:
        f_out.write(stdout)

    if delete_db:
        seqdb.unlink()

    return logfile


def run_hhsearch(seq_info, hhsearch_exe, hhsearch_db):
    logfile = "hhsearch.log"
    hhsearch_db = Path(hhsearch_db)
    cmd = [hhsearch_exe,
           '-i', seq_info.sequence_file,
           '-d', str(hhsearch_db.joinpath(hhsearch_db.stem)),
           '-o', logfile]
    run_cmd(cmd)
    return logfile


def run_phmmer_alphafold_api(seq_info):
    params = {
        'seqdb': 'alphafold',
        'seq': f'>Seq\n{seq_info.sequence}'
    }

    # post the seqrch request to the server
    results_url = requests.post(url='https://www.ebi.ac.uk/Tools/hmmer/search/phmmer', data=params).url

    # modify the range, format and presence of alignments in your results here
    res_params = {
        'output': 'json',
        'range': '1,10'
    }

    # send a GET request to the server
    data = requests.get(results_url, params=res_params)

    logfile = "phmmer_afdb.json"
    with open(logfile, 'w') as f_out:
        f_out.write(data.text)

    return logfile
