"""
Created on 18 Oct 2018

@author: jmht & hlasimpk & rmk65
"""
from Bio import SearchIO
from collections import OrderedDict
import json
import logging
import numpy as np
import os
from pathlib import Path
from pyjob.script import EXE_EXT
import requests
import uuid
import time, random

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
        self.data_created = None
        self.model_url = None
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


def find_hits(seq_info, search_engine=PHMMER, hhsearch_exe=None, hhsearch_db=None, afdb_seqdb=None, pdb_seqdb=None, phmmer_dblvl=95, use_api=False, max_hits=10, nproc=1):
    target_sequence = seq_info.sequence
    af2 = False
    dbtype = None
    if search_engine == PHMMER:
        if use_api:
            if phmmer_dblvl == "af2":
                logger.info("Attempting to search alphafold database search through 3DBeacons API..")
                try:
                    results = run_3dbeacons_alphafold_api(seq_info, max_hits=max_hits)
                    hits = _find_api_hits(results, max_hits=max_hits)
                    return hits
                except json.JSONDecodeError:
                    logger.debug("Phmmer API unavailable, running local phmmer search of AFDB")
                    af2 = True
        else:
            if phmmer_dblvl == "af2":
                logger.info("Running phmmer alphafold database search locally..")
                if afdb_seqdb is not None:
                    logger.info("Database file: %s" % afdb_seqdb)
                else:
                    logger.info("Using CCP4 afdb sequence file..")
                af2 = True
            else:
                logger.info("Running phmmer pdb database search locally..")
                if pdb_seqdb is not None:
                    logger.info("Database file: %s" % pdb_seqdb)
                else:
                    logger.info("Using CCP4 pdb sequence file..")
        logfile, dbtype = run_phmmer(seq_info, afdb_seqdb=afdb_seqdb, pdb_seqdb=pdb_seqdb, dblvl=phmmer_dblvl, nproc=nproc)
        searchio_type = 'hmmer3-text'
    elif search_engine == HHSEARCH:
        searchio_type = 'hhsuite2-text'
        logfile = run_hhsearch(seq_info, hhsearch_exe, hhsearch_db)
    else:
        raise RuntimeError(f"Unrecognised search_engine: {search_engine}")
    return _find_hits(logfile=logfile, searchio_type=searchio_type, target_sequence=target_sequence, af2=af2, max_hits=max_hits, dbtype=dbtype)


def _find_hits(logfile=None, searchio_type=None, target_sequence=None, af2=False, max_hits=10, dbtype=None):
    assert logfile and searchio_type and target_sequence

    if not af2:
        #startT=time.time()
        # Read in the header meta data from the PDB ALL database file
        from mrbump.tools import MRBUMP_utils
        gr = MRBUMP_utils.getPDBres()
        seqMetaDB=gr.readPDBALL()
        #print("Time to read sequence meta data: %.2lf seconds" % (time.time()-startT))

    hitDict = OrderedDict()
    if af2 or searchio_type == "hmmer3-text":
        if af2:
            fix_af_phmmer_log(logfile, "phmmer_af2_fixed.log")
            logfile = "phmmer_af2_fixed.log"
    
        # Read logfile with searchDB
        from mrparse.searchDB import phmmer  
    
        plog = open(logfile, "r")
        phmmerALNLog = plog.readlines()
        plog.close()
    
        phr=phmmer()
        phr.logfile=logfile
        if af2:
            phr.getPhmmerAlignments(targetSequence=target_sequence, phmmerALNLog=phmmerALNLog, PDBLOCAL=None, DB=dbtype, seqMetaDB=None)
        else:
            phr.getPhmmerAlignments(targetSequence=target_sequence, phmmerALNLog=phmmerALNLog, PDBLOCAL=None, DB='PDB', seqMetaDB=seqMetaDB)
        for hitname in (phr.resultsDict):
    
            sh = SequenceHit()
            sh.rank = phr.resultsDict[hitname].rank
            if af2:
                sh.pdb_id = phr.resultsDict[hitname].afdbName
            else:
                sh.pdb_id, sh.chain_id = phr.resultsDict[hitname].afdbName, phr.resultsDict[hitname].chainID
    
            sh.evalue = phr.resultsDict[hitname].evalue
    
            sh.query_start = phr.resultsDict[hitname].tarRange[0]
            sh.query_stop = phr.resultsDict[hitname].tarRange[1]
            sh.hit_start = int(phr.resultsDict[hitname].alnRange[0])
            sh.hit_stop = int(phr.resultsDict[hitname].alnRange[1])
            sh.target_alignment = phr.resultsDict[hitname].targetAlignment
            sh.alignment = phr.resultsDict[hitname].alignment

            hstart = sh.hit_start
            hstop = sh.hit_stop
            qstart, qstop = sh.query_start, sh.query_stop
    
            seq_ali = zip(range(qstart, qstop), sh.alignment)
            sh.seq_ali = [x[0] for x in seq_ali if x[1] != '-']
    
            local, overall = phr.resultsDict[hitname].localSEQID, phr.resultsDict[hitname].overallSEQID
            sh.local_sequence_identity = np.round(local)
            sh.overall_sequence_identity = np.round(overall)
    
            if af2:
                sh.score = phr.resultsDict[hitname].score
                hit_name = phr.resultsDict[hitname].afdbName # + "_" + str(phr.resultsDict[hitname].domainID)
                sh.search_engine = "phmmer"
            elif searchio_type == "hmmer3-text":
                sh.score = phr.resultsDict[hitname].score
                hit_name = phr.resultsDict[hitname].afdbName + "_" + phr.resultsDict[hitname].chainID # + "_" + str(phr.resultsDict[hitname].domainID)
                sh.search_engine = "phmmer"
            sh.name = hit_name
            if sh.rank <= max_hits:
                hitDict[hit_name] = sh
        
    else:
        try:
            io = SearchIO.read(logfile, searchio_type)
        except ValueError:
            logger.exception("ValueError while running Biopython")
            raise RuntimeError('Problem running Biopython SearchIO - you may need to update your version of Biopython.')
        except:
            logger.exception("Unexpected error while running Biopython")
            raise
    
        included = io.hit_filter(lambda x: x.is_included)
        for i, hit in enumerate(included):
            rank = i + 1
            for hsp in hit.hsps:
                sh = SequenceHit()
                sh.rank = rank
                if af2:
                    sh.pdb_id = hsp.hit_id.split("-")[1]
                else:
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
                sh.local_sequence_identity = np.round(local)
                sh.overall_sequence_identity = np.round(overall)
    
                sh.score = hit.score
                hit_name = hit.id + "_" + str(hsp.output_index)
                sh.search_engine = "hhsearch"
                sh.name = hit_name
                if sh.rank <= max_hits:
                    hitDict[hit_name] = sh

    return hitDict


def _find_api_hits(json_list, max_hits=10):
    hitDict = OrderedDict()
    for i, structure in enumerate(json_list):
        hsp_data = structure[0][0]
        data = structure[1]['summary']

        try:
            sh = SequenceHit()
            sh.rank = i + 1
            sh.pdb_id = data['model_identifier']
            sh.evalue = hsp_data['hsp_expect']

            sh.overall_sequence_identity = data['sequence_identity']
            sh.score = hsp_data['hsp_score']
            sh.date_created = data['created']
            sh.model_url = data['model_url']
            hit_name = data['entities'][0]['identifier']
            sh.name = hit_name
            sh.search_engine = "3D Beacons"
            if sh.rank <= max_hits:
                hitDict[hit_name] = sh
        except Exception:
            print(f"Issue with target {data['entities'][0]['identifier']}")
    return hitDict


def fix_af_phmmer_log(input_log, output_log):
    """Function to fix locally run alphafold phmmer results when there are duplicate entries"""
    last_entry = ""
    with open(output_log, 'w') as f_out:
        with open(input_log, 'r') as f_in:
            line = f_in.readline()
            while line:
                if "AFDB release_date" in line:
                    entry = line.split()[8]
                    if entry != last_entry:
                        last_entry = entry
                    else:
                        line = line.replace(entry, f"{entry}_{uuid.uuid4()}")
                f_out.write(line)
                line = f_in.readline()


def sort_hits_by_size(hits, ascending=False):
    reverse = not (ascending)
    return OrderedDict(sorted(hits.items(), key=lambda x: x[1].length, reverse=reverse))


def run_phmmer(seq_info, afdb_seqdb=None, pdb_seqdb=None, dblvl=95, nproc=1):
    logfile = f"phmmer_{dblvl}.log"
    alnfile = f"phmmerAlignment_{dblvl}.log"
    phmmerTblout = f"phmmerTblout_{dblvl}.log"
    phmmerDomTblout = f"phmmerDomTblout_{dblvl}.log"
    phmmerEXE = Path(os.environ["CCP4"], "libexec", "phmmer")
    delete_db = False
    dbtype = None
    if dblvl == "af2":
        if afdb_seqdb is not None:
            seqdb = afdb_seqdb
            dbtype = "AFDB"
        else:
            seqdb = Path(os.environ["CCP4"], "share", "mrbump", "data", "afdb.fasta")
            dbtype= "AFCCP4"
    else:
        sb = makeSeqDB.sequenceDatabase()
        if pdb_seqdb is not None:
            seq_protein_file=Path(os.environ["CCP4_SCR"], "pdb_seqres_protein_%s.txt" % random.randint(0, 9999999)) 
            get_seqres_protein(pdb_seqdb, seq_protein_file)
            if os.path.isfile(seq_protein_file):
                seqdb = seq_protein_file
            else:
                logger.exception("Failed to create PDB sequence file: %s" % seq_protein_file)
                raise
            dbtype= "PDB"
            delete_db = True
        else:
            seqdb = Path(sb.makePhmmerFasta(RLEVEL=dblvl))
            dbtype= "PDBCCP4"
            delete_db = True

    if afdb_seqdb is not None and dblvl == "af2":
        cmd = [str(phmmerEXE) + EXE_EXT,
           '--notextw',
           '--tblout', phmmerTblout,
           '--domtblout', phmmerDomTblout,
           '--F1', '1e-15',
           '--F2', '1e-15',
           '--cpu', str(nproc),
           '-A', alnfile,
           str(seq_info.sequence_file), str(seqdb)]
    else:
        cmd = [str(phmmerEXE) + EXE_EXT,
           '--notextw',
           '--tblout', phmmerTblout,
           '--domtblout', phmmerDomTblout,
           '--cpu', str(nproc),
           '-A', alnfile,
           str(seq_info.sequence_file), str(seqdb)]
    stdout = run_cmd(cmd)
    if os.name == 'nt':
        lines = stdout.split('\n')
        lines[0] = "# phmmer :: search a protein sequence against a protein database"
        stdout = "\n".join(lines)
    with open(logfile, 'w') as f_out:
        f_out.write(stdout)

    if delete_db:
        seqdb.unlink()

    return logfile, dbtype


def run_hhsearch(seq_info, hhsearch_exe, hhsearch_db):
    logfile = "hhsearch.log"
    hhsearch_db = Path(hhsearch_db)
    cmd = [hhsearch_exe,
           '-i', seq_info.sequence_file,
           '-d', str(hhsearch_db.joinpath(hhsearch_db.stem)),
           '-o', logfile]
    run_cmd(cmd)
    return logfile


def run_3dbeacons_alphafold_api(seq_info, max_hits=10):
    url = 'https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/sequence/search'
    headers = {
    'accept': 'application/json',
    'Content-Type': 'application/json'
    }
    data = {
    'sequence': seq_info.sequence,
    }

    response = requests.post(url, headers=headers, json=data)
    job_id = json.loads(response.text)['job_id']

    url = 'https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/sequence/result?job_id=' + job_id


    url = 'https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/sequence/result'
    params = {
        'job_id': job_id,
    }
    headers = {
        'accept': 'application/json'
    }

    response = requests.get(url, params=params, headers=headers)

    data = json.loads(response.text)

    alphafold_structures = []
    for entry in data:
        hsp_data = entry['hit_hsps']
        for structure in entry['summary']['structures']:
            provider = structure['summary']['provider']
            if provider == 'AlphaFold DB':
                alphafold_structures.append([hsp_data, structure])

    return alphafold_structures

def get_seqres_protein(pdbseqfile, outfile):
    """ extract the protein sequences from the full pdb_seqres.txt file """

    inf=open(pdbseqfile, "r")
    inflines=inf.readlines()
    inf.close()

    protlist=""
    count=0
    for line in inflines:
        if ">" in line and "mol:na" not in line:
            protlist+=line
            protlist+=inflines[count+1]
        count+=1

    outf=open(outfile, "w")
    outf.write(protlist)
    outf.close()
