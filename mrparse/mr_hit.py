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

from mrparse.mr_alphafold import PdbModelException
from mrparse.mr_util import run_cmd
from mrparse.searchDB import phmmer
from mrbump.seq_align.simpleSeqID import simpleSeqID
from mrbump.tools import makeSeqDB
from simbad.util.pdb_util import PdbStructure

PHMMER = 'phmmer'
HHSEARCH = 'hhsearch'

api_url = "http://localhost:8080/phmmer/v1/jobs"

logger = logging.getLogger(__name__)
logging.getLogger("requests").setLevel(logging.WARNING)


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


def find_hits(seq_info, search_engine=PHMMER, hhsearch_exe=None, hhsearch_db=None, afdb_seqdb=None, bfvd_seqdb=None, esm_seqdb=None, pdb_seqdb=None, phmmer_dblvl=95, max_hits=10, nproc=1, ccp4cloud=False, phmmer_exe=None):
    target_sequence = seq_info.sequence
    dbtype = None
    if search_engine == PHMMER:
        if phmmer_dblvl == "af2":
            logger.info("Running phmmer alphafold database search locally..")
            if afdb_seqdb is not None:
                logger.info("Database file: %s" % afdb_seqdb)
            else:
                logger.info("Using CCP4 afdb sequence file..")
            seqdb = afdb_seqdb
        elif phmmer_dblvl == "esmfold":
            if esm_seqdb is not None:
                logger.info("Running phmmer esmfold database search locally..")
                logger.info("Database file: %s" % esm_seqdb)
                seqdb = esm_seqdb
            else:
                return {}
        elif phmmer_dblvl == "bfvd":
            if bfvd_seqdb is not None:
                logger.info("Running phmmer bfvd database search locally..")
                logger.info("Database file: %s" % bfvd_seqdb)
                seqdb = bfvd_seqdb
            else:
                return {}
        else:
            logger.info("Running phmmer pdb database search locally..")
            if pdb_seqdb is not None:
                logger.info("Database file: %s" % pdb_seqdb)
                seqdb = pdb_seqdb
            else:
                logger.info("Using CCP4 pdb sequence file..")
                seqdb=None
        logfile, dbtype = run_phmmer(seq_info, seqdb=seqdb, dblvl=phmmer_dblvl, nproc=nproc, phmmer_exe=phmmer_exe)
        searchio_type = 'hmmer3-text'
    elif search_engine == HHSEARCH:
        searchio_type = 'hhsuite2-text'
        logfile = run_hhsearch(seq_info, hhsearch_exe, hhsearch_db)
    else:
        raise RuntimeError(f"Unrecognised search_engine: {search_engine}")
    return _find_hits(logfile=logfile, searchio_type=searchio_type, target_sequence=target_sequence, phmmer_dblvl=phmmer_dblvl, max_hits=max_hits, dbtype=dbtype)

# def _find_api_hits(seq_info, max_hits=10, database="af2"):
#     databases = {'af2': 'afdb', 'bfvd': 'bfvd', 'esmfold': 'esmatlas'}

#     job_id = str(uuid.uuid4())
#     data = {
#         "job_id": job_id,
#         "input_sequence": seq_info.sequence,
#         "database": databases[database],
#         "number_of_hits": max_hits,
#         "run_type": "mrparse"
#     }
#     json_file = str(Path(seq_info.sequence_file).parent.resolve().joinpath('phmmer.json'))
#     with open(json_file, 'w') as f:
#         json.dump(data, f)

#     api_file = str(Path(__file__).parent.resolve().joinpath('scripts', 'phmmer_api.py'))
    
#     cmd = ['ccp4-python', api_file, '-i', json_file, '-o', str(Path.cwd())]

#     run_cmd(cmd)
#     with open(Path.cwd().joinpath(f'{job_id}.log'), 'r') as f:
#         data = f.read()

#     results_dict = json.loads(data)

#     hitDict = OrderedDict()
#     for hit in results_dict:
#         sh = SequenceHit()
#         if ':' in hit:
#             name = hit.split(":")[1].rsplit("_", 1)[0]
#         else:
#             name = hit.rsplit("_", 1)[0]
#         result = results_dict[hit]
#         sh.rank = result['rank']
#         sh.pdb_id = name
#         sh.evalue = result['evalue']
#         sh.query_start = result['query_start']
#         sh.query_stop = result['query_stop']
#         sh.hit_start = result['hit_start']
#         sh.hit_stop = result['hit_stop']
#         sh.target_alignment = result['target_alignment']
#         sh.alignment = result['alignment']
#         seq_ali = zip(range(sh.query_start, sh.query_stop), sh.alignment)
#         sh.seq_ali = [x[0] for x in seq_ali if x[1] != '-']
#         sh.local_sequence_identity = result['local_sequence_identity']
#         sh.overall_sequence_identity = result['overall_sequence_identity']
#         sh.search_engine = 'phmmer'
#         sh.name = name
#         hitDict[sh.name] = sh

#     return hitDict

def _find_hits(logfile=None, searchio_type=None, target_sequence=None, phmmer_dblvl=95, max_hits=10, dbtype=None):
    assert logfile and searchio_type and target_sequence

    if phmmer_dblvl != "af2":
        # Read in the header meta data from the PDB ALL database file
        from mrbump.tools import MRBUMP_utils
        gr = MRBUMP_utils.getPDBres()
        seqMetaDB=gr.readPDBALL()

    hitDict = OrderedDict()
    if phmmer_dblvl == "af2" or searchio_type == "hmmer3-text":
        if phmmer_dblvl == "af2":
            fix_af_phmmer_log(logfile, "phmmer_af2_fixed.log")
            logfile = "phmmer_af2_fixed.log"
    
        # Read logfile with searchDB
        from mrparse.searchDB import phmmer  
    
        plog = open(logfile, "r")
        phmmerALNLog = plog.readlines()
        plog.close()
    
        phr=phmmer()
        phr.logfile=logfile
        if phmmer_dblvl == "af2":
            phr.getPhmmerAlignments(targetSequence=target_sequence, phmmerALNLog=phmmerALNLog, PDBLOCAL=None, DB=dbtype, seqMetaDB=None)
        else:
            phr.getPhmmerAlignments(targetSequence=target_sequence, phmmerALNLog=phmmerALNLog, PDBLOCAL=None, DB='PDB', seqMetaDB=seqMetaDB)
        for hitname in (phr.resultsDict):

            if phmmer_dblvl == "bfvd" or phmmer_dblvl == "esmfold":
                if ':' in hitname:
                    name = hitname.split(":")[1].rsplit("_", 1)[0]
                else:
                    name = hitname.rsplit("_", 1)[0]
            else:
                name = phr.resultsDict[hitname].afdbName

            sh = SequenceHit()
            sh.rank = phr.resultsDict[hitname].rank
            if phmmer_dblvl == "af2":
                sh.pdb_id = phr.resultsDict[hitname].afdbName
            else:
                sh.pdb_id, sh.chain_id = name, phr.resultsDict[hitname].chainID
    
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
    
            if phmmer_dblvl == "af2":
                sh.score = phr.resultsDict[hitname].score
                hit_name = phr.resultsDict[hitname].afdbName # + "_" + str(phr.resultsDict[hitname].domainID)
                sh.search_engine = "phmmer"
            elif phmmer_dblvl == "bfvd" or phmmer_dblvl == "esmfold":
                sh.score = phr.resultsDict[hitname].score
                hit_name = name
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
                if phmmer_dblvl == "af2":
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


def run_phmmer(seq_info, seqdb=None, dblvl=95, nproc=1, phmmer_exe=None):
    logfile = f"phmmer_{dblvl}.log"
    alnfile = f"phmmerAlignment_{dblvl}.log"
    phmmerTblout = f"phmmerTblout_{dblvl}.log"
    phmmerDomTblout = f"phmmerDomTblout_{dblvl}.log"
    if phmmer_exe is not None and dblvl == "af2":
        phmmerEXE=phmmer_exe
    else:
        phmmerEXE = Path(os.environ["CCP4"], "libexec", "phmmer")
    delete_db = False
    dbtype = None
    logger.debug("Running phmmer executable: %s" % phmmerEXE)
    logger.debug("Running phmmer database level: %s" % dblvl)
    if dblvl == "af2":
        if seqdb is not None or phmmer_exe is not None:
            dbtype = "AFDB"
        else:
            seqdb = Path(os.environ["CCP4"], "share", "mrbump", "data", "afdb.fasta")
            dbtype= "AFCCP4"
    else:
        sb = makeSeqDB.sequenceDatabase()
        if seqdb is not None:
            seq_protein_file=Path(
                os.environ["CCP4_SCR"], 
                "pdb_seqres_protein_%s.txt" % random.randint(0, 9999999)
            ) 
            get_seqres_protein(seqdb, seq_protein_file)
            if os.path.isfile(seq_protein_file):
                seqdb = seq_protein_file
            else:
                logger.exception(
                    "Failed to create PDB sequence file: %s" % seq_protein_file
                )
                raise
            dbtype= "PDB"
            delete_db = True
        else:
            seqdb = Path(sb.makePhmmerFasta(RLEVEL=dblvl))
            dbtype= "PDBCCP4"
            delete_db = True

    if dblvl == "af2":
        if os.name != 'nt':
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
               '--F1', '1e-15',
               '--F2', '1e-15',
               '-A', alnfile,
               str(seq_info.sequence_file), str(seqdb)]
    else:
        if os.name != 'nt':
            cmd = [str(phmmerEXE) + EXE_EXT,
                       '--notextw',
                       '--tblout', phmmerTblout,
                       '--domtblout', phmmerDomTblout,
                       '--cpu', str(nproc),
                       '-A', alnfile,
                       str(seq_info.sequence_file), str(seqdb)]
        else:
            cmd = [str(phmmerEXE) + EXE_EXT,
                       '--notextw',
                       '--tblout', phmmerTblout,
                       '--domtblout', phmmerDomTblout,
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
           '-i', str(seq_info.sequence_file),
           '-d', str(hhsearch_db.joinpath(hhsearch_db.stem)),
           '-o', logfile]
    run_cmd(cmd)
    return logfile


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
