'''
Created on 18 Oct 2018

@author: jmht
'''
from ample.util.sequence_util import Sequence
from collections import OrderedDict
from mrbump.seq_align.MRBUMP_phmmer import PHHit, Domains
from mrbump.seq_align.simpleSeqID import simpleSeqID
import os
from pyjob import cexec
from pyjob.script import EXE_EXT
from simbad.util.pdb_util import PdbStructure

import phaser

import pandas as pd

from Bio import SearchIO

# from Bio import AlignIO
# from Bio.SeqFeature import SeqFeature, FeatureLocation
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


class HomologData(object):
    def __init__(self):
        self.name = None
        self.eLLG = None
        self.frac_scat = None
        self.total_frac_scat = None
        self.total_frac_scat_known = None
        self.rmsd = None
        self.ncopies = None
        self.molecular_weight = None
        self.seqid = None
        self.pdb = None
        self.domain = None
        self.range = None

    def __str__(self):
        attrs = [k for k in self.__dict__.keys() if not k.startswith('_')]
        INDENT = "  "
        out_str = "Class: {}\nData:\n".format(self.__class__)
        for a in sorted(attrs):
            out_str += INDENT + "{} : {}\n".format(a, self.__dict__[a])
        return out_str


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
            homologs[hit.name] = HomologData()
            homologs[hit.name].name = hit.name
            homologs[hit.name].pdb = fpath
            homologs[hit.name].molecular_weight = float(pdb_struct.molecular_weight)
            homologs[hit.name].seqid = hits[hit.name].localSEQID / 100.0
            homologs[hit.name].domain = domain.ID
            homologs[hit.name].range = hits[hit.name].tarRange
            pdb_struct.save(fpath)
    return homologs


def ellg_data_from_phaser_log(fpath, homologs):
    with open(fpath) as fh:
        line = fh.readline()
        while line:
            # Get base homolog data
            if line.strip().startswith('eLLG: eLLG of ensemble alone'):
                fh.readline()
                while True:
                    line = fh.readline()
                    if not line.strip():
                        break
                    eLLG, rmsd, frac_scat, name = line.strip().split()
                    h = homologs[name]
                    h.eLLG = float(eLLG)
                    h.rmsd = float(rmsd)
                    h.frac_scat = float(frac_scat)
            # Get ncoopies
            if line.strip().startswith('Number of copies for eLLG target'):
                fh.readline()
                fh.readline()
                while True:
                    line = fh.readline()
                    if not line.strip():
                        break
                    _, _, total_frac_scat_known, total_frac_scat, ncopies, homolog = line.strip().split()
                    h = homologs[homolog]
                    h.total_frac_scat_known = float(total_frac_scat_known)
                    h.total_frac_scat = float(total_frac_scat)
                    h.ncopies = int(ncopies)
            line = fh.readline()
    return homologs


def calculate_ellg(homologs):
    """Stuff from : ccp4-src-2016-02-10/checkout/cctbx-phaser-dials-2015-12-22/phaser/phaser/CalcCCFromMRsolutions.py"""
    mrinput = phaser.InputMR_DAT()
    HKLIN = '../data/2uvo_pdbredo.mtz'
    mrinput.setHKLI(HKLIN)
    F = 'FP'
    SIGF = 'SIGFP'
    mrinput.setLABI_F_SIGF(F, SIGF)
    
    datrun = phaser.runMR_DAT(mrinput)
    
    if not datrun.Success():
        raise RuntimeError("NO SUCCESS")
    
    ellginput = phaser.InputMR_ELLG()
    ellginput.setSPAC_HALL(datrun.getSpaceGroupHall())
    ellginput.setCELL6(datrun.getUnitCell())
    ellginput.setREFL_DATA(datrun.getDATA())
    # Can't mute or no logfile!
    #ellginput.setMUTE(True)
    
    asu_mw = 72846.44

    # Should calculate MW without the search model so that the total MW will be correct when we add the search model
    ellginput.addCOMP_PROT_MW_NUM(asu_mw, 1)
    search_models = []
    for hname, d in homologs.items():
        ellginput.addENSE_PDB_ID(hname, d.pdb, d.seqid)
        search_models.append(hname)
    ellginput.addSEAR_ENSE_OR_ENSE_NUM(search_models, 1)
    
    runellg = phaser.runMR_ELLG(ellginput)
    
    """
    DIR runellg ['ErrorMessage', 'ErrorName', 'Failed', 'Failure', 'Success',  'concise', 'cpu_time', 
    'get_ellg_full_resolution', 'get_ellg_target_nres', 'get_map_chain_ellg', 'get_map_ellg_full_resolution', 
    'get_perfect_data_resolution', 'get_target_resolution', 'get_useful_resolution', 'git_branchname', 'git_commitdatetime', 
    'git_hash', 'git_rev', 'git_shorthash', 'git_tagname', 'git_totalcommits', 'logfile', 'loggraph', 'output_strings', 'process', 
    'run_time', 'setLevel', 'setPackageCCP4', 'setPackagePhenix', 'setPhenixCallback', 'setPhenixPackageCallback', 'set_callback', 
    'set_file_object', 'set_sys_stdout', 'summary', 'svn_revision', 'verbose', 'version_number', 'warnings']
    """
    
    stroutput = runellg.logfile()
    phaser_log = 'phaser1.log'
    with open(phaser_log, 'w') as w:
        w.write(stroutput)
    ellg_data_from_phaser_log(phaser_log, homologs)


class SearchModelFinder(object):
    def __init__(self, seqin):
        self.seqin = seqin
        self.find_homologs()
    
    def find_homologs(self, mock=True):
#         import pickle
#         with open('../tests/homologs.pkl') as f:
#             self.homologs = pickle.load(f)
#             return
        hits = find_hits(self.seqin)
        domains = DomainFinder().find_domains_from_hits(hits)
        homologs = get_homologs(hits, domains)
        if mock:
            ellg_data_from_phaser_log('phaser1.log', homologs)
        else:
            calculate_ellg(homologs)
        self.homologs = homologs
        
    def as_dataframe(self):
        homolog_dict = [h.__dict__ for _, h in self.homologs.items()]
        columns = ['name', 'domain', 'range', 'eLLG', 'ncopies', 'molecular_weight', 'rmsd', 'seqid',
                   'frac_scat', 'total_frac_scat', 'total_frac_scat_known', 'pdb']
        df = pd.DataFrame(homolog_dict, columns=columns)
        df.sort_values('eLLG', inplace=True, ascending=False)
        return df
    
    def as_html(self):
        return self.as_dataframe().to_html(index=False)


    

