"""
Created on 18 Oct 2018

@author: jmht
"""
import os
import logging
from mrparse import mr_homolog 
from mrparse import mr_alphafold
from mrparse import mr_bfvd
from mrparse import mr_esmatlas
from mrparse import mr_hit
from mrparse.mr_region import RegionFinder
from mrparse import mr_pfam
from mrparse.mr_util import now

logger = logging.getLogger(__name__)


class SearchModelFinder(object):
    def __init__(self, seq_info, **kwargs):
        self.seq_info = seq_info
        self.hkl_info = kwargs.get("hkl_info", None)
        self.pdb_dir = kwargs.get("pdb_dir", None)
        self.pdb_local = kwargs.get("pdb_local", None)
        self.search_engine = kwargs.get("search_engine", "phmmer")
        self.phmmer_dblvl = kwargs.get("phmmer_dblvl", 95)
        self.plddt_cutoff = kwargs.get("plddt_cutoff", 70)
        self.hhsearch_exe = kwargs.get("hhsearch_exe", None)
        self.hhsearch_db = kwargs.get("hhsearch_db", None)
        self.afdb_seqdb = kwargs.get("afdb_seqdb", None)
        self.bfvd_seqdb = kwargs.get("bfvd_seqdb", None)
        self.esm_seqdb = kwargs.get("esm_seqdb", None)
        self.pdb_seqdb = kwargs.get("pdb_seqdb", None)
        self.use_api = kwargs.get("use_api", False)
        self.max_hits = kwargs.get("max_hits", 10)
        self.database = kwargs.get("database", "all")
        self.nproc = kwargs.get("nproc", 1)
        self.phmmer_exe = kwargs.get("phmmer_exe", None)
        self.ccp4cloud = kwargs.get("ccp4cloud", False)
        self.hits = None
        self.af_model_hits = None
        self.bfvd_model_hits = None
        self.esm_model_hits = None
        self.regions = None
        self.af_model_regions = None
        self.bfvd_model_regions = None
        self.esm_model_regions = None
        self.homologs = {}
        self.af_models = {}
        self.bfvd_models = {}
        self.esm_models = {}

    def __call__(self):
        """Required so that we can use multiprocessing pool. We need to be able to pickle the object passed
        to the pool and instance methods don't work, so we add the object to the pool and define __call__
        https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map/6975654#6975654
        """
        if self.database in ["all", "pdb"]:
            logger.debug(f'SearchModelFinder started at {now()}')
            self.find_homolog_regions()
            logger.debug(f'SearchModelFinder homolog regions done at {now()}')
            self.prepare_homologs()
            logger.debug(f'SearchModelFinder homologs done at {now()}')
        if self.database in ["all", "afdb"]: 
            self.find_af2_model_regions()
            logger.debug(f'SearchModelFinder AF2 model regions done at {now()}')
            self.prepare_af2_models()
            logger.debug(f'SearchModelFinder AF2 models done at {now()}')
        if self.database in ['all', 'bfvd']:
            self.find_bfvd_model_regions()
            logger.debug(f'SearchModelFinder BFVD model regions done at {now()}')
            self.prepare_bfvd_models()
            logger.debug(f'SearchModelFinder BFVD models done at {now()}')
        if self.database in ["all", "esmfold"]:
            self.find_esm_model_regions()
            logger.debug(f'SearchModelFinder ESM model regions done at {now()}')
            self.prepare_esm_model()
            logger.debug(f'SearchModelFinder ESM models done at {now()}')
        return self
    
    def find_homolog_regions(self):
        self.hits = mr_hit.find_hits(self.seq_info, search_engine=self.search_engine, hhsearch_exe=self.hhsearch_exe,
                                     hhsearch_db=self.hhsearch_db, afdb_seqdb=self.afdb_seqdb, bfvd_seqdb=self.bfvd_seqdb,
                                     esm_seqdb=self.esm_seqdb, pdb_seqdb=self.pdb_seqdb,
                                     phmmer_dblvl=self.phmmer_dblvl, max_hits=self.max_hits,
                                     nproc=self.nproc, ccp4cloud=self.ccp4cloud, phmmer_exe=self.phmmer_exe)
        if not self.hits:
            logger.critical('SearchModelFinder PDB search could not find any hits!')
            return None
        self.regions = RegionFinder().find_regions_from_hits(self.hits)
        return self.regions

    def find_af2_model_regions(self):
        self.af_model_hits = mr_hit.find_hits(self.seq_info, search_engine="phmmer", hhsearch_exe=None,
                                              hhsearch_db=None, afdb_seqdb=self.afdb_seqdb, bfvd_seqdb=self.bfvd_seqdb,
                                              esm_seqdb=self.esm_seqdb, pdb_seqdb=self.pdb_seqdb,
                                              phmmer_dblvl="af2", max_hits=self.max_hits,
                                              nproc=self.nproc, ccp4cloud=self.ccp4cloud, phmmer_exe=self.phmmer_exe)
        if not self.af_model_hits:
            logger.critical('SearchModelFinder EBI Alphafold database search could not find any hits!')
            return None
        self.af_model_regions = RegionFinder().find_regions_from_hits(self.af_model_hits)
        return self.af_model_regions
    
    def find_bfvd_model_regions(self):
        self.bfvd_model_hits = mr_hit.find_hits(self.seq_info, search_engine="phmmer", hhsearch_exe=None,
                                              hhsearch_db=None, afdb_seqdb=self.afdb_seqdb, bfvd_seqdb=self.bfvd_seqdb,
                                              esm_seqdb=self.esm_seqdb, pdb_seqdb=self.pdb_seqdb,
                                              phmmer_dblvl="bfvd", max_hits=self.max_hits,
                                              nproc=self.nproc, ccp4cloud=self.ccp4cloud, phmmer_exe=self.phmmer_exe)
        if not self.bfvd_model_hits:
            logger.critical('SearchModelFinder Big Fantastic Virus Database search could not find any hits!')
            return None
        self.bfvd_model_regions = RegionFinder().find_regions_from_hits(self.bfvd_model_hits)
        return self.bfvd_model_regions
    
    def find_esm_model_regions(self):
        self.esm_model_hits = mr_hit.find_hits(self.seq_info, search_engine="phmmer", hhsearch_exe=None,
                                              hhsearch_db=None, afdb_seqdb=self.afdb_seqdb, bfvd_seqdb=self.bfvd_seqdb,
                                              esm_seqdb=self.esm_seqdb, pdb_seqdb=self.pdb_seqdb,
                                              phmmer_dblvl="esmfold", max_hits=self.max_hits,
                                              nproc=self.nproc, ccp4cloud=self.ccp4cloud, phmmer_exe=self.phmmer_exe)
        if not self.esm_model_hits:
            logger.critical('SearchModelFinder ESMfold Atlas database search could not find any hits!')
            return None
        self.esm_model_regions = RegionFinder().find_regions_from_hits(self.esm_model_hits)
        return self.esm_model_regions

    def prepare_homologs(self):
        if not self.hits and self.regions:
            return None
        self.homologs = mr_homolog.homologs_from_hits(self.hits, self.pdb_dir, self.pdb_local)
        if self.hkl_info:
            mr_homolog.calculate_ellg(self.homologs, self.hkl_info)
        return self.homologs

    def prepare_af2_models(self):
        if not self.af_model_hits and self.af_model_regions:
            return None
        self.af_models = mr_alphafold.models_from_hits(self.af_model_hits, self.plddt_cutoff)
        return self.af_models
    
    def prepare_bfvd_models(self):
        if not self.bfvd_model_hits and self.bfvd_model_regions:
            return None
        self.bfvd_models = mr_bfvd.models_from_hits(self.bfvd_model_hits, self.plddt_cutoff)
        return self.bfvd_models

    def prepare_esm_model(self):
        self.esm_models = mr_esmatlas.models_from_hits(self.esm_model_hits, self.plddt_cutoff)
        return self.esm_models

    def homologs_as_dicts(self):
        """Return a list of per homlog dictionaries serializable to JSON"""
        if not (self.regions and len(self.regions)):
            raise RuntimeError("No regions generated by SearchModelFinder")
        return [h.static_dict for h in self.homologs.values()]

    def af_models_as_dicts(self):
        """Return a list of per model dictionaries serializable to JSON"""
        if not (self.af_model_regions and len(self.af_model_regions)):
            raise RuntimeError("No regions generated by SearchModelFinder")
        return sorted([m.static_dict for m in self.af_models.values()], key=lambda k: k['sum_plddt'], reverse=True)[:20]
    
    def bfvd_models_as_dicts(self):
        """Return a list of per model dictionaries serializable to JSON"""
        if not (self.bfvd_model_regions and len(self.bfvd_model_regions)):
            raise RuntimeError("No regions generated by SearchModelFinder")
        return sorted([m.static_dict for m in self.bfvd_models.values()], key=lambda k: k['sum_plddt'], reverse=True)[:20]

    def esm_models_as_dicts(self):
        return [m.static_dict for m in self.esm_models.values()]

    def homologs_with_graphics(self):
        """List of homologs including PFAM graphics directives
        
        This needs to be done better - the PFAM graphics shouldn't be stored in the 
        list of homologs - this was just done because it made development quicker.
        The list of homologs and PFAM graphics needs to be kept separate
        """
        if not (self.regions and len(self.regions)):
            raise RuntimeError("No regions generated by SearchModelFinder")
        mr_pfam.add_pfam_dict_to_homologs(self.homologs, self.seq_info.nresidues)
        return self.homologs_as_dicts()

    def af_models_with_graphics(self):
        """List of models including PFAM graphics directives

        This needs to be done better - the PFAM graphics shouldn't be stored in the
        list of models - this was just done because it made development quicker.
        The list of models and PFAM graphics needs to be kept separate
        """
        if not (self.af_model_regions and len(self.af_model_regions)):
            raise RuntimeError("No regions generated by SearchModelFinder")
        mr_pfam.add_pfam_dict_to_models(self.af_models, self.seq_info.nresidues, database="EBI AlphaFold database")
        return self.af_models_as_dicts()
    
    def bfvd_models_with_graphics(self):
        """List of models including PFAM graphics directives

        This needs to be done better - the PFAM graphics shouldn't be stored in the
        list of models - this was just done because it made development quicker.
        The list of models and PFAM graphics needs to be kept separate
        """
        if not self.bfvd_model_regions:
            raise RuntimeError("No regions generated by SearchModelFinder")
        mr_pfam.add_pfam_dict_to_models(self.bfvd_models, self.seq_info.nresidues, database='BFVD database')
        return self.bfvd_models_as_dicts()

    def esm_models_with_graphics(self):
        """List of models including PFAM graphics directives

        This needs to be done better - the PFAM graphics shouldn't be stored in the
        list of models - this was just done because it made development quicker.
        The list of models and PFAM graphics needs to be kept separate
        """
        if not self.esm_models:
            raise RuntimeError("No regions generated by SearchModelFinder")
        mr_pfam.add_pfam_dict_to_models(self.esm_models, self.seq_info.nresidues, database='ESMfold Atlas database')
        return self.esm_models_as_dicts()
