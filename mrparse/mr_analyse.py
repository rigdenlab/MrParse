import json
import multiprocessing
import os
from pathlib import Path
import subprocess
import sys

from jinja2 import Environment, FileSystemLoader

from mrparse.mr_log import setup_logging
from mrparse.mr_util import make_workdir, now
from mrparse.mr_hkl import HklInfo
from mrparse.mr_search_model import SearchModelFinder
from mrparse.mr_sequence import Sequence, MultipleSequenceException, merge_multiple_sequences
from mrparse.mr_classify import MrClassifier
from mrparse.mr_version import __version__

THIS_DIR = Path(__file__).parent
HTML_DIR = THIS_DIR.joinpath('html')
HTML_TEMPLATE = HTML_DIR.joinpath('mrparse.html.jinja2')

HTML_OUT = 'mrparse.html'
HOMOLOGS_JS = 'homologs.json'
AF_MODELS_JS = 'af_models.json'
BFVD_MODELS_JS = 'bfvd_models.json'
ESM_MODELS_JS = 'esm_models.json'

logger = None


def run(seqin, **kwargs):
    # Get any input kwargs
    hklin = kwargs.get('hklin', None)
    run_serial = kwargs.get('run_serial', None)
    do_classify = kwargs.get('do_classify', None)
    pdb_dir = kwargs.get('pdb_dir', None)
    pdb_local = kwargs.get('pdb_local', None)
    phmmer_dblvl = kwargs.get('phmmer_dblvl', '95')
    plddt_cutoff = kwargs.get('plddt_cutoff', '70')
    search_engine = kwargs.get('search_engine', 'phmmer')
    deeptmhmm_exe = kwargs.get('deeptmhmm_exe', None)
    deepcoil_exe = kwargs.get('deepcoil_exe', None)
    hhsearch_exe = kwargs.get('hhsearch_exe', None)
    hhsearch_db = kwargs.get('hhsearch_db', None)
    afdb_seqdb = kwargs.get('afdb_seqdb', None)
    bfvd_seqdb = kwargs.get('bfvd_seqdb', None)
    esm_seqdb = kwargs.get('esm_seqdb', None)
    pdb_seqdb = kwargs.get('pdb_seqdb', None)
    ccp4cloud = kwargs.get('ccp4cloud', None)
    use_api = kwargs.get('use_api', None)
    max_hits = kwargs.get('max_hits', 10)
    database = kwargs.get('database', 'all')
    nproc = kwargs.get('nproc', 1)

    # Need to make a work directory first as all logs go into there
    work_dir = make_workdir()
    os.chdir(work_dir)
    global logger
    logger = setup_logging()
    program_name = Path(sys.argv[0]).parent
    logger.info(f"Running: {program_name}")
    logger.info(f"Version: {__version__}")
    logger.info(f"Program started at: {now()}")
    logger.info(f"Running from directory: {work_dir}")

    if not (seqin and Path(seqin).exists()):
        raise RuntimeError(f"Cannot find seqin file: {seqin}")
    logger.info(f"Running with seqin {Path(seqin).resolve()}")

    try:
        seq_info = Sequence(seqin)
        seq_info.write(os.path.join(work_dir, 'input.fasta'))
        seq_info.sequence_file = os.path.join(work_dir, 'input.fasta')
    except MultipleSequenceException:
        logger.info(f"Multiple sequences found seqin: {seqin}\n\nAttempting to merge sequences")
        seq_info = merge_multiple_sequences(seqin)
        logger.info(f"Merged sequence file: {seq_info.sequence_file}")

    hkl_info = None
    if hklin:
        if not os.path.isfile(hklin):
            raise RuntimeError(f"Cannot find hklin file: {hklin}")
        logger.info(f"Running with hklin {Path(hklin).resolve()}")
        hkl_info = HklInfo(hklin, seq_info=seq_info)

    if search_engine == "hhsearch":
        if not hhsearch_exe:
            raise RuntimeError("HHSearch executable needs to be defined with --hhsearch_exe")
        elif not hhsearch_db:
            raise RuntimeError("HHSearch database needs to be defined with --hhsearch_db")

    search_model_finder = SearchModelFinder(seq_info, hkl_info=hkl_info, pdb_dir=pdb_dir, phmmer_dblvl=phmmer_dblvl,
                                            plddt_cutoff=plddt_cutoff, search_engine=search_engine, hhsearch_exe=hhsearch_exe, 
                                            hhsearch_db=hhsearch_db, afdb_seqdb=afdb_seqdb, bfvd_seqdb=bfvd_seqdb, esm_seqdb=esm_seqdb, 
                                            pdb_seqdb=pdb_seqdb, max_hits=max_hits, database=database, nproc=nproc, pdb_local=pdb_local, 
                                            ccp4cloud=ccp4cloud)

    classifier = None
    if do_classify:
        classifier = MrClassifier(seq_info=seq_info, deeptmhmm_exe=deeptmhmm_exe, deepcoil_exe=deepcoil_exe)

    if run_serial:
        run_analyse_serial(search_model_finder, classifier, hkl_info, do_classify)
    else:
        search_model_finder, classifier, hkl_info = run_analyse_parallel(search_model_finder,
                                                                         classifier,
                                                                         hkl_info,
                                                                         do_classify)

    html_out = write_output_files(search_model_finder, seq_info=seq_info, hkl_info=hkl_info, classifier=classifier, ccp4cloud=ccp4cloud, database=database)
    logger.info(f"Wrote MrParse output file: {html_out}")

    if not ccp4cloud:
        opencmd = None
        if sys.platform.lower().startswith('linux'):
            opencmd = 'xdg-open'
        elif sys.platform.lower().startswith('darwin'):
            opencmd = 'open'
        if opencmd:
            subprocess.Popen([opencmd, html_out])
    return 0


def run_analyse_serial(search_model_finder, classifier, hkl_info, do_classify):
    if hkl_info:
        try:
            hkl_info()
        except Exception as e:
            logger.critical(f'HklInfo failed: {e}')
            logger.debug("Traceback is:", exc_info=sys.exc_info())
    try:
        search_model_finder()
    except Exception as e:
        logger.critical(f'SearchModelFinder failed: {e}')
        logger.debug("Traceback is:", exc_info=sys.exc_info())
    if do_classify:
        try:
            classifier()
        except Exception as e:
            logger.critical(f'MrClassifier failed: {e}')
            logger.debug("Traceback is:", exc_info=sys.exc_info())
    

def run_analyse_parallel(search_model_finder, classifier, hkl_info, do_classify):
    nproc = 3 if hkl_info else 2
    logger.info(f"Running on {nproc} processors.")
    pool = multiprocessing.Pool(nproc)
    smf_result = pool.apply_async(search_model_finder)

    if do_classify:
        mrc_result = pool.apply_async(classifier)
    if hkl_info:
        hklin_result = pool.apply_async(hkl_info)
    pool.close()
    logger.debug("Pool waiting")
    pool.join()
    logger.debug("Pool finished")
    if hkl_info:
        try:
            hkl_info = hklin_result.get()
        except Exception as e:
            logger.critical(f'HklInfo failed: {e}')
            logger.debug("Traceback is:", exc_info=sys.exc_info())
    try:
        search_model_finder = smf_result.get()
    except Exception as e:
        logger.critical(f'SearchModelFinder failed: {e}')
        logger.debug("Traceback is:", exc_info=sys.exc_info())
    if do_classify:
        try:
            classifier = mrc_result.get()
        except Exception as e:
            logger.critical(f'MrClassifier failed: {e}')
            logger.debug("Traceback is:", exc_info=sys.exc_info())
    
    return search_model_finder, classifier, hkl_info


def write_output_files(search_model_finder, seq_info=None, hkl_info=None, classifier=None, ccp4cloud=None, database="all"):
    # write out homologs for CCP4cloud
    # This code should be updated to separate the storing of homologs from the PFAM directives

    homologs_pfam = {}
    try:
        homologs = search_model_finder.homologs_as_dicts()
        with open(HOMOLOGS_JS, 'w') as w:
            w.write(json.dumps(homologs))
        homologs_pfam = search_model_finder.homologs_with_graphics()
        if ccp4cloud:
            for homolog in homologs_pfam:
                del homolog['pdb_file']
    except RuntimeError:
        logger.debug('No homologues found')

    af_models_pfam = {}
    try:
        af_models = search_model_finder.af_models_as_dicts()
        af_models_js_out = Path(AF_MODELS_JS).resolve()
        with open(af_models_js_out, 'w') as w:
            w.write(json.dumps(af_models))
        af_models_pfam = search_model_finder.af_models_with_graphics()
        if ccp4cloud:
            for af_model in af_models_pfam:
                del af_model['pdb_file']
    except RuntimeError:
        logger.debug('No models found')

    bfvd_models_pfam = {}
    try:
        bfvd_models = search_model_finder.bfvd_models_as_dicts()
        bfvd_models_js_out = Path(BFVD_MODELS_JS).resolve()
        with open(bfvd_models_js_out, 'w') as w:
            w.write(json.dumps(bfvd_models))
        bfvd_models_pfam = search_model_finder.bfvd_models_with_graphics()
        if ccp4cloud:
            for bfvd_model in bfvd_models_pfam:
                del bfvd_model['pdb_file']
    except RuntimeError:
        logger.debug('No models found')

    esm_models_pfam = {}
    try:
        esm_models = search_model_finder.esm_models_as_dicts()
        esm_models_js_out = Path(ESM_MODELS_JS).resolve()
        with open(esm_models_js_out, 'w') as w:
            w.write(json.dumps(esm_models))
        esm_models_pfam = search_model_finder.esm_models_with_graphics()
        if ccp4cloud:
            for esm_model in esm_models_pfam:
                del esm_model['pdb_file']
    except RuntimeError:
        logger.debug('No models found')

    results_dict = {'pfam': {'sequence': seq_info.sequence, 'homologs': homologs_pfam, 'af_models': af_models_pfam, 'bfvd_models': bfvd_models_pfam, 'esm_models': esm_models_pfam}}
    if classifier:
        results_dict['pfam'].update(classifier.pfam_dict())
    if hkl_info:
        results_dict['hkl_info'] = hkl_info.as_dict()
        if ccp4cloud:
            del results_dict['hkl_info']['hklin']
    results_json = json.dumps(results_dict)

    html_out = Path(HTML_OUT).resolve()
    render_template(HTML_TEMPLATE, html_out,
                    # kwargs appear as variables in the template
                    mrparse_html_dir=HTML_DIR,
                    results_json=results_json,
                    version=__version__)
    return html_out


def render_template(in_file_path, out_file_path, **kwargs):
    """
    Templates the given file with the keyword arguments.

    Parameters
    ----------
    in_file_path : Path
       The path to the template
    out_file_path : Path
       The path to output the templated file
    **kwargs : dict
       Variables to use in templating
    """
    env = Environment(loader=FileSystemLoader(in_file_path.parent), keep_trailing_newline=True)
    template = env.get_template(in_file_path.name)
    output = template.render(**kwargs)
    with open(str(out_file_path), "w") as f:
        f.write(output)
