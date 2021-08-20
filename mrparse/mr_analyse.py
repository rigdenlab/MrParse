import json
import multiprocessing
import os
import subprocess
import sys

from jinja2 import Environment, FileSystemLoader

from mr_log import setup_logging
from mr_util import make_workdir, now
from mr_hkl import HklInfo
from mr_search_model import SearchModelFinder
from mr_sequence import Sequence, MultipleSequenceException, merge_multiple_sequences
from mr_classify import MrClassifier
from mr_version import __version__

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
HTML_DIR = os.path.join(THIS_DIR, 'html')
HTML_TEMPLATE = os.path.join(HTML_DIR, 'mrparse.html.jinja2')
HTML_OUT = 'mrparse.html'
HOMOLOGS_JS = 'homologs.json'
MODELS_JS = 'models.json'

logger = None


def run(seqin, hklin=None, run_serial=False, do_classify=True, pdb_dir=None, db_lvl=None, tmhmm_exe=None,
        deepcoil_exe=None, ccp4cloud=None):
    # Need to make a work directory first as all logs go into there
    work_dir = make_workdir()
    os.chdir(work_dir)
    global logger
    logger = setup_logging()
    program_name = os.path.basename(sys.argv[0])
    logger.info("Running: %s", program_name)
    logger.info("Version: %s", __version__)
    logger.info("Program started at: %s", now())
    logger.info("Running from directory: %s", work_dir)

    if not (seqin and os.path.isfile(seqin)):
        raise RuntimeError("Cannot find seqin file: %s" % seqin)
    logger.info("Running with seqin %s", os.path.abspath(seqin))

    try:
        seq_info = Sequence(seqin)
    except MultipleSequenceException:
        logger.info("Multiple sequences found seqin: %s\n\nAttempting to merge sequences", seqin)
        seq_info = merge_multiple_sequences(seqin)
        logger.info("Merged sequence file: %s", seq_info.sequence_file)

    hkl_info = None
    if hklin:
        if not os.path.isfile(hklin):
            raise RuntimeError("Cannot find hklin file: %s" % hklin)
        logger.info("Running with hklin %s", os.path.abspath(hklin))
        hkl_info = HklInfo(hklin, seq_info=seq_info)

    search_model_finder = SearchModelFinder(seq_info, hkl_info=hkl_info, pdb_dir=pdb_dir, db_lvl=db_lvl)

    classifier = None
    if do_classify:
        classifier = MrClassifier(seq_info=seq_info, tmhmm_exe=tmhmm_exe, deepcoil_exe=deepcoil_exe)

    if run_serial:
        run_analyse_serial(search_model_finder, classifier, hkl_info, do_classify)
    else:
        search_model_finder, classifier, hkl_info = run_analyse_parallel(search_model_finder,
                                                                         classifier,
                                                                         hkl_info,
                                                                         do_classify)

    #     results_json = get_results_json(search_model_finder, hkl_info=hkl_info, classifier=classifier)
    #     html_out = write_output_files(results_json)
    html_out = write_output_files(search_model_finder, hkl_info=hkl_info, classifier=classifier, ccp4cloud=ccp4cloud)
    logger.info("Wrote MrParse output file: %s", html_out)

    # Display results in browser
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
    try:
        search_model_finder()
    except Exception as e:
        logger.critical('SearchModelFinder failed: %s' % e)
        logger.debug("Traceback is:", exc_info=sys.exc_info())
    if do_classify:
        try:
            classifier()
        except Exception as e:
            logger.critical('MrClassifier failed: %s' % e)
            logger.debug("Traceback is:", exc_info=sys.exc_info())
    if hkl_info:
        try:
            hkl_info()
        except Exception as e:
            logger.critical('HklInfo failed: %s' % e)
            logger.debug("Traceback is:", exc_info=sys.exc_info())


def run_analyse_parallel(search_model_finder, classifier, hkl_info, do_classify):
    nproc = 3 if hkl_info else 2
    logger.info("Running on %d processors." % nproc)
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
    try:
        search_model_finder = smf_result.get()
    except Exception as e:
        logger.critical('SearchModelFinder failed: %s' % e)
        logger.debug("Traceback is:", exc_info=sys.exc_info())
    if do_classify:
        try:
            classifier = mrc_result.get()
        except Exception as e:
            logger.critical('MrClassifier failed: %s' % e)
            logger.debug("Traceback is:", exc_info=sys.exc_info())
    if hkl_info:
        try:
            hkl_info = hklin_result.get()
        except Exception as e:
            logger.critical('HklInfo failed: %s' % e)
            logger.debug("Traceback is:", exc_info=sys.exc_info())
    return search_model_finder, classifier, hkl_info


def write_output_files(search_model_finder, hkl_info=None, classifier=None, ccp4cloud=None):
    # write out homologs for CCP4cloud
    # This code should be updated to separate the storing of homologs from the PFAM directives

    homologs_pfam = {}
    try:
        homologs = search_model_finder.homologs_as_dicts()
        homologs_js_out = os.path.abspath(HOMOLOGS_JS)
        with open(homologs_js_out, 'w') as w:
            w.write(json.dumps(homologs))
        homologs_pfam = search_model_finder.homologs_with_graphics()
        if ccp4cloud:
            for homolog in homologs_pfam:
                del homolog['pdb_file']
    except RuntimeError:
        logger.debug('No homologues found')

    models_pfam = {}
    try:
        models = search_model_finder.models_as_dicts()
        models_js_out = os.path.abspath(MODELS_JS)
        with open(models_js_out, 'w') as w:
            w.write(json.dumps(models))
        models_pfam = search_model_finder.models_with_graphics()
        if ccp4cloud:
            for model in models_pfam:
                del model['pdb_file']
    except RuntimeError:
        logger.debug('No models found')

    results_dict = {'pfam': {'homologs': homologs_pfam, 'models': models_pfam}}
    if classifier:
        results_dict['pfam'].update(classifier.pfam_dict())
    if hkl_info:
        results_dict['hkl_info'] = hkl_info.as_dict()
    results_json = json.dumps(results_dict)

    html_out = os.path.abspath(HTML_OUT)
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
    in_file_path : str
       The path to the template
    out_file_path : str
       The path to output the templated file
    **kwargs : dict
       Variables to use in templating
    """
    env = Environment(loader=FileSystemLoader(os.path.dirname(in_file_path)), keep_trailing_newline=True)
    template = env.get_template(os.path.basename(in_file_path))
    output = template.render(**kwargs)
    with open(out_file_path, "w") as f:
        f.write(output)
