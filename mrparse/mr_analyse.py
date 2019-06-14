import logging
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

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
HTML_DIR = os.path.join(THIS_DIR, '../html')
HTML_TEMPLATE = os.path.join(HTML_DIR, 'mrparse.html.jinja2')
JS_OUT = 'mrparse.js'
HTML_OUT = 'mrparse.html'

logger = None


def run(seqin, hklin=None, run_serial=False, do_classify=True, pdb_download_dir=None):
    # Need to make a work directory first as all logs go into there
    work_dir = make_workdir()
    os.chdir(work_dir)
    global logger
    logger = setup_logging()
    program_name = os.path.basename(sys.argv[0])
    logger.info("Starting: %s", program_name)
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

    search_model_finder = SearchModelFinder(seq_info, hkl_info=hkl_info, pdb_download_dir=pdb_download_dir)
    classifier = None
    if do_classify:
        classifier = MrClassifier(seq_info=seq_info)

    if run_serial:
        run_analyse_serial(search_model_finder, classifier, hkl_info, do_classify)
    else:
        search_model_finder, classifier, hkl_info = run_analyse_parallel(search_model_finder,
                                                                         classifier,
                                                                         hkl_info,
                                                                         do_classify)

    results_json = get_results_json(search_model_finder, hkl_info=hkl_info, classifier=classifier)
    html_out = write_output_files(results_json)
    logger.info("Wrote MrParse output file: %s", html_out)

    # Display results in browser
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


def get_results_json(search_model_finder, hkl_info=None, classifier=None):
    pfam_dict = {}
    pfam_dict.update(search_model_finder.pfam_dict())
    if classifier:
        pfam_dict.update(classifier.pfam_dict())
    data_dict = {'pfam': pfam_dict}
    if hkl_info:
        data_dict['hkl_info'] = hkl_info.as_dict()
    return json.dumps(data_dict)


def write_output_files(results_json):
    js_data = 'const mrparse_data = %s;\n' % results_json
    js_out = os.path.abspath(JS_OUT)
    with open(js_out, 'w') as w:
        w.write(js_data)
    html_out = os.path.abspath(HTML_OUT)
    render_template(HTML_TEMPLATE, html_out, mrparse_html_dir=HTML_DIR)
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