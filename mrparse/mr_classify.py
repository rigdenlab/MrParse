"""
Created on 18 Oct 2018

@author: jmht
"""
import logging
import threading
import sys

from mrparse.mr_deepcoil import CCPred
from mrparse.mr_deeptmhmm import TMPred
from mrparse.mr_jpred import JPred
from mrparse.mr_pfam import pfam_dict_from_annotation

logger = logging.getLogger(__name__)


class PredictorThread(threading.Thread):
    def __init__(self, classifier):
        super(PredictorThread, self).__init__()
        self.classifier = classifier
        self.exc_info = None
        self.exception = None

    def run(self):
        try:
            self.classifier.get_prediction()
        except Exception as e:
            self.exc_info = sys.exc_info()
            self.exception = e


class MrClassifier(object):
    def __init__(self, seq_info, do_ss_predictor=True, do_cc_predictor=True, do_tm_predictor=True, deeptmhmm_exe=None,
                 deepcoil_exe=None):
        self.seq_info = seq_info
        self.do_ss_predictor = do_ss_predictor
        self.do_cc_predictor = do_cc_predictor
        self.do_tm_predictor = do_tm_predictor
        self.deeptmhmm_exe = deeptmhmm_exe
        self.deepcoil_exe = deepcoil_exe
        self.ss_prediction = None
        self.classification_prediction = None

    def __call__(self):
        """Required so that we can use multiprocessing pool. We need to be able to pickle the object passed
        to the pool and instance methods don't work, so we add the object to the pool and define __call__
        https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map/6975654#6975654
        """
        self.get_prediction()
        return self

    @staticmethod
    def generate_consensus_classification(annotations):
        consensus = None
        lengths = [len(a) for a in annotations]
        assert lengths.count(lengths[0]) == len(lengths), f"Annotations have different lengths: {lengths}"
        for i, a in enumerate(annotations):
            if i == 0:
                consensus = a
                continue
            consensus = consensus + a
        return consensus

    def get_prediction(self):
        cc_predictor = None
        tm_predictor = None
        ss_predictor = None

        if self.do_cc_predictor:
            cc_predictor = CCPred(self.seq_info, self.deepcoil_exe)
            cc_thread = PredictorThread(cc_predictor)
            cc_thread.start()
            cc_thread.join()
            if cc_thread.exception:
                logger.warning(f"Coiled-Coil predictor raised an exception: {cc_thread.exception}")
                logger.debug("Traceback is:", exc_info=cc_thread.exc_info)
                self.do_cc_predictor = False
            logger.info('Coiled-Coil predictor finished')

        if self.do_tm_predictor:
            tm_predictor = TMPred(self.seq_info, self.deeptmhmm_exe)
            tm_thread = PredictorThread(tm_predictor)
            tm_thread.start()
            tm_thread.join()
            if tm_thread.exception:
                logger.warning(f"Transmembrane predictor raised an exception: {tm_thread.exception}")
                logger.debug("Traceback is:", exc_info=tm_thread.exc_info)
                self.do_tm_predictor = False
            logger.info('TM predictor finished')
        if self.do_ss_predictor:
            ss_predictor = JPred(seq_info=self.seq_info)
            ss_thread = PredictorThread(ss_predictor)
            ss_thread.start()
            ss_thread.join()
            if ss_thread.exception:
                logger.warning(f"JPred predictor raised error: {ss_thread.exception}")
                logger.debug("Traceback is:", exc_info=ss_thread.exc_info)
                self.do_ss_predictor = False
            logger.info('SS predictor finished')

        # Determine pediction
        if self.do_cc_predictor and cc_predictor and self.do_tm_predictor and tm_predictor:
            self.classification_prediction = self.generate_consensus_classification(
                [cc_predictor.prediction, tm_predictor.prediction])
        elif self.do_cc_predictor and cc_predictor:
            self.classification_prediction = cc_predictor.prediction
        elif self.do_tm_predictor and tm_predictor:
            self.classification_prediction = tm_predictor.prediction

        if self.do_ss_predictor and ss_predictor:
            self.ss_prediction = ss_predictor.prediction

    def pfam_dict(self):
        d = {}
        if self.classification_prediction:
            d['classification'] = pfam_dict_from_annotation(self.classification_prediction)
        if self.ss_prediction:
            d['ss_pred'] = pfam_dict_from_annotation(self.ss_prediction)
        return d
