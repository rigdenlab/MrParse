'''
Created on 18 Oct 2018

@author: jmht
'''
import logging
import threading


from mrparse.mr_deepcoil import CCPred
from mrparse.mr_topcons import TMPred
from mrparse.mr_jpred import JPred
from mrparse.mr_pfam import pfam_dict_from_annotation


logger = logging.getLogger(__name__)

class PredictorThread(threading.Thread):
    def __init__(self, classifier):
        super(PredictorThread, self).__init__()
        self.classifier = classifier
        self.exception = None
    
    def run(self):
        try:
            self.classifier.get_prediction()
        except Exception as e:
            self.exception = e


class MrClassifier(object):
    def __init__(self, seqin, topcons_dir=None):
        self.seqin = seqin
        self.topcons_dir = topcons_dir
        self.ss_prediction = None
        self.classification_prediction = None
    
    def __call__(self):
        """Required so that we can use multiprocessing pool. We need to be able to pickle the object passed
        to the pool and instance methods don't work, so we add the object to the pool and define __call__
        https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map/6975654#6975654
        """
        self.get_prediction()
        return self
        
    def get_prediction(self):
        cc_predictor = CCPred(self.seqin)
        tm_predictor = TMPred(self.seqin, topcons_dir=self.topcons_dir)
        ss_predictor = JPred(seqin=self.seqin)
        cc_thread = PredictorThread(cc_predictor)
        tm_thread = PredictorThread(tm_predictor)
        ss_thread = PredictorThread(ss_predictor)
        cc_thread.start()
        tm_thread.start()
        ss_thread.start()
        cc_thread.join()
        logger.info('Coiled-Coil predictor finished')
        tm_thread.join()
        logger.info('TM predictor finished')
        ss_thread.join()
        logger.info('SS predictor finished')
        
        if cc_thread.exception and tm_thread.exception:
            logger.critical("BOTH predictors raised expeptions: %s | %s" % (cc_thread.exception, tm_thread.exception))
        elif cc_thread.exception:
            logger.critical("Coiled-Coil predictor raised an exception: %s" % cc_thread.exception)
            self.classification_prediction = tm_predictor.prediction
        elif tm_thread.exception:
            logger.critical("Transmembrane predictor raised an exception: %s" % tm_thread.exception)
            self.classification_prediction = cc_predictor.prediction
        else:
            self.classification_prediction = self.generate_consensus_classification([cc_predictor.prediction, tm_predictor.prediction])
        
        if ss_thread.exception:
            logger.critical("JPred predictor raised error: %s" % ss_thread.exception)
        self.ss_prediction = ss_predictor.prediction
    
    def generate_consensus_classification(self, annotations):
        lengths = [len(a) for a in annotations]
        assert lengths.count(lengths[0]) == len(lengths), "Annotations have different lengths: %s" % lengths
        for i, a in enumerate(annotations):
            if i == 0:
                consensus = a
                continue
            consensus = consensus + a
        return consensus
    
    def pfam_dict(self):
        d = {}
        if self.classification_prediction:
            d['classification'] = pfam_dict_from_annotation(self.classification_prediction)
        if self.ss_prediction:
            d['ss_pred'] = pfam_dict_from_annotation(self.ss_prediction)
        return d