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
    def __init__(self):
        pass
        
    def get_classification(self, seqin, topcons_dir=None):
        cc_predictor = CCPred(seqin)
        tm_predictor = TMPred(seqin, topcons_dir=topcons_dir)
        cc_thread = PredictorThread(cc_predictor)
        tm_thread = PredictorThread(tm_predictor)
        cc_thread.start()
        tm_thread.start()
        cc_thread.join()
        tm_thread.join()
        
        if cc_thread.exception and tm_thread.exception:
            logger.critical("BOTH predictors raised expcetions: %s | %s" % (cc_thread.exception, tm_thread.exception))
        elif cc_thread.exception:
            logger.critical("Coiled-Coil predictor raised an exception: %s" % cc_thread.exception)
            return tm_predictor.prediction
        elif tm_thread.exception:
            logger.critical("Transmembrane predictor raised an exception: %s" % tm_thread.exception)
            return cc_predictor.prediction
        else:
            return self.generate_consensus_classification([cc_predictor.prediction, tm_predictor.prediction])
    
    def generate_consensus_classification(self, annotations):
        lengths = [len(a) for a in annotations]
        assert lengths.count(lengths[0]) == len(lengths), "Annotations have different lengths: %s" % lengths
        for i, a in enumerate(annotations):
            if i == 0:
                consensus = a
                continue
            consensus = consensus + a
        return consensus
    
    def get_annotation_pfam_dict(self, seqin):
        class_a = self.get_classification(seqin)
        class_d = pfam_dict_from_annotation(class_a)
          
        ss_pred = JPred().secondary_structure_prediction(seqin=seqin)
        sspred_d = pfam_dict_from_annotation(ss_pred)
    
        return { 'ss_pred' : sspred_d,
               'classification' :class_d}