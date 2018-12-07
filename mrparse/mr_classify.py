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
    
    def run(self):
        self.classifier.get_prediction()


class MrClassifier(object):
    def __init__(self):
        self.cc_pred = None
        self.tm_pred = None
        self.ss_pred = None
        
    def get_classification(self, seqin, topcons_dir=None):
        self.cc_pred = CCPred(seqin)
        self.tm_pred = TMPred(seqin, topcons_dir=topcons_dir)
        cc_thread = PredictorThread(self.cc_pred)
        tm_thread = PredictorThread(self.tm_pred)
        cc_thread.start()
        tm_thread.start()
        cc_thread.join()
        tm_thread.join()
        return self.generate_consensus_classification([self.cc_pred.prediction, self.tm_pred.prediction])
    
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