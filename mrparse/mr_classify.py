'''
Created on 18 Oct 2018

@author: jmht
'''
from mrparse.mr_deepcoil import coiled_coil_prediction, CC
from mrparse.mr_topcons import Topcons, TM
from mrparse.mr_annotation import NULL_ANNOTATION


class MrClassifier(object):
    """Generate a consensus Transmembrane / Coiled-Coil protein classification"""
    
    def __init__(self):
        self.annotations = []
        pass
        
    def get_predictions(self, seqin, topcons_dir=None):
        cc_pred = coiled_coil_prediction(seqin)
        tm_pred = Topcons().transmembrane_prediction(seqin, topcons_dir=topcons_dir)
        self.annotations = [cc_pred, tm_pred]
        return self.generate_consensus_classification(cc_pred, tm_pred)
    
    @staticmethod
    def generate_consensus_classification(cc, tm):
        assert len(cc) == len(tm), "Different lengths: %d : %d" % (len(cc), len(tm))
        classification = ''
        # For now just use prediction and leave probabiltiies
        for i in range(cc.length):
            if cc.annotation[i] == CC and tm.annotation[i] == TM:
                classification += NULL_ANNOTATION.symbol
            elif cc.annotation[i] == CC:
                classification += cc.annotation[i].symbol
            elif tm.annotation[i] == TM:
                classification += tm.annotation[i].symbol
            else:
                classification += NULL_ANNOTATION.symbol
        return classification
