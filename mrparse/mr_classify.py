'''
Created on 18 Oct 2018

@author: jmht
'''
from mrparse import mr_deepcoil
from mrparse import mr_topcons
from mrparse.mr_annotation import NULL_SYMBOL


class MrClassifier(object):
    """Generate a consensus Transmembrane / Coiled-Coil protein classification"""
    
    def __init__(self):
        self.annotations = []
        pass
        
    def get_predictions(self, seqin, topcons_dir=None):
        cc_pred = mr_deepcoil.coiled_coil_prediction(seqin)
        tm_pred = mr_topcons.Topcons().transmembrane_prediction(seqin, topcons_dir=topcons_dir)
        self.annotations = [cc_pred, tm_pred]
        return self.generate_consensus_classification(cc_pred, tm_pred)
    
    @staticmethod
    def generate_consensus_classification(cc, tm):
        assert cc.length == tm.length, "Different lengths: %d : %d" % (cc.length, tm.length)
        classification = ''
        # For now just use prediction and leave probabiltiies
        for i in range(cc.length):
            if cc.annotation[i] in cc.symbols and tm.annotation[i] in tm.symbols:
                classification += NULL_SYMBOL
            elif cc.annotation[i] in cc.symbols:
                classification += cc.annotation[i]
            elif tm.annotation[i] in tm.symbols:
                classification += tm.annotation[i]
            else:
                classification += NULL_SYMBOL
        return classification
