'''
Created on 18 Oct 2018

@author: jmht

Run jpred
Run deepcoil
Run topcons

Consolidate results and generate prediction:
* get single array with prediction (ideally as probability from each predictor)
* merge into an array with the final prediction
* generate the tuples of (start, stop, type)
* generate pfam json


'''
from mrparse import mr_deepcoil
from mrparse import mr_topcons
from mrparse import mr_jpred
from mrparse import mr_pfam
from mrparse.mr_sequence import read_fasta, SequenceChunk, get_sequence_chunks, \
     UNKNOWN_SYMBOL, TM_SYMBOL, CC_SYMBOL, HELIX_SYMBOL, BSHEET_SYMBOL

    
def SequenceClassification(object):
    def __init__(self):
        self.source = None
        self.symbol = None


class MrClassifier(object):
    
    def __init__(self):
        self.classification = None
        self.sspred = None
        self.seqlen = None
        
    def get_predictions(self, seqin, topcons_dir=None, jpred_dir=None):
        assert topcons_dir and jpred_dir
        self.seqlen =len(read_fasta(seqin))
        cc_pred = mr_deepcoil.coiled_coil_prediction(seqin)
        tm_pred = mr_topcons.Topcons().transmembrane_prediction(seqin)
        self.classification = self.generate_consensus_classification(cc_pred, tm_pred)
        self.sspred = mr_jpred.secondary_structure_prediction(jpred_dir)
        return
    
    @staticmethod
    def generate_consensus_classification(cc_pred, tm_pred):
        assert len(cc_pred) == len(tm_pred)
        classification = ''
        THRESHOLD = 0.6
        for i, (cc, tm) in enumerate(zip(cc_pred, tm_pred)):
            if cc <= THRESHOLD and tm <= THRESHOLD:
                classification += UNKNOWN_SYMBOL
            else:
                if cc > tm:
                    classification += CC_SYMBOL
                elif tm > cc:
                    classification += TM_SYMBOL
                else:
                    # If they're the same, use whatever came last
                    if i > 0:
                        classification += classification[i-1]
                    else:
                        classification += UNKNOWN_SYMBOL
        return classification
                
    
    def pfam_data(self):
        assert self.classification and self.sspred
        classification_chunks = get_sequence_chunks(self.classification, markers=[CC_SYMBOL, TM_SYMBOL], source='Deepcoil')
        sspred_chunks = get_sequence_chunks(self.sspred, markers=[HELIX_SYMBOL, BSHEET_SYMBOL], source='Jpred')
        classification_dict = mr_pfam.pfam_classification_dict(classification_chunks, self.seqlen)
        sspred_dict = mr_pfam.pfam_classification_dict(sspred_chunks, self.seqlen)
        return {'classification' : classification_dict,
                'ss_pred' : sspred_dict} 
        

# seqin = 'foo'
# classifier = MrClassifier()
# classfier.get_predictions(seqin)
# classifier.pfam_json()