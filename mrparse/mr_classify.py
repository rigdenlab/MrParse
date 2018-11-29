'''
Created on 18 Oct 2018

@author: jmht
'''
from mrparse.mr_deepcoil import coiled_coil_prediction
from mrparse.mr_topcons import Topcons


def get_annotation(seqin, topcons_dir=None):
    cc_pred = coiled_coil_prediction(seqin)
    tm_pred = Topcons().transmembrane_prediction(seqin, topcons_dir=topcons_dir)
    return generate_consensus_classification([cc_pred, tm_pred])


def generate_consensus_classification(annotations):
    lengths = [len(a) for a in annotations]
    assert lengths.count(lengths[0]) == len(lengths), "Annotations have different lengths: %s" % lengths
    for i, a in enumerate(annotations):
        if i == 0:
            consensus = a
            continue
        consensus = consensus + a
    return consensus
