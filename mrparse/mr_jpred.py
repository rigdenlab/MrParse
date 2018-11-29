'''
Created on 14 Nov 2018

@author: jmht
'''
import os

from mrparse.mr_annotation import AnnotationSymbol, SequenceAnnotation

HELIX = AnnotationSymbol()
HELIX.symbol = 'H'
HELIX.stype = 'Alpha Helix'
HELIX.name = 'helix'

SHEET = AnnotationSymbol()
SHEET.symbol = 'E'
SHEET.stype = 'B-sheet'
SHEET.name = 'b-sheet'


def parse_jpred_output(jpred_rundir):
    if not os.path.isdir(jpred_rundir):
        raise RuntimeError("Cannot find directory:%s" % jpred_rundir)
    out_concise = [f for f in os.listdir(jpred_rundir) if f.endswith('.concise')][0]
    out_concise = os.path.join(jpred_rundir, out_concise)
    cc_28 = None
    ss_pred = None
    with open(out_concise) as f:
        line = f.readline()
        while line:
            prefix = 'Lupas_28:'
            if line.startswith(prefix):
                line = line.strip().replace(prefix,'')
                cc_28 = "".join(line.split(","))
            prefix = 'jnetpred:'
            if line.startswith(prefix):
                line = line.strip().replace(prefix,'')
                ss_pred = "".join(line.split(","))
            line = f.readline()
    assert cc_28 and ss_pred
    return ss_pred, cc_28


def create_annotation(annotation):
    ann = SequenceAnnotation()
    ann.source = 'Jpred server'
    ann.annotation = annotation
    ann.library_add_annotation(HELIX)
    ann.library_add_annotation(SHEET)
    ann.scores = [1.0] * len(annotation)
    return ann


def secondary_structure_prediction(jpred_rundir):
    ss_pred, cc_28 = parse_jpred_output(jpred_rundir)
    return create_annotation(ss_pred)


