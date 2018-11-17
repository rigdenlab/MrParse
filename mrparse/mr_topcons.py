'''
Created on 16 Nov 2018

@author: jmht
'''
import os


from mrparse.mr_sequence import TM_SYMBOL

def parse_topcons_output(topcons_dir):
    assert os.path.isdir(topcons_dir), "Cannot find directory: %s" % topcons_dir
    results_file = os.path.join(topcons_dir, 'query.result.txt')
    with open(results_file) as fh:
        line = fh.readline()
        while line:
            if line.startswith('TOPCONS predicted topology:'):
                prediction = fh.readline().strip()
            if line.startswith('Predicted TOPCONS reliability'):
                fh.readline()
                line = fh.readline().strip()
                scores = []
                while line:
                    try:
                        seqid, score = line.split()
                    except ValueError:
                        break
                    scores.append((int(seqid), float(score)))
                    line = fh.readline().strip()
            line = fh.readline()
    return prediction, scores

def XXcalculate_probabilities(prediction, scores):
    DEFAULT_PROBABILITY = 50.0
    seqid, score = scores.pop(0)
    probabilities = []
    for i, pred in enumerate(prediction):
        if pred == TM_SYMBOL:
            if seqid == i + 1:  # Assum seqids start from 1
                prob = score
            else:
                prob = DEFAULT_PROBABILITY
        else:
            prob = 0.0
        probabilities.append(prob)
        if seqid == i + 1:
            try:
                seqid, score = scores.pop(0)
            except IndexError:
                # End of list
                pass
    return probabilities


def transmembrane_prediction(topcons_dir):
    prediction, scores = parse_topcons_output(topcons_dir)
    return [1.0 if p == TM_SYMBOL else 0.0 for p in prediction]
