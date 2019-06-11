import set_mrparse_path
import conftest


from mrparse.mr_topcons import TMPred, TM


def test_parse():
    tc = TMPred(None)
    results_dir = '/opt/MrParse/data/Q13586/topcons'
    prediction, scores = tc.parse_topcons_output(results_dir)
    annotation = tc.create_annotation(prediction, scores)

    assert len(annotation) == 685
    assert annotation.annotation[212] == TM.symbol, annotation.annotation[212] 
    assert annotation.annotation[232] == TM.symbol, annotation.annotation[232] 
    assert annotation[212].score == annotation.scores[212]

    assert len(annotation.scores) == 685
    assert annotation.scores[212] > 0.6, annotation.scores[212] 
    assert annotation.scores[232] > 0.6, annotation.scores[232] 


def test_run():
    seqin = '/opt/MrParse/data/Q13586.fasta'
    tc = TMPred(seqin)
    annotation = tc.get_prediction()
    assert len(annotation) == 685
    assert annotation.annotation[212] == TM.symbol, annotation.annotation[212] 
    assert annotation.annotation[232] == TM.symbol, annotation.annotation[232] 
    assert annotation[212].score == annotation.scores[212]

    assert len(annotation.scores) == 685
    assert annotation.scores[212] > 0.6, annotation.scores[212] 
    assert annotation.scores[232] > 0.6, annotation.scores[232] 
