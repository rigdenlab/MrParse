#!/usr/bin/env ccp4-python
import set_mrparse_path
import pytest
import gemmi
import logging
from mrparse.mr_alphafold import download_model, models_from_hits, calculate_quality_threshold, \
    calculate_quality_h_score, calculate_avg_plddt, calculate_sum_plddt, get_afdb_version
from simbad.util.pdb_util import PdbStructure

logging.basicConfig(level=logging.DEBUG)

@pytest.mark.skip(reason="Broken test, method in flux so will fix later")
def test_2uvoA_models(get_2uvo_alphafold_test_hits):
    hits = get_2uvo_alphafold_test_hits
    models = models_from_hits(hits, plddt_cutoff=0)

    m = models['Q0JF21_1']
    assert m.length == 166
    assert abs(m.molecular_weight - 17033) < 0.001
    assert m.model_id == "Q0JF21_1"
    assert m.h_score == 92
    assert int(m.avg_plddt) == 95
    assert int(m.sum_plddt) == 16094

@pytest.mark.skip(reason="Broken test, method in flux so will fix later")
def test_calculate_quality_scores(get_2uvo_alphafold_test_hits):
    """Test the quality scoring methods"""
    hits = get_2uvo_alphafold_test_hits
    hit = hits['Q0JF21_1']

    database_version = get_afdb_version()
    pdb_name = f"AF-{hit.pdb_id}-F1-model_{database_version}.pdb"
    pdb_struct = PdbStructure()
    pdb_string = download_model(pdb_name)
    pdb_struct.structure = gemmi.read_pdb_string(pdb_string)
    seqid_range = range(hit.hit_start, hit.hit_stop + 1)
    pdb_struct.select_residues(to_keep_idx=seqid_range)

    threshold = calculate_quality_threshold(pdb_struct.structure)
    threshold_2 = calculate_quality_threshold(pdb_struct.structure, plddt_threshold=95)
    h_score = calculate_quality_h_score(pdb_struct.structure)
    avg_plddt = calculate_avg_plddt(pdb_struct.structure)
    sum_plddt = calculate_sum_plddt(pdb_struct.structure)
    assert threshold == 100
    assert int(threshold_2) == 70
    assert h_score == 92
    assert int(avg_plddt) == 95
    assert int(sum_plddt) == 16094


if __name__ == '__main__':
    import sys

    pytest.main([__file__] + sys.argv[1:])
