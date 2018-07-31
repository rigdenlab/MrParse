#!/usr/bin/env ccp4-python
import os
import sys
from mmtbx.scaling import twin_analyses
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
import iotbx.bioinformatics
from mmtbx.scaling import matthews

hklin = os.path.join(os.environ["CEXAM"],"data", "1vr7_lr_i.mtz")
#hklin = os.path.expandvars("$CCP4/lib/py2/site-packages/phaser/tutorial/beta_blip_P3221.mtz")

xray_data_server =  reflection_file_utils.reflection_file_server(
  crystal_symmetry = None,
  force_symmetry = True,
  reflection_files=[])

miller_obs = xray_data_server.get_xray_data(
  file_name = hklin,
  labels = 'IMEAN,SIGIMEAN',
  ignore_all_zeros = True,
  parameter_scope = 'L1',
  parameter_name='obs_labels')


n_residues = None
# sequence_file = None
# seq_comp = iotbx.bioinformatics.composition_from_sequence_file(file_name=sequence_file,
#                                                                log=text_out)
# if (seq_comp is not None) :
#   n_residues = seq_comp.n_residues

matthews_results = matthews.matthews_rupp(
  crystal_symmetry = miller_obs,
  n_residues = n_residues,
  n_bases = None)

print("GOT ",matthews_results.n_copies)


# reflection_file = reflection_file_reader.any_reflection_file(file_name=hklin)
# mas =  reflection_file.file_content().as_miller_arrays()
# 
# #mas[0].has_twinning()
# # This seems the best way to check for twinning.
# mas[0].analyze_intensity_statistics().has_twinning()
# #laws = twin_analyses.twin_laws(miller_array=mas[0])


