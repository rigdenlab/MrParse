"""
Created on 18 Oct 2018

@author: jmht
"""
import logging
from pathlib import Path

from ample.util.ample_util import filename_append
from cctbx.crystal import symmetry
from mmtbx.scaling.matthews import matthews_rupp
from mrbump.ccp4.MRBUMP_ctruncate import Ctruncate
from simbad.parsers import mtz_parser

logger = logging.getLogger(__name__)


class HklInfo(object):
    def __init__(self, hklin, seq_info=None):
        self.hklin = hklin
        self.seq_info = seq_info
        if not Path(hklin).exists():
            raise RuntimeError(f"Cannot find hklin file: {hklin}")
        self.name = Path(hklin).stem
        self.input_mtz_obj = mtz_parser.MtzParser(hklin)
        self.input_mtz_obj.parse()
        self.space_group = self.input_mtz_obj.spacegroup_symbol
        self.resolution = self.input_mtz_obj.resolution
        self.cell_parameters = self.input_mtz_obj.cell.parameters
        self.predicted_solvent_content = None
        self.predicted_ncopies = None
        self.molecular_weight = None
        self.has_ncs = False
        self.has_twinning = False
        self.has_anisotropy = False
        if self.seq_info:
            self.molecular_weight = self.seq_info.molecular_weight
            self.calculate_matthews_probabilties()

    def __call__(self):
        """Required so that we can use multiprocessing pool. We need to be able to pickle the object passed
        to the pool and instance methods don't work, so we add the object to the pool and define __call__
        https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map/6975654#6975654
        """
        self.check_pathologies()
        return self

    def calculate_matthews_probabilties(self):
        crystal_symmetry = symmetry(unit_cell=self.cell_parameters, space_group_symbol=self.space_group)
        result = matthews_rupp(crystal_symmetry, n_residues=self.seq_info.nresidues)
        self.predicted_solvent_content = result.solvent_content
        self.predicted_ncopies = result.n_copies

    def check_pathologies(self):
        """Function to check crystal pathologies"""
        hklin = self.hklin
        hklout = filename_append(filename=hklin, directory=Path.cwd(), astr='fixcols')
        ctr_colin = None
        ctr_colin_sig = None
        plus_minus = None
        ctr = Ctruncate()
        ctr.debug = False
        log_file = Path(hklout).with_suffix('.log')
        ctr.setlogfile(str(log_file))
        if self.input_mtz_obj.f:
            input_f = True
        else:
            input_f = False

        if self.input_mtz_obj.f or self.input_mtz_obj.i:
            plus_minus = False
            if self.input_mtz_obj.i:
                ctr_colin = self.input_mtz_obj.i
                ctr_colin_sig = self.input_mtz_obj.sigi
            else:
                ctr_colin = self.input_mtz_obj.f
                ctr_colin_sig = self.input_mtz_obj.sigf
        elif self.input_mtz_obj.i_plus:
            plus_minus = True
            ctr_colin = []
            ctr_colin_sig = []
            ctr_colin.append(self.input_mtz_obj.i_plus)
            ctr_colin.append(self.input_mtz_obj.i_minus)
            ctr_colin_sig.append(self.input_mtz_obj.sigi_plus)
            ctr_colin_sig.append(self.input_mtz_obj.sigi_minus)
        elif self.input_mtz_obj.f_plus:
            plus_minus = True
            ctr_colin = []
            ctr_colin_sig = []
            ctr_colin.append(self.input_mtz_obj.f_plus)
            ctr_colin.append(self.input_mtz_obj.f_minus)
            ctr_colin_sig.append(self.input_mtz_obj.sigf_plus)
            ctr_colin_sig.append(self.input_mtz_obj.sigf_minus)

        if self.input_mtz_obj.i and self.input_mtz_obj.free:
            ctr.ctruncate(hklin, hklout, ctr_colin, ctr_colin_sig, colout="from_MRPARSE",
                          colinFREE=self.input_mtz_obj.free, USEINTEN=True, INPUTF=input_f, PLUSMINUS=plus_minus)
        if self.input_mtz_obj.i and not self.input_mtz_obj.free:
            ctr.ctruncate(hklin, hklout, ctr_colin, ctr_colin_sig, colout="from_MRPARSE",
                          USEINTEN=True, INPUTF=input_f, PLUSMINUS=plus_minus)
        elif self.input_mtz_obj.free:
            ctr.ctruncate(hklin, hklout, ctr_colin, ctr_colin_sig, colout="from_MRPARSE",
                          colinFREE=self.input_mtz_obj.free,
                          USEINTEN=False, PLUSMINUS=plus_minus)
        else:
            ctr.ctruncate(hklin, hklout, ctr_colin, ctr_colin_sig, colout="from_MRPARSE",
                          USEINTEN=False, PLUSMINUS=plus_minus)
        self.has_ncs = ctr.NCS
        self.has_twinning = ctr.TWIN
        self.has_anisotropy = ctr.ANISO

        logger.info('Updating HKLIN to %s', hklout)
        self.hklin = hklout
        self.input_mtz_obj = mtz_parser.MtzParser(hklout)
        self.input_mtz_obj.parse()
        return

    def as_dict(self):
        attrs = ['hklin', 'name', 'space_group', 'resolution', 'cell_parameters', 'has_ncs', 'has_twinning',
                 'has_anisotropy']
        d = {}
        for k, v in self.__dict__.items():
            if k in attrs:
                d[k] = v
        return d

    def as_html(self):
        return """
<table border="1">
  <thead>
    <tr style="text-align: right;">
      <th>name</th>
      <th>Resolution</th>
      <th>Space Group</th>
      <th>Has NCS?</th>
      <th>Has Twinning?</th>
      <th>Has Anisotropy?</th>
      <th>File Path</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>{name}</td>
      <td>{resolution:5.3F}</td>
      <td>{space_group}</td>
      <td>{has_ncs}</td>
      <td>{has_twinning}</td>
      <td>{has_anisotropy}</td>
      <td>{hklin}</td>
    </tr>
  </tbody>
</table>
""".format(**self.__dict__)

    def __str__(self):
        ostr = f"HKL Info for file {self.hklin}\n"
        ostr += f"Space Group: {self.space_group}\n"
        ostr += f"Resolution: {self.resolution}\n"
        ostr += f"Has NCS?: {self.has_ncs}\n"
        ostr += f"Has Twinning?: {self.has_twinning}\n"
        ostr += f"Has Anisotropy?: {self.has_anisotropy}\n"
        return ostr
