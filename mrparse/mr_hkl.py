"""
Created on 18 Oct 2018

@author: jmht
"""
import os

from ample.util.ample_util import filename_append
from cctbx.crystal import symmetry
from mmtbx.scaling.matthews import matthews_rupp
from mrbump.ccp4.MRBUMP_ctruncate import Ctruncate

from iotbx import reflection_file_reader
from iotbx.reflection_file_utils import looks_like_r_free_flags_info


class GetLabels(object):
    """Class to get the column labels for input mtz file
    Attributes
    ----------
    f : str
        f column label
    sigf : str
        fp column label
    i : str
        i column label
    sigi : str
        sigi column label
    fplus : str
        f(+) column label
    sigfplus : str
        sigf(+) column label
    fminus : str
        f(-) column label
    sigfminus : str
        sigf(-) column label
    iplus : str
        i(+) column label
    sigiplus : str
        sigi{+} column label
    iminus : str
        i(-) column label
    sigiminus : str
        sigi(-) column label
    dano : str
        dano column label
    sigdano : str
        sigdano column label
    free : str
        free column label
    """

    def __init__(self, mtz_file):
        self.f = None
        self.sigf = None
        self.i = None
        self.sigi = None
        self.fplus = None
        self.sigfplus = None
        self.fminus = None
        self.sigfminus = None
        self.iplus = None
        self.sigiplus = None
        self.iminus = None
        self.sigiminus = None
        self.dano = None
        self.sigdano = None
        self.free = None

        self.run(mtz_file)

    def run(self, mtz_file):
        reflection_file = reflection_file_reader.any_reflection_file(file_name=mtz_file)
        if not reflection_file.file_type() == "ccp4_mtz":
            msg = "File is not of type ccp4_mtz: {0}".format(mtz_file)
            LOG.critical(msg)
            raise RuntimeError(msg)

        miller_arrays = reflection_file.as_miller_arrays()

        for m_a in miller_arrays:
            if looks_like_r_free_flags_info(m_a.info()) and not self.free:
                self.free = m_a.info().labels[0]
            elif self.check_anomalous(m_a):
                if self.check_for_dano_labels(m_a):
                    if len(m_a.info().labels) == 5:
                        self.f, self.sigf, self.dano, self.sigdano, isym = m_a.info().labels
                    elif len(m_a.info().labels) == 4:
                        self.f, self.sigf, self.dano, self.sigdano = m_a.info().labels
                    elif len(m_a.info().labels) == 2:
                        self.dano, self.sigdano = m_a.info().labels
                    else:
                        LOG.debug("Unexpected number of columns found in anomalous miller array")
                elif self.check_for_plus_minus_labels(m_a):
                    if m_a.is_xray_amplitude_array():
                        self.fplus, self.sigfplus, self.fminus, self.sigfminus = m_a.info().labels
                    elif m_a.is_xray_intensity_array():
                        self.iplus, self.sigiplus, self.iminus, self.sigiminus = m_a.info().labels
                    else:
                        LOG.debug("Type of anomalous miller array unknown")
                else:
                    LOG.debug("Type of anomalous miller array unknown")
            elif m_a.is_xray_intensity_array() and len(m_a.info().labels) == 2 and not self.i:
                self.i, self.sigi = m_a.info().labels
            elif m_a.is_xray_amplitude_array() and len(m_a.info().labels) == 2 and not self.f:
                self.f, self.sigf = m_a.info().labels

    def check_anomalous(self, miller_array):
        if miller_array.anomalous_flag():
            return True
        elif miller_array.info().type_hints_from_file == "anomalous_difference":
            return True
        # Check for anomalous miller arrays which aren't properly labeled
        elif self.check_for_dano_labels(miller_array):
            return True
        elif self.check_for_plus_minus_labels(miller_array):
            return True
        return False

    @staticmethod
    def check_for_dano_labels(miller_array):
        return any(["DANO" in i.upper() or "DP" == i.upper() for i in miller_array.info().labels])

    @staticmethod
    def check_for_plus_minus_labels(miller_array):
        return any(["(+)" in i for i in miller_array.info().labels])


def crystal_data(mtz_file):
    """Set crystallographic parameters from mtz file
    Parameters
    ----------
    mtz_file : str
       The path to the mtz file
    Returns
    -------
    space_group : str
       The space group
    resolution : str
       The resolution
    cell_parameters : tuple
       The cell parameters
    """

    reflection_file = reflection_file_reader.any_reflection_file(file_name=mtz_file)
    content = reflection_file.file_content()
    space_group = content.space_group_name().replace(" ", "")
    resolution = content.max_min_resolution()[1]
    cell_parameters = content.crystals()[0].unit_cell_parameters()
    return space_group, resolution, cell_parameters


class HklInfo(object):
    def __init__(self, hklin, seq_info=None):
        self.hklin = hklin
        self.seq_info = seq_info
        if not os.path.isfile(hklin):
            raise RuntimeError("Cannot find hklin file: %s" % hklin)
        self.name = os.path.splitext(os.path.basename(hklin))[0]
        self.labels = GetLabels(hklin)
        self.space_group, self.resolution, self.cell_parameters = crystal_data(self.hklin)
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
        """DOC TODO"""
        hklin = self.hklin
        hklout = filename_append(filename=hklin, directory=os.getcwd(), astr='fixcols')
        ctr_colin = None
        ctr_colin_sig = None
        plus_minus = None
        mtz_obj = self.labels
        ctr = Ctruncate()
        ctr.debug = False
        log_file = hklout.rsplit(".", 1)[0] + '.log'
        ctr.setlogfile(log_file)
        if mtz_obj.f:
            input_f = True
        else:
            input_f = False
    
        if mtz_obj.f or mtz_obj.i:
            plus_minus = False
            if mtz_obj.i:
                ctr_colin = mtz_obj.i
                ctr_colin_sig = mtz_obj.sigi
            else:
                ctr_colin = mtz_obj.f
                ctr_colin_sig = mtz_obj.sigf
        elif mtz_obj.iplus:
            plus_minus = True
            ctr_colin = []
            ctr_colin_sig = []
            ctr_colin.append(mtz_obj.fplus)
            ctr_colin.append(mtz_obj.fminus)
            ctr_colin_sig.append(mtz_obj.sigfplus)
            ctr_colin_sig.append(mtz_obj.sigfminus)
        elif mtz_obj.fplus:
            plus_minus = True
            ctr_colin = []
            ctr_colin_sig = []
            ctr_colin.append(mtz_obj.fplus)
            ctr_colin.append(mtz_obj.fminus)
            ctr_colin_sig.append(mtz_obj.sigfplus)
            ctr_colin_sig.append(mtz_obj.sigfminus)
    
        if mtz_obj.i and mtz_obj.free:
            ctr.ctruncate(hklin, hklout, ctr_colin, ctr_colin_sig, colout="from_SIMBAD", colinFREE=mtz_obj.free,
                          USEINTEN=True, INPUTF=input_f, PLUSMINUS=plus_minus)
        elif mtz_obj.i and not mtz_obj.free:
            ctr.ctruncate(hklin, hklout, ctr_colin, ctr_colin_sig, colout="from_SIMBAD", USEINTEN=True, INPUTF=input_f,
                          PLUSMINUS=plus_minus)
        elif mtz_obj.free:
            ctr.ctruncate(hklin, hklout, ctr_colin, ctr_colin_sig, colout="from_SIMBAD", colinFREE=mtz_obj.free,
                          USEINTEN=False, PLUSMINUS=plus_minus)
        else:
            ctr.ctruncate(hklin, hklout, ctr_colin, ctr_colin_sig, colout="from_SIMBAD", USEINTEN=False,
                          PLUSMINUS=plus_minus)
        self.has_ncs = ctr.NCS
        self.has_twinning = ctr.TWIN
        self.has_anisotropy = ctr.ANISO
        os.unlink(hklout)
        return

    def as_dict(self):
        attrs = ['hklin', 'name', 'space_group', 'resolution', 'cell_parameters', 'has_ncs', 'has_twinning', 'has_anisotropy']
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
        ostr = "HKL Info for file %s\n" % self.hklin
        ostr += "Space Group: %s\n" % self.space_group
        ostr += "Resolution: %s\n" % self.resolution
#         ostr += "Cell Parameters: %s\n" % self.cell_parameters
        ostr += "Has NCS?: %s\n" % self.has_ncs
        ostr += "Has Twinning?: %s\n" % self.has_twinning
        ostr += "Has Anisotropy?: %s\n" % self.has_anisotropy
        return ostr
