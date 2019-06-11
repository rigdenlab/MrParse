'''
Created on 18 Oct 2018

@author: jmht
'''
import os

from mrbump.ccp4.MRBUMP_ctruncate import Ctruncate
from simbad.util import mtz_util
from simbad.util.matthews_prob import MatthewsProbability
from ample.util.ample_util import filename_append

from Bio import SeqIO
from Bio.SeqUtils import molecular_weight
from Bio.Alphabet import IUPAC

class HklInfo(object):
    def __init__(self, hklin, seqin=None):
        self.hklin = hklin
        self.seqin = seqin
        if not os.path.isfile(hklin):
            raise RuntimeError("Cannot find hklin file: %s" % hklin)
        if seqin and not os.path.isfile(seqin):
            raise RuntimeError("Cannot find seqin file: %s" % seqin)
        self.name = os.path.splitext(os.path.basename(hklin))[0]
        self.labels = mtz_util.GetLabels(hklin)
        self.space_group, self.resolution, self.cell_parameters = mtz_util.crystal_data(self.hklin)
        self.predicted_solvent_content = None
        self.predicted_ncopies = None
        self.molecular_weight = None
        self.has_ncs = False
        self.has_twinning = False
        self.has_anisotropy = False
        if self.seqin:
            self.calculate_matthews_probabilties()
    
    def __call__(self):
        """Required so that we can use multiprocessing pool. We need to be able to pickle the object passed
        to the pool and instance methods don't work, so we add the object to the pool and define __call__
        https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map/6975654#6975654
        """
        self.check_pathologies()
        return self
    
    def calculate_matthews_probabilties(self):
        mp = MatthewsProbability(self.cell_parameters, self.space_group)
        seq_record = SeqIO.read(self.seqin, "fasta", alphabet=IUPAC.protein)
        self.molecular_weight = molecular_weight(seq_record.seq, 'protein')
        nres = len(seq_record)
        self.predicted_solvent_content, self.predicted_ncopies = mp._calculate(nres)
    
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
