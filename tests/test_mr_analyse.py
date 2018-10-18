import set_mrparse_path

from mrparse.mr_analyse import HklInfo

hklin = '../data/2uvo_pdbredo.mtz'
hkl_info = HklInfo(hklin)

print hkl_info.to_str()