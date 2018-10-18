import set_mrparse_path

from mrparse.mr_search_model import DomainFinder

dfinder = DomainFinder()

hklin = '../data/2uvoA.fasta'
print dfinder.find_domains(hklin)