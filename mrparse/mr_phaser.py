'''
Created on 9 Nov 2018

@author: jmht

Code from: /opt/ccp4/ccp4-7.0/lib/py2/site-packages/phaser/phaser/pipeline/brunett.py

See also:
phaser/pipeline/test_jobs.py
phaser/pipeline/phaser_ai.py
phaser/pipeline/calculation.py



def read_any_format_data(hklin, labin, crystal_symmetry, resolution, no_tncs, logger):

see def run

# Create current state
problem = mr_object.Problem(
    hklin = params.hklin,
    labin = params.labin,
    crystal_symmetry = params.crystal_symmetry,
    resolution = params.resolution_cutoff,
    no_tncs = params.no_tncs,
    )

for e in models:
    problem.data.write( obj = e )


a = mr_object.Assembly(
    definitions = definitions,
    symmetry = mr_object.PointGroup(),
    strategy = mr_object.SearchStrategy(
        sweep = assembly.local_search.sweep,
        extent = assembly.local_search.extent,
        ),
    )
problem.data.write( obj = a )

    
# Main code lives in mr_object.py

    
'''
import gzip
import os
import pickle
import sys

import phaser
from phaser import mmt
from phaser.pipeline import mr_object
from phaser import output
from phaser import tbx_utils


class CrystalSymmetry():
    def __init__(self):
        self.space_group = None
        self.unit_cell = None
        

def get_number_of_residues(file_name):

    if not os.path.isfile( file_name ):
        raise RuntimeError("Cannot find PDB file: %s" % file_name)

    import iotbx.pdb
    root = iotbx.pdb.input( file_name ).construct_hierarchy()

    if not root.models():
        raise RuntimeError("PDB file empty: %s" % file_name)

    return len( list( root.models()[0].residue_groups() ) )


 
# hklin = '/opt/MrParse/data/2uvo_pdbredo.mtz'
# labin = 'FP,SIGFP'
# fasta = '/opt/MrParse/data/2uvoA.fasta'
# name = '2uvo'
hklin = '/opt/ample.git/examples/toxd-example/input/toxd.mtz'
labin = 'FTOXD3,SIGFTOXD3'
fasta = '/opt/ample.git/examples/toxd-example/input/toxd_.fasta'
search_pdb = '/opt/ample.git/examples/toxd-example/input/1dtx_model.pdb'
identity = 1.0
name = 'toxd'

stoichiometry = 1 # number of copies of this component in the ASU
crystal_symmetry = CrystalSymmetry()
resolution_cutoff = None
no_tncs = True

read = True
if read:
    logger = output.SingleStream(stream=sys.stdout, level=0, gui=False)
    problem = mr_object.Problem(
        hklin = hklin,
        labin = labin,
        crystal_symmetry = crystal_symmetry,
        resolution = resolution_cutoff,
        no_tncs = no_tncs,
        )
     
    problem.load( logger = logger )
    problem.enpickle(name)
else:
    with gzip.open('%s_data.pkl.gz' % name) as f:
        problem = pickle.load(f)

print(problem)

components = []
models = []


mtype_str = 'protein'
mtype = mmt.by_name( name = mtype_str )
logger.info( msg = "Component type: %s" % mtype.name )

seqfile = tbx_utils.SequenceObject.from_file( file_name = fasta )

if len( seqfile.object ) != 1:
    raise RuntimeError("Sequence file '%s' contain multiple sequences" % fasta)

sequence = mr_object.SequenceData(
    mtype = mtype,
    sequence = seqfile.object[0],
    )

from phaser.pipeline import domain_analysis

sequence_component = mr_object.SequenceComponent(
    sequence = sequence,
    selection = domain_analysis.ChainSequenceSelection.from_selection(
        selection = [ True ] * len( sequence )
        )
    )
logger.info( msg = "Molecular weight for this component: %.1f" % sequence_component.mw )
components.append( ( sequence_component, stoichiometry ) )


num_residues = get_number_of_residues( file_name = search_pdb )
rms = phaser.rms_estimate().rms(identity, num_residues)
rmsds = [rms] 
pdbs = [search_pdb] 


composition = mr_object.Composition( components = [sequence_component] )
model = mr_object.ModelCollection(
    pdbs = pdbs,
    rmsds = rmsds,
    trim = False,
    composition = composition
    )

problem.data.write( obj = model )


logger.info( msg = "Guessing number of copies based on Matthews coefficient" )
from mmtbx.scaling import matthews
count_probabilities = matthews.number_table(
    components = [
        matthews.component(
            mw = comp.mw * stoich,
            rho_spec = mmt.PROTEIN.specific_volume
            )
        for ( comp, stoich ) in components
        ],
    density_calculator = matthews.density_calculator(
        crystal = problem.crystal()
        )
    )
count = max( count_probabilities, key = lambda p: p[1] )[0]


# solution_progress = mr_object.Case(
#                 composition = composition,
#                 space_group_hall = problem.space_group_hall,
#                 problem = problem
#                 )
