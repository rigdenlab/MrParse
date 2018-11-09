'''
Created on 9 Nov 2018

@author: jmht

Code from: /opt/ccp4/ccp4-7.0/lib/py2/site-packages/phaser/phaser/pipeline/brunett.py



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


with holder( creator = creator ) as manager:
    queue = engine.Adapter( manager = manager, info = info, qslots = qslots )
    state = SYMMETRY_OPTIONS_FOR[ params.symmetry_exploration ](
        composition = mr_object.Composition(
            components = reduce(
                operator.add,
                [
                    [ comp ] * ( stoich * params.composition.count )
                    for ( comp, stoich ) in components
                    ]
                )
            ),
        problem = problem,
        queue = queue,
        output = logger,
        superposition = superposition.Proxy(
            method = superposition.simple_ssm_superposition
            ),
        template_equivalence = params.template_equivalence,
        identity_to_rms = calculate_rmsd_with_phaser,
        )

    # Create solution process control
    ( level1, level2, postprocessors ) = MODE_SETUP_FOR[ params.mode ]( params = params )

    plan = phaser_ai.Plan(
        level1 = level1,
        level2 = level2,
        postprocessors = postprocessors,
        criterion = lambda peak: params.significant_peak_threshold <= peak.tfz,
        b_factor_refinement = params.b_factor_refinement,
        threshold_postrefinement = params.post_refinement_cutoff,
        )

    # Solve
    state.solve( plan = plan )
    manager.shutdown()
    manager.join()
    
# Main code lives in mr_object.py

    
'''
import gzip
import pickle
import sys


from phaser.pipeline import mr_object
from phaser import output

class CrystalSymmetry():
    def __init__(self):
        self.space_group = None
        self.unit_cell = None
 
 
hklin = '/opt/MrParse/data/2uvo_pdbredo.mtz'
labin = 'FP,SIGFP'
crystal_symmetry = CrystalSymmetry()
resolution_cutoff = None
no_tncs = True

load = False
if load:
    logger = output.SingleStream(stream=sys.stdout, level=0, gui=False)
    problem = mr_object.Problem(
        hklin = hklin,
        labin = labin,
        crystal_symmetry = crystal_symmetry,
        resolution = resolution_cutoff,
        no_tncs = no_tncs,
        )
     
    problem.load( logger = logger )
    problem.enpickle('foo')
else:
    with gzip.open('foo_data.pkl.gz') as f:
        problem = pickle.load(f)

print(problem)

