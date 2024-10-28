import argparse
import os
from pathlib import Path
import configparser as ConfigParser

from mrparse.mr_version import __version__


class FilePathAction(argparse.Action):
    """Class to handle paths to files or directories.
    
    We cd into a work directory so relative paths to files don't work.
    We set absolulte paths here.
    """

    def __call__(self, parser, namespace, values, option_string=None):
        if isinstance(values, str):
            values = os.path.abspath(values)
        setattr(namespace, self.dest, values)


def mrparse_argparse(parser):
    """Parse MrParse command line arguments"""
    # Read command line arguments
    sg = parser.add_argument_group("Basic options")
    sg.add_argument('-hkl', '--hklin', action=FilePathAction, help='MTZ/CIF Crystal Data file')
    sg.add_argument('--do_classify', action='store_true',
                    help='Run the SS/TM/CC classifiers - requires internet access.')
    sg.add_argument('--pdb_dir', action=FilePathAction, help='Directory of PDB files')
    sg.add_argument('--pdb_local', action=FilePathAction, help='Path to locally stored PDB mirror')
    sg.add_argument('--phmmer_dblvl', help='Redundancy level of PDB database used by Phmmer', default='95',
                    choices=['50', '70', '90', '95', '100'])
    sg.add_argument('--plddt_cutoff', help='Removes residues from AFDB models below this pLDDT threshold',
                    default='70', choices=['50', '70', '90', 'None'])
    sg.add_argument('--run_serial', action='store_true', help='Run on a single processor')
    sg.add_argument('-seq', '--seqin', action=FilePathAction, help='Sequence file')
    sg.add_argument('--search_engine', help="Select search engine", default="phmmer",
                    choices=['phmmer', 'hhsearch'])
    sg.add_argument('--deeptmhmm_exe', action=FilePathAction,
                    help="Location of DeepTMHMM executable for transmembrane classification")
    sg.add_argument('--deepcoil_exe', action=FilePathAction,
                    help="Location of Deepcoil executable for coiled-coil classification")
    sg.add_argument('--hhsearch_exe', action=FilePathAction,
                    help="Location of hhsearch executable")
    sg.add_argument('--hhsearch_db', help="Location of hhsearch database")
    sg.add_argument('--afdb_seqdb', help="Location of alternative alphafold sequence database. To search the entire alphafold database, download the latest sequence listing here: https://ftp.ebi.ac.uk/pub/databases/alphafold/sequences.fasta and set path to file using this option. Note: very large file size!")
    sg.add_argument('--pdb_seqdb', help="Location of alternative PDB sequence database. To search the entire PDB database, download the latest sequence listing here: https://ftp.pdbj.org/pub/pdb/derived_data/pdb_seqres.txt and set path to file using this option.")
    sg.add_argument('--ccp4cloud', action='store_true', help="Specify running through CCP4Cloud")
    sg.add_argument('--use_api', action='store_true', default=True, help='Run alphafold database search using EBI API database search')
    sg.add_argument('--max_hits', required=False, type=int, choices=range(1,101), metavar="[1-100]", default=10, help='Maximum number of models to download and prepare for each database search')
    sg.add_argument('--nproc', required=False, type=int, default=1, help='Number of cores to use in phmmer search')
    sg.add_argument('--database', help='Database to search', default='all', choices=['all', 'pdb', 'afdb', 'bfvd', 'esmfold'])
    sg.add_argument('-v', '--version', action='version', version='%(prog)s version: ' + __version__)


def parse_command_line():
    """Parse MrParse command line arguments"""
    # Read config file, check for local config file for documentation
    if Path(__file__).joinpath("..", "data", "mrparse.config").exists():
        config_file = Path(__file__).joinpath("..", "data", "mrparse.config")
    else:
        config_file = Path(os.environ["CCP4"], "share", "mrparse", "data", "mrparse.config")
    defaults = {}
    config = ConfigParser.SafeConfigParser()
    config.read(str(config_file))
    defaults.update(dict(config.items("Defaults")))
    defaults.update(dict(config.items("Executables")))
    defaults.update(dict(config.items("Databases")))

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    mrparse_argparse(parser)
    parser.set_defaults(**defaults)
    args = parser.parse_args()

    # Add executables and databases to config file so that it only needs to be specified once
    update_config = False
    if args.deeptmhmm_exe != defaults['deeptmhmm_exe']:
        config.set('Executables', 'deeptmhmm_exe', args.deeptmhmm_exe)
        update_config = True
    if args.deepcoil_exe != defaults['deepcoil_exe']:
        config.set('Executables', 'deepcoil_exe', args.deepcoil_exe)
        update_config = True
    if args.hhsearch_exe != defaults['hhsearch_exe']:
        config.set('Executables', 'hhsearch_exe', args.hhsearch_exe)
        update_config = True
    if args.hhsearch_db != defaults['hhsearch_db']:
        config.set('Databases', 'hhsearch_db', args.hhsearch_db)
        update_config = True
    #if args.afdb_seqdb != defaults['afdb_seqdb']:
    #    config.set('Databases', 'afdb_seqdb', args.afdb_seqdb)
    #    update_config = True
    if update_config:
        with open(config_file, 'w') as f:
            config.write(f)

    return args
