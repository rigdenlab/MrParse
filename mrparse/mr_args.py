import argparse
import ConfigParser
import os

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


def parse_command_line():

    # Read config file
    config_file = os.path.join(os.environ["CCP4"], "share", "mrparse", "data", "mrparse.config")
    defaults = {}
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)
    defaults.update(dict(config.items("Defaults")))
    defaults.update(dict(config.items("Executables")))

    # Read command line arguments
    parser = argparse.ArgumentParser(description="MrParse Molecular Replacement Search Model analysis")
    parser.set_defaults(**defaults)
    parser.add_argument('-hkl', '--hklin', action=FilePathAction, help='MTZ/CIF Crystal Data file')
    parser.add_argument('--do_classify', action='store_true', help='Run the SS/TM/CC classifiers - requires internet access.')
    parser.add_argument('--pdb_dir', action=FilePathAction, help='Directory of PDB files')
    parser.add_argument('--db_lvl', help='Redundancy level of PDB database', default='95', choices=['95', '100'])
    parser.add_argument('--run_serial', action='store_true', help='Run on a single processor')
    parser.add_argument('-seq', '--seqin', action=FilePathAction, required=True, help='Sequence file')
    parser.add_argument('--tmhmm_exe', action=FilePathAction,
                        help="Location of TMHMM executable for transmembrane classification")
    parser.add_argument('--deepcoil_exe', action=FilePathAction,
                        help="Location of Deepcoil executable for coiled-coil classification")
    parser.add_argument('--ccp4cloud', action='store_true', help="specify running through CCP4Cloud")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version: ' + __version__)
    args = parser.parse_args()

    # Add tmhmm & deepcoil exe to config file so that it only needs to be specified once
    update_config = False
    if args.tmhmm_exe != defaults['tmhmm_exe']:
        config.set('Executables', 'tmhmm_exe', args.tmhmm_exe)
        update_config = True
    if args.deepcoil_exe != defaults['deepcoil_exe']:
        config.set('Executables', 'deepcoil_exe', args.deepcoil_exe)
        update_config = True
    if update_config:
        with open(config_file, 'w') as f:
            config.write(f)

    return args
