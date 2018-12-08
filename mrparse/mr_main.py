#!/usr/bin/env ccp4-python
# encoding: utf-8
'''
mrparse.main -- shortdesc

mrparse.main is a description

It defines classes_and_methods

@author:     user_name

@copyright:  2018 organization_name. All rights reserved.

@license:    license

@contact:    user_email
@deffield    updated: Updated
'''

import os
import sys
import traceback

sys.path.insert(0,'/opt/MrParse')


from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

import mr_analyse
DEBUG=True


def main():
    '''Command line options.'''
    program_name = os.path.basename(sys.argv[0])
    try:
        # Setup argument parser
        parser = ArgumentParser()
        parser.add_argument('-hklin', help='MTZ/CIF Crystal Data file')
        parser.add_argument('-seqin', help='Sequence file')

        # Process arguments
        args = parser.parse_args()
        mr_analyse.run(args.seqin, hklin=args.hklin)

        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception as e:
        if DEBUG:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_traceback,
                                      file=sys.stdout)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2

if __name__ == "__main__":
    sys.exit(main())