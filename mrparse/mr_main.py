#!/usr/bin/env ccp4-python
# encoding: utf-8
from argparse import ArgumentParser
import logging.config
import json
import os
import sys

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
logging_json = os.path.join(THIS_DIR, 'logging.json')
with open(logging_json, 'rt') as f:
    config = json.load(f)
logging.config.dictConfig(config)
logger = logging.getLogger()

MRPARSE_DIR = os.path.join(THIS_DIR, '..')
sys.path.insert(0, MRPARSE_DIR)
from mrparse import mr_analyse

def main():
    '''Command line options.'''
    program_name = os.path.basename(sys.argv[0])
    logger.info("Starting: %s", program_name)
    try:
        parser = ArgumentParser()
        parser.add_argument('--hklin', help='MTZ/CIF Crystal Data file')
        parser.add_argument('--seqin', required=True, help='Sequence file')
        parser.add_argument('--multip', type=bool, help='Run on multiple processors')
        args = parser.parse_args()
        return mr_analyse.run(args.seqin, hklin=args.hklin, multiprocessing=args.multip)
    except KeyboardInterrupt:
        logger.critcal("Interrupted by keyboard!")
        return 0
    except Exception as e:
        indent = len(program_name) * " "
        logger.error(program_name + ": " + repr(e) + "\n", exc_info=True)
        logger.error(indent + "  for help use --help")
        return 2

if __name__ == "__main__":
    sys.exit(main())
