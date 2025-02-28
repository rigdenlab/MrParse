#!/usr/bin/env ccp4-python
# encoding: utf-8
import os
import sys
import traceback
from mrparse import mr_analyse
from mrparse import mr_args


if "CCP4" not in os.environ:
    raise RuntimeError("Cannot find CCP4 installation - please make sure CCP4 is installed and the setup scripts have been run!")


def main():
    args = mr_args.parse_command_line()
    try:
        return mr_analyse.run(args.seqin,
                              hklin=args.hklin,
                              run_serial=args.run_serial,
                              do_classify=args.do_classify,
                              pdb_dir=args.pdb_dir,
                              pdb_local=args.pdb_local,
                              phmmer_dblvl=args.phmmer_dblvl,
                              plddt_cutoff=args.plddt_cutoff,
                              search_engine=args.search_engine,
                              deeptmhmm_exe=args.deeptmhmm_exe,
                              deepcoil_exe=args.deepcoil_exe,
                              hhsearch_exe=args.hhsearch_exe,
                              hhsearch_db=args.hhsearch_db,
                              afdb_seqdb=args.afdb_seqdb,
                              bfvd_seqdb=args.bfvd_seqdb,
                              esm_seqdb=args.esmfold_seqdb,
                              pdb_seqdb=args.pdb_seqdb,
                              ccp4cloud=args.ccp4cloud,
                              use_api=args.use_api,
                              max_hits=args.max_hits,
                              database=args.database,
                              phmmer_exe=args.phmmer_exe,
                              nproc=args.nproc)
    except KeyboardInterrupt:
        sys.stderr.write("Interrupted by keyboard!")
        return 0
    except Exception as e:
        sys.stderr.write("Error running mrparse: {}\n".format(e))
        traceback.print_exc(file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
