#!/bin/bash

afdb_split_dir=<path to AFDB sequence split files>

seqin="${@:(-2):1}"
phmmer_args="${@:1:$#-1}"

ccp4-python -m mrbump.seq_align.phmmer_par --seqin ${seqin} -n 16 -f ${afdb_split_dir} -p

$CCP4/libexec/phmmer ${phmmer_args} afdb_sequences.fasta

