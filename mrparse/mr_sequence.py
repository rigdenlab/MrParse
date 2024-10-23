"""
Created on 17 Nov 2018

@author: jmht & hlasimpk
"""

import copy
from pathlib import Path

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import molecular_weight
from Bio.SeqRecord import SeqRecord

SUFFIX_TO_TYPE = {'fasta': 'fasta',
                  'pir': 'pir',
                  'seq': 'fasta',
                  'fastq': 'fastq'}


class MultipleSequenceException(Exception):
    pass


class Sequence(object):
    """Class for handling sequence data"""
    def __init__(self, seq_file=None, sequence=None, sequence_type=None):
        self.sequence_file = seq_file
        self._molecular_weight = None
        self._bio_seq = None
        self._bio_seq_record = None
        if seq_file:
            self._read_sequence_file(seq_file, sequence_type)
        elif sequence:
            self._bio_seq = Seq(sequence)
            self._bio_seq_record = SeqRecord(self._bio_seq, id="")
        self.nresidues = len(self._bio_seq)
        self.remove_dashes()
        self.sequence = str(self._bio_seq)

    def __len__(self):
        if isinstance(self.nresidues, int):
            return self.nresidues
        return 0

    def _read_sequence_file(self, seq_file, sequence_type):
        if sequence_type is None:
            sequence_type = self.sequence_type_from_filename(seq_file)
            if not sequence_type:
                raise RuntimeError(f"Cannot determine sequence type from file: {seq_file}")

        try:
            self._bio_seq_record = SeqIO.read(seq_file, sequence_type)
            self._bio_seq = self._bio_seq_record.seq

        except ValueError:
            raise MultipleSequenceException

    def remove_dashes(self):
        """Remove dashes from sequence and replace with Xs"""
        sequence_str = str(self._bio_seq)
        sequence_str = sequence_str.replace('-', 'X')
        self._bio_seq = Seq(sequence_str)
        self._bio_seq_record.seq = self._bio_seq

    @staticmethod
    def sequence_type_from_filename(seq_file):
        suffix = Path(seq_file).suffix
        suffix = suffix.lstrip('.').lower()
        try:
            return SUFFIX_TO_TYPE[suffix]
        except KeyError:
            return None

    @property
    def molecular_weight(self):
        if self._molecular_weight is None:
            self._calculate_molecular_weight()
        return self._molecular_weight

    def _calculate_molecular_weight(self):
        self._molecular_weight = molecular_weight(self._bio_seq, 'protein')

    def write(self, seq_file, sequence_type=None, description=None):
        """Write sequence out to file seq_file of type sequence_type.

        Parameters
        ----------
        seq_file : str
           The filename of the file
        sequence_type : str
           The type of the sequence (recogniseable to Biopython e.g. 'fasta') 
        description : str
           The text to put on the first line of the file if different from that already in the record
        """
        if sequence_type is None:
            sequence_type = self.sequence_type_from_filename(seq_file)
            if not sequence_type:
                raise RuntimeError(f"Cannot determine sequence type from file: {seq_file}")
        if description:
            seq_record = copy.copy(self._bio_seq_record)
            seq_record.id = description
            seq_record.description = description
        else:
            seq_record = self._bio_seq_record
        SeqIO.write(seq_record, seq_file, sequence_type)


def merge_multiple_sequences(seq_file):
    """
    Function to merge multiple sequences from a fasta file

    Parameters
    ----------
    seq_file : str
        The filename of the file

    Returns
    -------
    seq : :obj:
        mrparse.mr_sequence.Sequence object

    Raises
    ------
    RuntimeError
    """

    sequence_type = Sequence.sequence_type_from_filename(seq_file)
    if not sequence_type:
        raise RuntimeError(f"Cannot determine sequence type from file: {seq_file}")

    sequence = ""
    identifier = []
    previous_seqs = []
    for seq in SeqIO.parse(seq_file, sequence_type):
        if seq.seq in previous_seqs:
            continue
        sequence += seq.seq
        identifier.append(seq.id)
        previous_seqs.append(seq.seq)

    bio_seq_record = SeqRecord(sequence)
    bio_seq_record.id = "||".join(identifier)

    file_name = Path(seq_file).stem
    merged_seq_file = Path.cwd().joinpath(f'{file_name}_merged.fasta')
    SeqIO.write(bio_seq_record, str(merged_seq_file), "fasta")

    return Sequence(merged_seq_file)
