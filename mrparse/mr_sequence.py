"""
Created on 17 Nov 2018

@author: jmht & hlasimpk
"""

import copy
import os

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import molecular_weight
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

SUFFIX_TO_TYPE = {'fasta': 'fasta'}


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
            self._bio_seq = Seq(sequence, IUPAC.protein)
            self._bio_seq_record = SeqRecord(self._bio_seq)
        self.nresidues = len(self._bio_seq)
        self.sequence = str(self._bio_seq)
        self.sequence = str(self._bio_seq)

    def __len__(self):
        if isinstance(self.nresidues, int):
            return self.nresidues
        return 0

    def _read_sequence_file(self, seq_file, sequence_type):
        if sequence_type is None:
            sequence_type = self.sequence_type_from_filename(seq_file)
            if not sequence_type:
                raise RuntimeError("Cannot determine sequence type from file: {}".format(seq_file))

        try:
            self._bio_seq_record = SeqIO.read(seq_file, sequence_type, alphabet=IUPAC.protein)
            self._bio_seq = self._bio_seq_record.seq
        except ValueError:
            raise MultipleSequenceException

    @staticmethod
    def sequence_type_from_filename(seq_file):
        _, suffix = os.path.splitext(seq_file)
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
            sequence_type = self._sequence_type_from_filename(seq_file)
            if not sequence_type:
                raise RuntimeError("Cannot determine sequence type from file: {}".format(seq_file))
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
        raise RuntimeError("Cannot determine sequence type from file: {}".format(seq_file))

    sequence = ""
    identifier = []
    previous_seqs = []
    for seq in SeqIO.parse(seq_file, sequence_type, alphabet=IUPAC.protein):
        if seq.seq in previous_seqs:
            continue
        sequence += seq.seq
        identifier.append(seq.id)
        previous_seqs.append(seq.seq)

    bio_seq_record = SeqRecord(sequence)
    bio_seq_record.id = "||".join(identifier)

    file_name = os.path.basename(seq_file).split('.')[0]
    merged_seq_file = os.path.join(os.getcwd(), '{}_merged.fasta'.format(file_name))
    SeqIO.write(bio_seq_record, merged_seq_file, "fasta")

    return Sequence(merged_seq_file)
