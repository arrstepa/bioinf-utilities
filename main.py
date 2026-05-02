import os
from abc import ABC
from collections import Counter

from Bio import SeqIO, SeqUtils
from Bio.SeqRecord import SeqRecord


class BiologicalSequence(ABC):

    """
    Abstract base class for nucleic and aminoacid sequences

    !!! subclasses must have an _alphabet (set of allowed characters) in order to work correctly
    This class supports len(seq), indexing and slicing 
    and it also checks _alphabet
    """

    _alphabet = set()

    def __init__(self, sequence: str):
        super().__init__()
        self.sequence = sequence.upper()
        if not self._check_seq():
            raise ValueError(f"Sequence contains incorrect symbols, check with current alphabet: {self._alphabet}")

    def __len__(self):
        return len(self.sequence)
        
    def __getitem__(self, key):
        return self.__class__(self.sequence[key])
    
    def __str__(self):
        return self.sequence
    
    def _check_seq(self):
        return set(self.sequence) <= self._alphabet



class NucleicAcidSequence(BiologicalSequence):
    """
    Common base class for DNA and RNA sequences
    Subclasses must have _complement - otherwise complementary sequence cannot be found
    This class allows user to make complement sequence, reverse sequence or both
    """
    

    _alphabet = set()
    _complement = {}


    def complement(self):

        if not self._complement:
            raise NotImplementedError('Complement is not defined! Use DNASequence or RNASequence only')
        
        table_of_complement = str.maketrans("".join(self._complement.keys()), "".join(self._complement.values()))

        return self.__class__(self.sequence.translate(table_of_complement))
    
    def reverse(self):
        return self.__class__(self.sequence[::-1])
    
    def reverse_complement(self):
        return self.complement().reverse()



class DNASequence(NucleicAcidSequence):
    """
    Class describing DNA sequences - also provides user with the ability to transcribe DNA to RNA sequence
    """
    
    _alphabet = set('ATGC')
    _complement = {"A": "U", "T": "A", "G": "C", "C": "G"}

    def transcribe(self):
        return RNASequence(self.sequence.replace('T', 'U'))


class RNASequence(NucleicAcidSequence):
    """
    Class describing RNA sequences
    """
    _alphabet = set('AUGC')
    _complement = {"A": "U", "U": "A", "G": "C", "C": "G"}


class AminoAcidSequence(BiologicalSequence):
    """
    Class describing protein sequences
    """
    
    _alphabet = set('ACDEFGHIKLMNPQRSTVWY')

    def aa_composition(self):
        """Counts amino acid composition in the provided sequence"""
        return Counter(self.sequence)







def _formate_bounds(n):

    if isinstance(n, (int, float)):
        return (0, n)
    return tuple(n)

def _is_sequence_good(seqrec, gc_bounds, length_bounds, quality_threshold):
    
    gc_low, gc_hi = _formate_bounds(gc_bounds)
    len_low, len_hi = _formate_bounds(length_bounds)
    
    seq_len = len(seqrec.seq)
    gc = SeqUtils.GC(seqrec.seq)
    quality = seqrec.letter_annotations.get('phred_quality', [])

    return (len_low <= seq_len <= len_hi) and (gc_low <= gc <= gc_hi) and (sum(quality) / len(quality) >= quality_threshold)


def filter_fastq(input_fastq, output_fastq=None, gc_bounds=(0,100), length_bounds=(0, 2**32), quality_threshold=0):
    """
    A set of tools for filtrating fastq files using BioPython

    Arguments:
    input_fastq: str
    output_fastq: str
    gc_bound: int/tuple
    length_bounds: int/tuple
    quality_threshold: int

    Returns path to output_fastq: str
    """

    good_seqrecs = (rec for rec in SeqIO.parse(input_fastq, 'fastq') if _is_sequence_good(rec, gc_bounds, length_bounds, quality_threshold))
    
    SeqIO.write(good_seqrecs, output_fastq, 'fastq')

    return output_fastq

    
def convert_multiline_fasta_to_oneline(input_fasta):
    """
    A tool for converting multiline fasta files to oneline format

    Arguments:
    input_fastq: str

    Returns str
    """
    oneline_fasta = []
    current_seq = ''
    with open(input_fasta) as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    oneline_fasta.append(current_seq)
                    current_seq = ''
                oneline_fasta.append(line)
            else:
                current_seq += line
        if current_seq:
            oneline_fasta.append(current_seq)
    output_fasta = os.path.join(os.path.dirname(input_fasta), 'output_fasta.txt')
    with open(output_fasta, 'a') as file:
        for line in oneline_fasta:
            file.write(line + "\n")

filter_fastq("/Users/arrrstepa/bioinf-utilities/test.fastq", "/Users/arrrstepa/bioinf-utilities/output.fastq", gc_bounds=(0,100), length_bounds=(0,100), quality_threshold=0)