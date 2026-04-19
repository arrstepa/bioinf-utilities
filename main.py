from modules.filter_fastq_module import *
from modules.dna_rna_tools_module import *
import os


def run_dna_rna_tools(*args):
    """
    A set of tools for working with DNA/RNA sequences

    Arguments:
    *args: str

    Returns str/bool/list
    """
    procedure = args[-1]
    sequences = list(args[:-1])
    sequences_only_na = drop_not_na(list(args[:-1]))

    if procedure == 'is_nucleic_acid':
        results = process_sequences(sequences, is_nucleic_acid)
    elif procedure == 'transcribe':
        results = process_sequences(sequences_only_na, transcribe)
    elif procedure == 'reverse':
        results = process_sequences(sequences_only_na, reverse)
    elif procedure == 'complement':
        results = process_sequences(sequences_only_na, complement)
    elif procedure == 'reverse_complement':
        results = process_sequences(sequences_only_na, reverse_complement)
    else:
        print('Unknown procedure')
        return

    if len(results) == 1:
        print(results[0])
    else:
        print(results)


def filter_fastq(input_fastq, gc_bounds, length_bounds, quality_threshold):
    """
    A set of tools for filtrating fastq files

    Arguments:
    input_fastq: str
    gc_bound: int/tuple
    length_bounds: int/tuple
    quality_threshold: int

    Returns dict
    """
    output_fastq = os.path.dirname(input_fastq)
    seqs = file_to_dict(input_fastq)
    filtrated_sequences = dict()
    for name, (sequence, quality) in seqs.items():
        if (gc_filter(sequence, bounds(gc_bounds)) and
                length_filter(sequence, bounds(length_bounds)) and
                quality_filter(quality, quality_threshold)):
            filtrated_sequences[name]=(sequence,quality)
    return filtrated_fastq_to_file(filtrated_sequences, output_fastq)


filter_fastq("/Users/arrrstepa/bioinf2025_python/bioinf-utilities/example_fastq.fastq", (0,100), 1, 0)