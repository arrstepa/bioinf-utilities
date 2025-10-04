from modules.filter_fastq_module import *
from modules.dna_rna_tools_module import *


def run_dna_rna_tools(*args):
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


def filter_fastq(seqs, gc_bounds, length_bounds, quality_threshold):
    filtrated_sequences = []
    for name, (sequence, quality) in seqs.items():
        if (gc_filter(sequence, bounds(gc_bounds)) and
                length_filter(sequence, bounds(length_bounds)) and
                quality_filter(quality, quality_threshold)):
            filtrated_sequences.append(sequence)
    return filtrated_sequences
