def is_dna(sequence):
    dna_nucleotides = set('ATGCatgc')
    unique_nucleotides = set(sequence)
    return unique_nucleotides <= dna_nucleotides


def is_rna(sequence):
    rna_nucleotides = set('AUGCaugc')
    unique_nucleotides = set(sequence)
    return unique_nucleotides <= rna_nucleotides


def is_nucleic_acid(sequence):
    return is_dna(sequence) != is_rna(sequence)


def transcribe(sequence):
    if is_dna(sequence):
        return sequence.replace('T', 'U').replace('t', 'u')
    else:
        return 'not DNA'


def reverse(sequence):
    return sequence[::-1]


def complement(sequence):
    if is_dna(sequence):
        complementary_rule = str.maketrans('ATGCatgc', 'TACGtacg')
        return sequence.translate(complementary_rule)
    else:
        complementary_rule = str.maketrans('AUGCaugc', 'UACGuacg')
        return sequence.translate(complementary_rule)


def reverse_complement(sequence):
    return reverse(complement(sequence))


def drop_not_na(list_of_seq):
    return [seq for seq in list_of_seq if is_nucleic_acid(seq)]


def process_sequences(list_of_seq, func):
    return [func(seq) for seq in list_of_seq]