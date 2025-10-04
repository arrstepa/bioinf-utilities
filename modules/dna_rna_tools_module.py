def is_dna(sequence):
    """
    Checks if sequence is DNA

    Arguments:
    sequence: str

    Returns bool
    """
    dna_nucleotides = set('ATGCatgc')
    unique_nucleotides = set(sequence)
    return unique_nucleotides <= dna_nucleotides


def is_rna(sequence):
    """
    Checks if sequence is RNA

    Arguments:
    sequence: str

    Returns bool
    """
    rna_nucleotides = set('AUGCaugc')
    unique_nucleotides = set(sequence)
    return unique_nucleotides <= rna_nucleotides


def is_nucleic_acid(sequence):
    """
    Checks if sequence is nucleic acid

    Arguments:
    sequence: str

    Returns bool
    """
    return is_dna(sequence) != is_rna(sequence)


def transcribe(sequence):
    """
    Replaces T/t with U/u, returns mRNA for sense DNA chain

    Arguments:
    sequence: str

    Returns str
    Raises warning if sequence not DNA
    """
    if is_dna(sequence):
        return sequence.replace('T', 'U').replace('t', 'u')
    else:
        return 'not DNA'


def reverse(sequence):
    """
    Returns backwards sequence

    Arguments:
    sequence: str

    Returns str
    """
    return sequence[::-1]


def complement(sequence):
    """
    Returns the chain according to the complementary rule: A=T(U),G=C

    Arguments:
    sequence: str

    Returns str
    """
    if is_dna(sequence):
        complementary_rule = str.maketrans('ATGCatgc', 'TACGtacg')
        return sequence.translate(complementary_rule)
    else:
        complementary_rule = str.maketrans('AUGCaugc', 'UACGuacg')
        return sequence.translate(complementary_rule)


def reverse_complement(sequence):
    """
    Returns backwards and complementary sequence

    Arguments:
    sequence: str

    Returns str
    """
    return reverse(complement(sequence))


def drop_not_na(list_of_seq):
    """
    Removes sequences that are not nucleic acids

    Arguments:
    list_of_seq: list

    Returns list
    """
    return [seq for seq in list_of_seq if is_nucleic_acid(seq)]


def process_sequences(list_of_seq, func):
    """
    Applies required function to the list of sequences

    Arguments:
    list_of_seq: list
    func: function

    Returns list
    """
    return [func(seq) for seq in list_of_seq]