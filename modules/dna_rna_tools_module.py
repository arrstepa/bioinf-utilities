def is_nucleic_acid(list_of_seq):
    dna_nucleotides = set('ATGCatgc')
    rna_nucleotides = set('AUGCaugc')
    list_of_results = []
    for sequence in list_of_seq:
        unique_nucleotides = set(sequence)
        is_dna = (unique_nucleotides <= dna_nucleotides)
        is_rna = (unique_nucleotides <= rna_nucleotides)
        list_of_results.append(is_dna != is_rna)
    return list_of_results


def transcribe(list_of_seq):
    list_of_results = []
    dna_nucleotides = set('ATGCatgc')
    for sequence in list_of_seq:
        unique_nucleotides = set(sequence)
        is_dna = (unique_nucleotides <= dna_nucleotides)
        if is_dna:
            transcribed = sequence.replace('T', 'U').replace('t', 'u')
            list_of_results.append(transcribed)
        else:
            list_of_results.append('not DNA')
    return list_of_results


def reverse(list_of_seq):
    return [seq[::-1] for seq in list_of_seq]


def complement(list_of_seq):
    list_of_results = []
    dna_nucleotides = set('ATGCatgc')
    for sequence in list_of_seq:
        unique_nucleotides = set(sequence)
        is_dna = (unique_nucleotides <= dna_nucleotides)
        if is_dna:
            complementary_rule_dna = str.maketrans('ATGCatgc', 'TACGtacg')
            list_of_results.append(sequence.translate(complementary_rule_dna))
        else:
            complementary_rule_rna = str.maketrans('AUGCaugc', 'UACGuacg')
            list_of_results.append(sequence.translate(complementary_rule_rna))
    return list_of_results


def reverse_complement(list_of_complements):
    return [seq[::-1] for seq in list_of_complements]


def drop_not_na(list_of_seq):
    existing_seq_only = [
        seq
        for seq, is_nucleic in zip(list_of_seq, is_nucleic_acid(list_of_seq))
        if is_nucleic
    ]
    return existing_seq_only
