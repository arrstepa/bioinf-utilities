import os
import argparse
import logging

logging.basicConfig(
    filename='filterfastq.log',
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(message)s'
)


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


def file_to_dict(filename):
    """
    Reads file and return dictionary with names as keys
    and tuple of sequence with quality string as value

    Arguments:
    filename: str

    Returns dict
    """
    reads = dict()
    with open(filename) as file:
        while True:
            name = file.readline().strip()
            if not name:
                break
            sequence = file.readline().strip()
            file.readline()
            quality = file.readline().strip()
            reads[name] = (sequence, quality)
    return reads


def filtrated_fastq_to_file(filtrated_sequences, output_fastq):
    """
    Creates an output file in directory /filtered with filtrated sequences

    Arguments:
    filtrated_sequences: dict
    output_fastq: str
    Returns str
    """


    file_path = os.path.join(output_fastq, 'filtered', 'output_fastq.txt')
    os.makedirs(os.path.join(output_fastq, 'filtered'), exist_ok=True)
   
    with open(file_path, 'w') as file:
        for name, sequence in filtrated_sequences.items():
            file.write(f"{name}\n{sequence[0]}\n{name.replace('@', '+')}\n{sequence[1]}")

    return file_path


def gc_perc(seq):
    """
    Computes the percentage of G and C in the sequence

    Arguments:
    seq: str

    Returns int
    """
    if len(seq) == 0:
        return 0
    else:
        return ((seq.count('G') + seq.count('C')) / len(seq)) * 100


def gc_filter(seq, bounds=(0,100)):
    """
    Removes sequences containing G and C outside the specified range

    Arguments:
    seq: str
    bounds: tuple

    Returns bool
    """
    return bounds[0] <= gc_perc(seq) <= bounds[1]


def length_filter(seq, bounds=(0, 2**32)):
    """
    Removes sequences of length out of the specified range

    Arguments:
    seq: str
    bounds: tuple

    Returns bool
    """
    return bounds[0] <= len(seq) <= bounds[1]


def quality_decode(quality=0):
    """
    Decodes quality string according to ASCII

    Arguments:
    quality: int

    Returns list
    """
    return [ord(x) - 33 for x in quality]


def quality_filter(quality, threshold=0):
    """
    Removes sequences of quality out of the specified range

    Arguments:
    quality: str
    threshold: int

    Returns bool
    """
    mean_quality = sum(quality_decode(quality))/len(quality)
    return mean_quality >= threshold


def bounds(x):
    """
    If a single boundary number is given, it turns it into an interval

    Arguments:
    x: int/tuple

    Returns tuple
    """
    if type(x) == int:
        boundaries = (0, x)
    elif type(x) == tuple:
        boundaries = x
    return boundaries



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
    gc_bounds: int/tuple
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


def main():

    logging.info('Filtering in process...')
    try:
        parser = argparse.ArgumentParser(description='A set of tools for filtrating fastq files')
        
        parser.add_argument('input_fastq', help='Path to FASTQ file')
        parser.add_argument('--gcbounds', required=True, help='GC content bounds, e.g. 20-30')
        parser.add_argument('--length', required=True, help='Length bounds, e.g. 0-2000')
        parser.add_argument('--quality', type=int, required=True, help='Quality threshold, e.g. 50')

    
        args = parser.parse_args()

        gc_bounds = tuple(map(int, args.gcbounds.split('-')))
        length_bounds = tuple(map(int, args.length.split('-')))

        
        result = filter_fastq(
            args.input_fastq,
            gc_bounds,
            length_bounds,
            args.quality
        )
        logging.info(f'Filtered FASTQ was saved at: {result}')
    
    except Exception as e:
        logging.error(f'Filtration failed: {e}', exc_info=True)
        raise

    finally:
        logging.info('Filtration was finished')


if __name__ == '__main__':
    main()
