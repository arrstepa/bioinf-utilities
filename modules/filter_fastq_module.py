import os

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
    os.makedirs('filtered', exist_ok=True)
    file_path = os.path.join(output_fastq, 'filtered', 'output_fastq.txt')
    with open(file_path, 'a') as file:
        for name, sequence in filtrated_sequences.items():
            file.write(f"{name}\n{sequence[0]}\n{name.replace('@', '+')}\n{sequence[1]}")


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
