from collections import defaultdict


def file_to_dict(filename):
    output_fastq = dict()
    with open(filename) as file:
        while True:
            name = file.readline().strip()
            if not name:
                break
            sequence = file.readline().strip()
            file.readline()
            quality = file.readline().strip()
            reads[name] = (sequence, quality)
    return output_fastq


def filtrated_fastq_to_file(filtrated_sequences):
    pass


def gc_perc(seq):
    """
    Computes the percentage of G and C in the sequence

    Arguments:
    seq: str

    Returns int
    """
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
