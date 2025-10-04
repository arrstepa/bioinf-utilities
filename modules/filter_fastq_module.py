def gc_perc(seq):
    return ((seq.count('G') + seq.count('C')) / len(seq)) * 100


def gc_filter(seq, bounds=(0,100)):
    return bounds[0] <= gc_perc(seq) <= bounds[1]


def length_filter(seq, bounds=(0, 2**32)):
    return bounds[0] <= len(seq) <= bounds[1]


def quality_decode(quality):
    return [ord(x) - 33 for x in quality]


def quality_filter(quality, threshold=0):
    mean_quality = sum(quality_decode(quality))/len(quality)
    return mean_quality >= threshold


def bounds(x):
    if type(x) == int:
        boundaries = (0, x)
    elif type(x) == tuple:
        boundaries = x
    return boundaries
