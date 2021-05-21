import itertools as itt

import bgreference


CB = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def mut_key_generator():

    subs = ['CA', 'CG', 'CT', 'TA', 'TC', 'TG']
    for s in sorted(subs):
        for c in sorted(itt.product({'A', 'C', 'G', 'T'}, repeat=2)):
            yield c[0] + s[0] + c[1] + '>' + s[1]


KEYS = list(mut_key_generator())


def pyrimidine_center(key):
    """if key of the form ABC>D is not pyrimidine-centered, applies reverse complement"""

    if key[1] not in ['C', 'T']:
        return CB[key[2]] + CB[key[1]] + CB[key[0]] + '>' + CB[key[4]]
    return key


def get_triplet_change(hash):

    chr_, pos_, alt_ = tuple(hash.split('.'))
    pos_ = int(pos_)
    ref_triplet = bgreference.hg38(chr_, start=pos_-1, size=3)
    ntchange = pyrimidine_center(ref_triplet + '>' + alt_)

    return ntchange


def compute_mutrate_vector(gene, ttype, mutrate_table=None):

    t = mutrate_table[(mutrate_table['gene'] == gene) & (mutrate_table['ttype'] == ttype)][KEYS]
    assert(len(t) > 0)
    return t.values


def set_mutrate(gene, ttype, hashes, mutrate_table=None):

    assert(mutrate_table is not None)

    t = mutrate_table[(mutrate_table['gene'] == gene) & (mutrate_table['ttype'] == ttype)][KEYS]
    assert(len(t) > 0)

    mutrate_dict = dict(zip(KEYS, list(t.values[0])))
    triplet_changes_vector = map(get_triplet_change, hashes)
    return list(map(lambda x: mutrate_dict[x], triplet_changes_vector))


def zip_vectors(binary_vector, num_vector):

    ones = [[b for _ in range(int(a))] for a, b in zip(binary_vector, num_vector) if a > 0]
    ones = list(itt.chain(*ones))
    zeros = [b for a, b in zip(binary_vector, num_vector) if a == 0]
    return ones, zeros
