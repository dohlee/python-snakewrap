import itertools
import os

def get_common_prefix(strings, include_ext=False, full_path=False):
    if not full_path:
        strings = [os.path.basename(x) for x in strings]

    if not include_ext:
        strings = [os.path.splitext(x)[0] for x in strings]

    all_same = lambda x: all(x[0] == y for y in x)
    prefix_tuples = itertools.takewhile(all_same, zip(*strings))

    return ''.join(x[0] for x in prefix_tuples)