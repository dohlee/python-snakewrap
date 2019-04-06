import itertools
import os
import re

def get_common_prefix(strings, include_ext=False, full_path=False):
    if not full_path:
        strings = [os.path.basename(x) for x in strings]

    if not include_ext:
        strings = [os.path.splitext(x)[0] for x in strings]

    all_same = lambda x: all(x[0] == y for y in x)
    prefix_tuples = itertools.takewhile(all_same, zip(*strings))

    return ''.join(x[0] for x in prefix_tuples)

def get_core_name(s):
    return os.path.splitext(os.path.basename(s))[0]

def extract_with_regex(s, pattern):
    found = re.search(pattern, s)
    if found is not None:
        return found.group(1)
    return None