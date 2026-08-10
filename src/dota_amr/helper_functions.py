import gzip
from typing import List

def get_arg_names(primers_file: str) -> List[str]:
    """Load the list of all arg names, from the primers file"""
    arg_names = []
    with open(primers_file, 'r') as f:
        i = 0
        line = f.readline()
        while line != "":
            if i != 0 and i != 1: # skip ove header row and 16s primer
                arg_names.append(line.split(",")[0].strip())
            line = f.readline()
            i += 1
    return arg_names 

def open_maybe_gzip(path: str, mode: str = "rt"):
    """Open plain text or .gz file based on extension."""
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    else:
        return open(path, mode)
    