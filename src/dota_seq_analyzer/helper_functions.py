import csv
import gzip
import os
from typing import Dict, List


# 2026-08-10: Create parent directories for every command-line output path.
# Reason: DoTA-Seq Analyzer separates temporary, report, and figure files without requiring manual mkdir commands.
def ensure_output_directories(*output_paths: str) -> None:
    for output_path in output_paths:
        parent = os.path.dirname(os.path.abspath(output_path))
        os.makedirs(parent, exist_ok=True)

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


# 2026-08-10: Read optional target-analysis modes from the primer CSV.
# Reason: blank targets use standard detection while only SSR and inversion targets need PV analysis.
def get_target_modes(primers_file: str) -> Dict[str, str]:
    """Return non-taxonomy target names mapped to blank, ssr, or inv modes."""
    with open(primers_file, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return {}

        target_column = "Target" if "Target" in reader.fieldnames else "Primer"
        mode_column = "Mode" if "Mode" in reader.fieldnames else None
        legacy_column = "Sub-ARG" if "Sub-ARG" in reader.fieldnames else None
        modes = {}
        for row_number, row in enumerate(reader):
            target = (row.get(target_column) or "").strip()
            if row_number == 0 and target.casefold() == "16s":
                continue
            if mode_column:
                mode = (row.get(mode_column) or "").strip().casefold()
            elif legacy_column:
                # 2026-08-10: Treat legacy single/family files as having no PV mode.
                # Reason: old primer panels remain usable without assigning SSR or inversion analysis.
                mode = ""
            else:
                mode = ""
            modes[target] = mode
        return modes

def open_maybe_gzip(path: str, mode: str = "rt"):
    """Open plain text or .gz file based on extension."""
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    else:
        return open(path, mode)
