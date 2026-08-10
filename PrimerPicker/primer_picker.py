#!/usr/bin/env python3
"""Pick multiplex PCR primer pools for DoTA-seq.

This is a command-line version of the original three-notebook workflow:

1. Generate candidate primer pairs for each FASTA target with primer3_core.
2. Score primer-pair dimerization with ntthal.
3. Use simulated annealing to choose low-dimer primer pools.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import os
import pickle
import random
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


DEFAULT_P5 = "CTGCGTGTCTCCGACTCAGACT"
DEFAULT_P7 = "CAAGCAGAAGACGGCATACGAGAT"
DEFAULT_16S_F = "cctacgggaggcagcagt"
DEFAULT_16S_R = "ggactaccagggtatctaatcctgt"


@dataclass
class Target:
    name: str
    sequence: str


@dataclass
class PrimerPair:
    id: str
    gene: str
    fseq: str
    rseq: str
    ftm: str | float
    rtm: str | float
    amplicon_len: str | int


@dataclass
class Gene:
    name: str
    primers: list[PrimerPair]


def parse_fasta(path: Path, exclude_terms: Iterable[str]) -> list[Target]:
    targets: list[Target] = []
    current_name: str | None = None
    sequence_parts: list[str] = []
    exclude_terms = [term for term in exclude_terms if term]

    def flush() -> None:
        nonlocal current_name, sequence_parts
        if current_name is None:
            return
        if not any(term in current_name for term in exclude_terms):
            sequence = "".join(sequence_parts).replace("-", "").replace(" ", "").strip()
            if sequence:
                targets.append(Target(current_name, sequence))
        current_name = None
        sequence_parts = []

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                current_name = line[1:].split()[0]
            else:
                sequence_parts.append(line)
        flush()

    if not targets:
        raise ValueError(f"No usable FASTA records found in {path}")
    return targets


def primer3_settings(num_primers: int, product_size_range: str) -> str:
    return "\n".join(
        [
            "PRIMER_TASK=generic",
            "PRIMER_PICK_LEFT_PRIMER=1",
            "PRIMER_PICK_INTERNAL_OLIGO=0",
            "PRIMER_PICK_RIGHT_PRIMER=1",
            "PRIMER_OPT_TM=60",
            "PRIMER_MIN_TM=55",
            "PRIMER_MAX_TM=63",
            "PRIMER_OPT_SIZE=20",
            "PRIMER_MIN_SIZE=18",
            "PRIMER_MAX_SIZE=25",
            "PRIMER_MAX_END_GC=2",
            "PRIMER_MAX_GC=50",
            "PRIMER_SALT_DIVALENT=2",
            "PRIMER_SALT_MONOVALENT=100",
            f"PRIMER_PRODUCT_SIZE_RANGE={product_size_range}",
            f"PRIMER_NUM_RETURN={num_primers}",
            "PRIMER_EXPLAIN_FLAG=0",
            "P3_FILE_FLAG=0",
            "=",
        ]
    )


def parse_boulder_output(stdout: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for line in stdout.splitlines():
        if not line or line == "=":
            continue
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        parsed[key] = value
    return parsed


def generate_primers(
    targets: list[Target],
    primer3_core: Path,
    num_primers: int,
    product_size_range: str,
) -> list[Gene]:
    settings = primer3_settings(num_primers, product_size_range)
    genes: list[Gene] = []

    for target in targets:
        boulder_input = (
            f"SEQUENCE_ID={target.name}\n"
            f"SEQUENCE_TEMPLATE={target.sequence}\n"
            f"{settings}\n"
        )
        result = subprocess.run(
            [str(primer3_core)],
            input=boulder_input,
            text=True,
            capture_output=True,
            check=False,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"primer3_core failed for {target.name} with exit code "
                f"{result.returncode}:\n{result.stderr.strip()}"
            )

        parsed = parse_boulder_output(result.stdout)
        errors = [value for key, value in parsed.items() if key.endswith("ERROR")]
        if errors:
            raise RuntimeError(f"Primer3 reported an error for {target.name}: {errors}")

        primers: list[PrimerPair] = []
        for index in range(num_primers):
            left_key = f"PRIMER_LEFT_{index}_SEQUENCE"
            right_key = f"PRIMER_RIGHT_{index}_SEQUENCE"
            if left_key not in parsed or right_key not in parsed:
                break
            primers.append(
                PrimerPair(
                    id=f"{target.name}-{index}",
                    gene=target.name,
                    fseq=parsed[left_key],
                    rseq=parsed[right_key],
                    ftm=parsed.get(f"PRIMER_LEFT_{index}_TM", ""),
                    rtm=parsed.get(f"PRIMER_RIGHT_{index}_TM", ""),
                    amplicon_len=parsed.get(f"PRIMER_PAIR_{index}_PRODUCT_SIZE", ""),
                )
            )

        if not primers:
            explain = parsed.get("PRIMER_PAIR_EXPLAIN", "no explanation from Primer3")
            raise RuntimeError(f"No primer pairs generated for {target.name}: {explain}")
        genes.append(Gene(target.name, primers))

    return genes


def add_fixed_primers(genes: list[Gene], include_adapters: bool, include_16s: bool) -> int:
    fixed: list[Gene] = []
    if include_adapters:
        pair = PrimerPair("P5/P7-0", "P5/P7", DEFAULT_P5, DEFAULT_P7, 0, 0, 0)
        fixed.append(Gene("P5/P7", [pair]))
    if include_16s:
        pair = PrimerPair("16S-0", "16S", DEFAULT_16S_F, DEFAULT_16S_R, 0, 0, 500)
        fixed.append(Gene("16S", [pair]))
    genes[:0] = fixed
    return len(fixed)


def write_primer_outputs(genes: list[Gene], outdir: Path) -> None:
    fasta_path = outdir / "primers.fasta"
    tsv_path = outdir / "primers.tsv"
    pickle_path = outdir / "primers.pickle"

    with fasta_path.open("w") as handle:
        for gene in genes:
            for primer in gene.primers:
                handle.write(f">{primer.id}-f\n{primer.fseq}\n")
                handle.write(f">{primer.id}-r\n{primer.rseq}\n")

    with tsv_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["primer_id", "gene", "forward", "reverse", "forward_tm", "reverse_tm", "amplicon_len"])
        for gene in genes:
            for primer in gene.primers:
                writer.writerow(
                    [
                        primer.id,
                        primer.gene,
                        primer.fseq,
                        primer.rseq,
                        primer.ftm,
                        primer.rtm,
                        primer.amplicon_len,
                    ]
                )

    with pickle_path.open("wb") as handle:
        pickle.dump(genes, handle)


def truncate_for_ntthal(sequence: str) -> str:
    return sequence[-60:] if len(sequence) > 60 else sequence


def parse_ntthal_delta_g(stdout: str) -> float:
    parts = stdout.split()
    if len(parts) > 13:
        try:
            return float(parts[13])
        except ValueError:
            pass
    match = re.search(r"(-?\d+(?:\.\d+)?)", stdout)
    if match:
        return float(match.group(1))
    raise ValueError(f"Could not parse ntthal output: {stdout!r}")


def run_ntthal(ntthal: str, s1: str, s2: str, alignment: str) -> float:
    result = subprocess.run(
        [
            ntthal,
            "-mv",
            "200",
            "-dv",
            "3",
            "-n",
            "0.2",
            "-t",
            "25",
            "-d",
            "50",
            "-s1",
            s1,
            "-s2",
            s2,
            "-a",
            alignment,
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"ntthal failed: {result.stderr.strip()}")
    return parse_ntthal_delta_g(result.stdout)


def dimer_score(ntthal: str, s1: str, s2: str, cache: dict[tuple[str, str, str], float]) -> float:
    scores = []
    for alignment in ("END1", "END2"):
        key = (s1, s2, alignment)
        if key not in cache:
            cache[key] = run_ntthal(ntthal, s1, s2, alignment)
        scores.append(cache[key])
    return min(scores)


def score_gene_against_later_genes(
    i: int,
    genes: list[Gene],
    ntthal: str,
) -> dict[frozenset[str], float]:
    cache: dict[tuple[str, str, str], float] = {}
    scores: dict[frozenset[str], float] = {}

    for primer in genes[i].primers:
        seqs_a = [truncate_for_ntthal(primer.fseq), truncate_for_ntthal(primer.rseq)]
        homo_a = [dimer_score(ntthal, seq, seq, cache) for seq in seqs_a]

        for later_gene in genes[i + 1 :]:
            for other in later_gene.primers:
                seqs_b = [truncate_for_ntthal(other.fseq), truncate_for_ntthal(other.rseq)]
                delta_gs = list(homo_a)
                delta_gs.extend(dimer_score(ntthal, seq, seq, cache) for seq in seqs_b)
                for seq1, seq2 in itertools.permutations(seqs_a + seqs_b, 2):
                    delta_gs.append(dimer_score(ntthal, seq1, seq2, cache))
                scores[frozenset((primer.id, other.id))] = min(delta_gs)

    return scores


def calculate_dimer_scores(genes: list[Gene], ntthal: Path, jobs: int) -> dict[frozenset[str], float]:
    ntthal_text = str(ntthal)
    combined: dict[frozenset[str], float] = {}
    workers = max(1, jobs)

    if workers == 1:
        for i in range(len(genes)):
            combined.update(score_gene_against_later_genes(i, genes, ntthal_text))
        return combined

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(score_gene_against_later_genes, i, genes, ntthal_text): i
            for i in range(len(genes))
        }
        for future in as_completed(futures):
            combined.update(future.result())
    return combined


def write_dimer_scores(scores: dict[frozenset[str], float], outdir: Path) -> None:
    with (outdir / "dimerization-deltaGs.pickle").open("wb") as handle:
        pickle.dump(scores, handle)
    with (outdir / "dimerization-deltaGs.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["primer_id_1", "primer_id_2", "delta_g"])
        for key, delta_g in sorted(scores.items(), key=lambda item: tuple(sorted(item[0]))):
            primer_1, primer_2 = sorted(key)
            writer.writerow([primer_1, primer_2, delta_g])


def matrix_scores(primer_set: list[PrimerPair], dimer_scores: dict[frozenset[str], float]) -> list[list[float | None]]:
    size = len(primer_set)
    matrix: list[list[float | None]] = [[None for _ in range(size)] for _ in range(size)]
    for i in range(size):
        for j in range(i):
            score = dimer_scores[frozenset((primer_set[i].id, primer_set[j].id))]
            matrix[i][j] = score
            matrix[j][i] = score
    return matrix


def matrix_values(matrix: list[list[float | None]]) -> list[float]:
    return [value for row in matrix for value in row if value is not None]


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def score_summary(matrix: list[list[float | None]]) -> tuple[float, float, float]:
    values = matrix_values(matrix)
    avg_score = mean(values)
    min_score = min(values)
    rank_score = avg_score + min_score
    return rank_score, avg_score, min_score


def simulated_annealing(
    genes: list[Gene],
    dimer_scores: dict[frozenset[str], float],
    fixed_count: int,
    runtime_minutes: float,
    temperature: float,
    beta: float,
    seed: int | None,
) -> tuple[list[list[float | None]], list[PrimerPair], list[float], list[float], list[float], int]:
    rng = random.Random(seed)
    primer_set = [rng.choice(gene.primers) for gene in genes]
    current_matrix = matrix_scores(primer_set, dimer_scores)
    best_matrix = [row[:] for row in current_matrix]
    best_set = list(primer_set)
    best_rank, _, _ = score_summary(best_matrix)
    min_history: list[float] = []
    avg_history: list[float] = []
    temp_history: list[float] = []
    end_time = time.time() + runtime_minutes * 60
    cycles = 0

    changeable_positions = list(range(fixed_count, len(genes)))
    if not changeable_positions:
        return best_matrix, best_set, min_history, avg_history, temp_history, cycles

    while time.time() < end_time:
        candidate_set = list(primer_set)
        row_sums = [
            sum(value for value in row if value is not None)
            for row in current_matrix
        ]
        if rng.random() >= 0.2:
            change_index = min(changeable_positions, key=lambda idx: row_sums[idx])
        else:
            change_index = rng.choice(changeable_positions)

        candidate_set[change_index] = rng.choice(genes[change_index].primers)
        candidate_matrix = matrix_scores(candidate_set, dimer_scores)

        current_rank, current_avg, current_min = score_summary(current_matrix)
        candidate_rank, candidate_avg, candidate_min = score_summary(candidate_matrix)
        current_energy = -current_avg
        candidate_energy = -candidate_avg

        if current_energy > candidate_energy:
            primer_set = candidate_set
            current_matrix = candidate_matrix
        elif temperature > 0 and math.exp(-(candidate_energy - current_energy) / temperature) > rng.random():
            primer_set = candidate_set
            current_matrix = candidate_matrix

        if candidate_rank > best_rank:
            best_rank = candidate_rank
            best_matrix = candidate_matrix
            best_set = candidate_set

        cycles += 1
        if cycles % 10 == 0:
            min_history.append(current_min)
            avg_history.append(current_avg)
            temp_history.append(temperature)
        if cycles % 100 == 0:
            temperature *= beta

    return best_matrix, best_set, min_history, avg_history, temp_history, cycles


def anneal_worker(
    run_index: int,
    genes: list[Gene],
    dimer_scores: dict[frozenset[str], float],
    fixed_count: int,
    runtime_minutes: float,
    temperature: float,
    beta: float,
    seed: int | None,
) -> dict[str, object]:
    run_seed = None if seed is None else seed + run_index
    matrix, primer_set, min_history, avg_history, temp_history, cycles = simulated_annealing(
        genes,
        dimer_scores,
        fixed_count,
        runtime_minutes,
        temperature,
        beta,
        run_seed,
    )
    _, avg_score, min_score = score_summary(matrix)
    return {
        "matrix": matrix,
        "primer_set": primer_set,
        "avg_score": avg_score,
        "min_score": min_score,
        "rank_score": avg_score + min_score,
        "min_history": min_history,
        "avg_history": avg_history,
        "temp_history": temp_history,
        "cycles": cycles,
    }


def run_annealing(
    genes: list[Gene],
    dimer_scores: dict[frozenset[str], float],
    fixed_count: int,
    iterations: int,
    jobs: int,
    runtime_minutes: float,
    temperature: float,
    beta: float,
    seed: int | None,
) -> list[dict[str, object]]:
    workers = max(1, jobs)
    if workers == 1:
        return [
            anneal_worker(i, genes, dimer_scores, fixed_count, runtime_minutes, temperature, beta, seed)
            for i in range(iterations)
        ]

    results: list[dict[str, object]] = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(
                anneal_worker,
                i,
                genes,
                dimer_scores,
                fixed_count,
                runtime_minutes,
                temperature,
                beta,
                seed,
            )
            for i in range(iterations)
        ]
        for future in as_completed(futures):
            results.append(future.result())
    return results


def write_top_sets(results: list[dict[str, object]], outdir: Path, top_n: int) -> None:
    results.sort(key=lambda result: float(result["rank_score"]), reverse=True)

    with (outdir / "top-primer-sets.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "set_rank",
                "primer_id",
                "gene",
                "forward",
                "reverse",
                "amplicon_len",
                "set_avg_delta_g",
                "set_min_delta_g",
            ]
        )
        for rank, result in enumerate(results[:top_n], start=1):
            for primer in result["primer_set"]:
                assert isinstance(primer, PrimerPair)
                writer.writerow(
                    [
                        rank,
                        primer.id,
                        primer.gene,
                        primer.fseq,
                        primer.rseq,
                        primer.amplicon_len,
                        result["avg_score"],
                        result["min_score"],
                    ]
                )

    with (outdir / "all-primer-sets.pickle").open("wb") as handle:
        pickle.dump(results, handle)


def resolve_executable(value: str) -> str:
    expanded = os.path.expanduser(value)
    if os.sep in expanded:
        path = Path(expanded).resolve()
        if not path.exists():
            raise FileNotFoundError(f"Executable not found: {path}")
        return str(path)

    found = shutil.which(expanded)
    if not found:
        raise FileNotFoundError(
            f"Executable '{value}' was not found in PATH. "
            f"Install Primer3 or pass its location with --primer3-core/--ntthal."
        )
    return found


def readable_path(value: str) -> Path:
    return Path(os.path.expanduser(value)).resolve()


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than 0")
    return parsed


def positive_float(value: str) -> float:
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than 0")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate and pool multiplex PCR primers from a FASTA file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("targets", type=readable_path, help="Input FASTA file of target sequences.")
    parser.add_argument("--outdir", type=readable_path, default=Path("primer_picker_results"), help="Output directory.")
    parser.add_argument("--primer3-core", default="primer3_core", help="primer3_core command or path.")
    parser.add_argument("--ntthal", default="ntthal", help="ntthal command or path.")
    parser.add_argument("--num-primers", type=positive_int, default=100, help="Candidate primer pairs per target.")
    parser.add_argument("--product-size-range", default="400-500", help="Primer3 PRIMER_PRODUCT_SIZE_RANGE.")
    parser.add_argument("--exclude", action="append", default=[], help="Skip FASTA records whose names contain this text.")
    parser.add_argument("--no-adapters", action="store_true", help="Do not include fixed P5/P7 adapter primers.")
    parser.add_argument("--no-16s", action="store_true", help="Do not include fixed 16S primers.")
    parser.add_argument("--dimer-jobs", type=positive_int, default=max(1, os.cpu_count() or 1), help="Parallel ntthal jobs.")
    parser.add_argument("--anneal-jobs", type=positive_int, default=max(1, os.cpu_count() or 1), help="Parallel annealing jobs.")
    parser.add_argument("--anneal-iterations", type=positive_int, default=200, help="Independent annealing runs.")
    parser.add_argument("--anneal-minutes", type=positive_float, default=0.2, help="Minutes per annealing run.")
    parser.add_argument("--temperature", type=positive_float, default=600.0, help="Initial simulated annealing temperature.")
    parser.add_argument("--beta", type=positive_float, default=0.997, help="Cooling multiplier applied every 100 cycles.")
    parser.add_argument("--top-n", type=positive_int, default=10, help="Number of primer pools to write.")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducible annealing.")
    parser.add_argument(
        "--skip-primer-generation",
        action="store_true",
        help="Reuse outdir/primers.pickle instead of running primer3_core.",
    )
    parser.add_argument(
        "--skip-dimer-scoring",
        action="store_true",
        help="Reuse outdir/dimerization-deltaGs.pickle instead of running ntthal.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    if not args.skip_primer_generation and not args.targets.exists():
        raise FileNotFoundError(f"Input FASTA not found: {args.targets}")
    primer3_core = None if args.skip_primer_generation else resolve_executable(args.primer3_core)
    ntthal = None if args.skip_dimer_scoring else resolve_executable(args.ntthal)

    primers_pickle = args.outdir / "primers.pickle"
    if args.skip_primer_generation:
        print(f"Loading primers from {primers_pickle}")
        with primers_pickle.open("rb") as handle:
            genes = pickle.load(handle)
        fixed_count = sum(1 for gene in genes if len(gene.primers) == 1 and gene.name in {"P5/P7", "16S"})
    else:
        print(f"Reading targets from {args.targets}")
        targets = parse_fasta(args.targets, args.exclude)
        print(f"Generating primers for {len(targets)} targets")
        assert primer3_core is not None
        genes = generate_primers(targets, Path(primer3_core), args.num_primers, args.product_size_range)
        fixed_count = add_fixed_primers(genes, not args.no_adapters, not args.no_16s)
        write_primer_outputs(genes, args.outdir)

    total_pairs = sum(len(gene.primers) for gene in genes)
    print(f"Using {total_pairs} primer pairs across {len(genes)} primer groups")
    if len(genes) < 2:
        raise ValueError("At least two primer groups are required for dimer scoring and pooling")

    dimer_pickle = args.outdir / "dimerization-deltaGs.pickle"
    if args.skip_dimer_scoring:
        print(f"Loading dimer scores from {dimer_pickle}")
        with dimer_pickle.open("rb") as handle:
            dimer_scores = pickle.load(handle)
    else:
        print(f"Scoring primer dimers with {args.dimer_jobs} job(s); this is usually the slow step")
        started = time.time()
        assert ntthal is not None
        dimer_scores = calculate_dimer_scores(genes, Path(ntthal), args.dimer_jobs)
        write_dimer_scores(dimer_scores, args.outdir)
        print(f"Dimer scoring finished in {time.time() - started:.1f} seconds")

    print(
        f"Running {args.anneal_iterations} annealing iteration(s) "
        f"for {args.anneal_minutes} minute(s) each"
    )
    results = run_annealing(
        genes,
        dimer_scores,
        fixed_count,
        args.anneal_iterations,
        args.anneal_jobs,
        args.anneal_minutes,
        args.temperature,
        args.beta,
        args.seed,
    )
    write_top_sets(results, args.outdir, args.top_n)

    best = max(results, key=lambda result: float(result["rank_score"]))
    print(f"Best avg deltaG: {float(best['avg_score']):.4f}")
    print(f"Best min deltaG: {float(best['min_score']):.4f}")
    print(f"Wrote outputs to {args.outdir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
