#!/usr/bin/env python
from __future__ import annotations

import argparse
import math
import pathlib
import re
from collections import Counter

import pyarrow.parquet as pq


DEFAULT_DATASET_DIR = pathlib.Path("/Users/searle.brian/Documents/huggingface/data/prospect-ptms-charge/data")
ALLOWED_SEQUENCE_RE = re.compile(r"^\[\]-([A-Z](?:\[UNIMOD:4\])?)+-\[\]$")
CAM_TOKEN = "[UNIMOD:4]"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Optimize the Scout seed-charge heuristic "
            "z=min(B+1,max(1,round(1+a*B+b*L))) "
            "against PROSPECT PTMs charge_state_dist labels."
        )
    )
    parser.add_argument(
        "--dataset-dir",
        type=pathlib.Path,
        default=DEFAULT_DATASET_DIR,
        help="Directory containing the train/val/test parquet shards.",
    )
    parser.add_argument(
        "--alpha-start",
        type=float,
        default=0.30,
        help="Starting B coefficient.",
    )
    parser.add_argument(
        "--alpha-stop",
        type=float,
        default=0.90,
        help="Inclusive ending B coefficient.",
    )
    parser.add_argument(
        "--alpha-step",
        type=float,
        default=0.01,
        help="B coefficient step size.",
    )
    parser.add_argument(
        "--beta-start",
        type=float,
        default=0.00,
        help="Starting L coefficient.",
    )
    parser.add_argument(
        "--beta-stop",
        type=float,
        default=0.08,
        help="Inclusive ending L coefficient.",
    )
    parser.add_argument(
        "--beta-step",
        type=float,
        default=0.001,
        help="L coefficient step size.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.01,
        help="Charge-state support threshold. Charges strictly above this are treated as present.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=8192,
        help="Parquet row batch size.",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=10,
        help="Number of best parameter pairs to report.",
    )
    return parser.parse_args()


def frange(start: float, stop: float, step: float) -> list[float]:
    count = int(round((stop - start) / step))
    values = []
    for i in range(count + 1):
        value = start + (i * step)
        values.append(round(value, 10))
    return values


def is_unmodified_or_cam_only(sequence: str) -> bool:
    if not ALLOWED_SEQUENCE_RE.match(sequence):
        return False
    cleaned = sequence.replace(CAM_TOKEN, "")
    return "[UNIMOD:" not in cleaned


def extract_residues(sequence: str) -> str:
    middle = sequence[3:-3]
    return re.sub(r"\[UNIMOD:4\]", "", middle)


def charge_mask(charge_state_dist: list[float], threshold: float) -> int:
    mask = 0
    for index, probability in enumerate(charge_state_dist):
        if probability > threshold:
            mask |= 1 << index
    return mask


def estimate_charge(basic_count: int, length: int, alpha: float, beta: float, max_charge: int) -> int:
    heuristic_charge = round(1.0 + (alpha * basic_count) + (beta * length))
    bounded = min(basic_count + 1, max(1, heuristic_charge))
    return min(max_charge, bounded)


def load_group_counts(
    dataset_dir: pathlib.Path,
    threshold: float,
    batch_size: int,
) -> tuple[Counter[tuple[int, int, int]], int, int]:
    grouped_counts: Counter[tuple[int, int, int]] = Counter()
    total_rows = 0
    eligible_rows = 0
    for parquet_path in sorted(dataset_dir.glob("*.parquet")):
        parquet_file = pq.ParquetFile(parquet_path)
        for batch in parquet_file.iter_batches(
                batch_size=batch_size,
                columns=["modified_sequence", "charge_state_dist"]):
            sequences = batch.column("modified_sequence").to_pylist()
            distributions = batch.column("charge_state_dist").to_pylist()
            total_rows += len(sequences)
            for sequence, distribution in zip(sequences, distributions):
                if not is_unmodified_or_cam_only(sequence):
                    continue
                residues = extract_residues(sequence)
                basic_count = residues.count("K") + residues.count("R")
                mask = charge_mask(distribution, threshold)
                if mask == 0:
                    continue
                grouped_counts[(basic_count, len(residues), mask)] += 1
                eligible_rows += 1
    return grouped_counts, total_rows, eligible_rows


def evaluate_parameter_grid(
    grouped_counts: Counter[tuple[int, int, int]],
    alpha_values: list[float],
    beta_values: list[float],
) -> list[tuple[int, float, float]]:
    results: list[tuple[int, float, float]] = []
    max_charge = 6
    for alpha in alpha_values:
        for beta in beta_values:
            misses = 0
            for (basic_count, length, mask), count in grouped_counts.items():
                seed_charge = estimate_charge(basic_count, length, alpha, beta, max_charge)
                if mask & (1 << (seed_charge - 1)):
                    continue
                misses += count
            results.append((misses, alpha, beta))
    results.sort(key=lambda item: (item[0], abs(item[1] - 0.6), abs(item[2] - 0.03), item[1], item[2]))
    return results


def main() -> None:
    args = parse_args()
    alpha_values = frange(args.alpha_start, args.alpha_stop, args.alpha_step)
    beta_values = frange(args.beta_start, args.beta_stop, args.beta_step)
    grouped_counts, total_rows, eligible_rows = load_group_counts(args.dataset_dir, args.threshold, args.batch_size)
    if eligible_rows == 0:
        raise SystemExit("No eligible unmodified/CAM-only peptides found.")

    results = evaluate_parameter_grid(grouped_counts, alpha_values, beta_values)
    best_misses, best_alpha, best_beta = results[0]
    current_misses = next(
        misses for misses, alpha, beta in results if math.isclose(alpha, 0.6) and math.isclose(beta, 0.03)
    )

    print(f"dataset_dir={args.dataset_dir}")
    print(f"total_rows={total_rows}")
    print(f"eligible_rows={eligible_rows}")
    print(f"threshold={args.threshold}")
    print(f"current_alpha=0.60 current_beta=0.03 current_misses={current_misses} current_hit_rate={(eligible_rows - current_misses) / eligible_rows:.6f}")
    print(f"best_alpha={best_alpha:.3f} best_beta={best_beta:.3f} best_misses={best_misses} best_hit_rate={(eligible_rows - best_misses) / eligible_rows:.6f}")
    print("top_candidates:")
    for misses, alpha, beta in results[:args.top_k]:
        hit_rate = (eligible_rows - misses) / eligible_rows
        print(f"  alpha={alpha:.3f} beta={beta:.3f} misses={misses} hit_rate={hit_rate:.6f}")


if __name__ == "__main__":
    main()
