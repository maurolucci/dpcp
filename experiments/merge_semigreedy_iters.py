#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge semigreedy .iter files into wide CSV matrices by variant and metric."
        )
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path("out") / "semigreedy2s",
        help="Directory containing v2 and v3 folders (default: out/semigreedy2s).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("out"),
        help="Directory where consolidated files will be written (default: out).",
    )
    return parser.parse_args()


def parse_iter_file(file_path: Path) -> tuple[dict[int, float], dict[int, float]]:
    value_by_iter: dict[int, float] = {}
    time_by_iter: dict[int, float] = {}

    with file_path.open("r", newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or len(row) < 3:
                continue

            try:
                iter_idx = int(float(row[0]))
                value = float(row[1])
                elapsed = float(row[2])
            except ValueError:
                continue

            value_by_iter[iter_idx] = value
            time_by_iter[iter_idx] = elapsed

    return value_by_iter, time_by_iter


def discover_iter_files(variant_dir: Path) -> list[Path]:
    iter_subdir = variant_dir / "iter"
    if iter_subdir.is_dir():
        return sorted(iter_subdir.glob("*.iter"))
    return sorted(variant_dir.glob("*.iter"))


def write_matrix(output_path: Path, matrix: dict[str, dict[int, float]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    instances = sorted(matrix)

    all_iters: set[int] = set()
    for values in matrix.values():
        all_iters.update(values.keys())
    sorted_iters = sorted(all_iters)

    with output_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["iter", *instances])

        for iter_idx in sorted_iters:
            row = [iter_idx]
            for instance in instances:
                value = matrix[instance].get(iter_idx)
                row.append("" if value is None else value)
            writer.writerow(row)


def consolidate_variant(base_dir: Path, output_dir: Path, variant: str) -> None:
    variant_dir = base_dir / variant
    if not variant_dir.is_dir():
        raise FileNotFoundError(f"Variant directory not found: {variant_dir}")

    iter_files = discover_iter_files(variant_dir)
    if not iter_files:
        raise FileNotFoundError(
            f"No .iter files found for {variant}. Looked in {variant_dir / 'iter'} and {variant_dir}."
        )

    value_matrix: dict[str, dict[int, float]] = {}
    time_matrix: dict[str, dict[int, float]] = {}

    for file_path in iter_files:
        instance_name = file_path.stem
        value_by_iter, time_by_iter = parse_iter_file(file_path)
        value_matrix[instance_name] = value_by_iter
        time_matrix[instance_name] = time_by_iter

    write_matrix(output_dir / f"semigreedy2s-{variant}-value.iter", value_matrix)
    write_matrix(output_dir / f"semigreedy2s-{variant}-time.iter", time_matrix)


def main() -> int:
    args = parse_args()

    consolidate_variant(args.base_dir, args.output_dir, "v2")
    consolidate_variant(args.base_dir, args.output_dir, "v3")

    print("Consolidated files written:")
    print(args.output_dir / "semigreedy2s-v2-value.iter")
    print(args.output_dir / "semigreedy2s-v2-time.iter")
    print(args.output_dir / "semigreedy2s-v3-value.iter")
    print(args.output_dir / "semigreedy2s-v3-time.iter")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
