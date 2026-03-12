#!/usr/bin/env python3
"""Fetch an OpenTree induced subtree for the species in a CSV file."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path


PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_SPECIES_CSV = (
    PROJECT
    / "output"
    / "runs"
    / "mixed_cellpose_yolo_full_dataset_v1_bgclean"
    / "linked_species_stats"
    / "species_overview.csv"
)
DEFAULT_OUTPUT_DIR = (
    PROJECT
    / "output"
    / "runs"
    / "mixed_cellpose_yolo_full_dataset_v1_bgclean"
    / "phylogenetic_analysis"
    / "opentree"
)

TNRS_URL = "https://api.opentreeoflife.org/v3/tnrs/match_names"
SUBTREE_URL = "https://api.opentreeoflife.org/v3/tree_of_life/induced_subtree"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-csv", type=Path, default=DEFAULT_SPECIES_CSV)
    parser.add_argument("--species-column", default="species")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def require_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


def read_species(path: Path, column: str) -> list[str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"CSV is empty: {path}")
    if column not in rows[0]:
        raise ValueError(f"Column {column!r} not found in {path}")
    species = []
    seen = set()
    for row in rows:
        value = str(row.get(column, "")).strip()
        if value and value not in seen:
            species.append(value)
            seen.add(value)
    if not species:
        raise ValueError(f"No species found in {path} column {column!r}")
    return species


def expand_species_name(label: str) -> str:
    label = str(label).strip()
    if label.startswith("D. "):
        return "Desmognathus " + label.split(" ", 1)[1]
    return label


def run_curl_json(url: str, payload: dict[str, object]) -> dict[str, object]:
    result = subprocess.run(
        [
            "curl",
            "-sS",
            "-X",
            "POST",
            url,
            "-H",
            "content-type: application/json",
            "--data",
            json.dumps(payload),
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    return json.loads(result.stdout)


def choose_match(name: str, matches: list[dict[str, object]]) -> dict[str, object]:
    if not matches:
        raise ValueError(f"No OpenTree matches for {name}")

    def key(match: dict[str, object]) -> tuple[int, int, float]:
        return (
            0 if not match.get("is_approximate_match") else 1,
            0 if not match.get("is_synonym") else 1,
            -float(match.get("score", 0.0)),
        )

    ranked = sorted(matches, key=key)
    return ranked[0]


def main() -> None:
    args = parse_args()
    require_exists(args.species_csv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    species_labels = read_species(args.species_csv, args.species_column)
    opentree_names = [expand_species_name(name) for name in species_labels]

    tnrs = run_curl_json(TNRS_URL, {"names": opentree_names})
    unmatched = list(tnrs.get("unmatched_names", []))
    if unmatched:
        raise ValueError(f"OpenTree unmatched species: {unmatched}")

    results = tnrs.get("results", [])
    if len(results) != len(species_labels):
        raise ValueError(
            f"Expected {len(species_labels)} TNRS results, got {len(results)}"
        )

    match_rows: list[dict[str, object]] = []
    ott_ids: list[int] = []
    for original_label, expanded_name, result in zip(species_labels, opentree_names, results, strict=True):
        match = choose_match(original_label, list(result.get("matches", [])))
        taxon = dict(match.get("taxon", {}))
        ott_id = int(taxon["ott_id"])
        ott_ids.append(ott_id)
        match_rows.append(
            {
                "input_species": original_label,
                "opentree_query_name": expanded_name,
                "matched_name": str(match.get("matched_name", "")),
                "taxon_name": str(taxon.get("name", "")),
                "unique_name": str(taxon.get("unique_name", "")),
                "ott_id": ott_id,
                "is_synonym": bool(match.get("is_synonym", False)),
                "is_approximate_match": bool(match.get("is_approximate_match", False)),
                "score": float(match.get("score", 0.0)),
                "rank": str(taxon.get("rank", "")),
                "tax_sources": "|".join(str(x) for x in taxon.get("tax_sources", [])),
            }
        )

    subtree = run_curl_json(SUBTREE_URL, {"ott_ids": ott_ids})
    newick = str(subtree.get("newick", "")).strip()
    if not newick:
        raise ValueError("OpenTree induced subtree response was missing newick")

    (args.output_dir / "opentree_tnrs.json").write_text(json.dumps(tnrs, indent=2))
    (args.output_dir / "opentree_subtree.json").write_text(json.dumps(subtree, indent=2))
    (args.output_dir / "opentree_subtree.tre").write_text(newick + ("\n" if not newick.endswith("\n") else ""))
    with (args.output_dir / "opentree_matches.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(match_rows[0].keys()))
        writer.writeheader()
        writer.writerows(match_rows)

    print(f"Matched species: {len(match_rows)}")
    print(f"Wrote matches: {args.output_dir / 'opentree_matches.csv'}")
    print(f"Wrote tree: {args.output_dir / 'opentree_subtree.tre'}")


if __name__ == "__main__":
    main()
