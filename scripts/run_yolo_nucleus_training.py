#!/usr/bin/env python3
"""Train a YOLO segmentation model on the manual nucleus dataset."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
from typing import Callable, Iterable, Iterator, TypeVar

PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_DATASET = PROJECT / "output" / "yolo_nucleus_dataset_round1"
DEFAULT_MODEL = PROJECT / "output" / "yolo_models" / "yolo26n-seg.pt"
DEFAULT_OUTPUT = PROJECT / "output" / "yolo_nucleus_training_round1"
_T = TypeVar("_T")
_U = TypeVar("_U")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-dir", type=Path, default=DEFAULT_DATASET)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--imgsz", type=int, default=1024)
    parser.add_argument("--batch", type=int, default=4)
    parser.add_argument("--workers", type=int, default=0)
    parser.add_argument("--device", default=None)
    parser.add_argument("--patience", type=int, default=15)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def choose_device(requested: str | None) -> str:
    import torch

    if requested:
        return requested
    return "0" if torch.cuda.is_available() else "cpu"


def require_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


def best_metrics(save_dir: Path) -> dict[str, float | int]:
    results_csv = save_dir / "results.csv"
    if not results_csv.exists():
        return {}

    rows = list(csv.DictReader(results_csv.open()))
    if not rows:
        return {}

    parsed_rows: list[dict[str, float]] = []
    for row in rows:
        parsed: dict[str, float] = {}
        for key, value in row.items():
            try:
                parsed[key] = float(value)
            except (TypeError, ValueError):
                continue
        parsed_rows.append(parsed)

    best = max(parsed_rows, key=lambda row: row.get("metrics/mAP50-95(M)", float("-inf")))
    return {
        "epochs_completed": int(len(parsed_rows)),
        "best_epoch_by_mask_map50_95": int(best.get("epoch", -1)),
        "best_mask_precision": float(best.get("metrics/precision(M)", 0.0)),
        "best_mask_recall": float(best.get("metrics/recall(M)", 0.0)),
        "best_mask_map50": float(best.get("metrics/mAP50(M)", 0.0)),
        "best_mask_map50_95": float(best.get("metrics/mAP50-95(M)", 0.0)),
        "best_box_map50": float(best.get("metrics/mAP50(B)", 0.0)),
        "best_box_map50_95": float(best.get("metrics/mAP50-95(B)", 0.0)),
    }


class SerialPool:
    """Minimal ThreadPool-compatible shim that avoids semaphore-backed workers."""

    def __init__(self, *_args, **_kwargs) -> None:
        pass

    def __enter__(self) -> "SerialPool":
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        return False

    def imap(self, func: Callable[[_T], _U], iterable: Iterable[_T]) -> Iterator[_U]:
        for item in iterable:
            yield func(item)

    def map(self, func: Callable[[_T], _U], iterable: Iterable[_T]) -> list[_U]:
        return [func(item) for item in iterable]

    def close(self) -> None:
        return None

    def join(self) -> None:
        return None

    def terminate(self) -> None:
        return None


def patch_ultralytics_sandbox_pools() -> None:
    import ultralytics.data.base as yolo_base
    import ultralytics.data.dataset as yolo_dataset
    import ultralytics.data.utils as yolo_utils

    yolo_base.ThreadPool = SerialPool
    yolo_dataset.ThreadPool = SerialPool
    yolo_utils.ThreadPool = SerialPool


def main() -> None:
    args = parse_args()
    require_exists(args.dataset_dir)
    require_exists(args.model)
    data_yaml = args.dataset_dir / "dataset.yaml"
    require_exists(data_yaml)

    config_dir = Path("/tmp/Ultralytics")
    config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("YOLO_CONFIG_DIR", str(config_dir))

    from ultralytics import YOLO

    patch_ultralytics_sandbox_pools()

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    device = choose_device(args.device)
    model = YOLO(str(args.model))
    results = model.train(
        data=str(data_yaml),
        task="segment",
        imgsz=args.imgsz,
        epochs=args.epochs,
        batch=args.batch,
        workers=args.workers,
        device=device,
        patience=args.patience,
        seed=args.seed,
        project=str(args.output_dir.parent),
        name=args.output_dir.name,
        exist_ok=True,
        pretrained=True,
        cache=False,
    )

    save_dir = Path(results.save_dir)
    summary = {
        "dataset_dir": str(args.dataset_dir.resolve()),
        "data_yaml": str(data_yaml.resolve()),
        "base_model": str(args.model.resolve()),
        "device": device,
        "epochs": int(args.epochs),
        "imgsz": int(args.imgsz),
        "batch": int(args.batch),
        "workers": int(args.workers),
        "save_dir": str(save_dir.resolve()),
        "best_weights": str((save_dir / "weights" / "best.pt").resolve()),
        "last_weights": str((save_dir / "weights" / "last.pt").resolve()),
    }
    summary.update(best_metrics(save_dir))
    (save_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
