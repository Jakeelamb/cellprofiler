#!/usr/bin/env python3
"""Train a Cellpose model while honoring training patch size (`bsize`)."""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

from cellpose import io, models, train
import torch


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dir", type=Path, required=True)
    parser.add_argument("--test-dir", type=Path, default=None)
    parser.add_argument("--mask-filter", default="_masks")
    parser.add_argument("--pretrained-model", default="cpsam")
    parser.add_argument("--model-name-out", required=True)
    parser.add_argument("--learning-rate", type=float, default=1e-5)
    parser.add_argument("--weight-decay", type=float, default=0.1)
    parser.add_argument("--n-epochs", type=int, default=100)
    parser.add_argument("--train-batch-size", type=int, default=1)
    parser.add_argument("--bsize", type=int, default=256)
    parser.add_argument("--nimg-per-epoch", type=int, default=0)
    parser.add_argument("--nimg-test-per-epoch", type=int, default=0)
    parser.add_argument("--min-train-masks", type=int, default=5)
    parser.add_argument("--save-every", type=int, default=10)
    parser.add_argument("--save-each", action="store_true")
    parser.add_argument("--look-one-level-down", action="store_true")
    parser.add_argument("--img-filter", default="")
    parser.add_argument("--use-gpu", action="store_true")
    parser.add_argument("--gpu-device", default=None)
    parser.add_argument("--channel-axis", type=int, default=None)
    parser.add_argument("--no-norm", action="store_true")
    parser.add_argument("--keep-bfloat16", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def install_dtype_guard(net: object, enabled: bool) -> callable:
    if not enabled:
        return lambda: None
    cls = type(net)
    dtype_prop = cls.dtype

    def patched_setter(self, value):
        if (
            getattr(self, "_keep_bfloat16_training", False)
            and dtype_prop.fget(self) == torch.bfloat16
            and value == torch.float32
        ):
            return
        dtype_prop.fset(self, value)

    cls.dtype = property(dtype_prop.fget, patched_setter)
    setattr(net, "_keep_bfloat16_training", True)

    def restore() -> None:
        setattr(net, "_keep_bfloat16_training", False)
        cls.dtype = dtype_prop

    return restore


def ensure_channel_axis(images: list) -> list:
    fixed = []
    for image in images:
        if getattr(image, "ndim", 0) == 2:
            fixed.append(image[..., None])
        else:
            fixed.append(image)
    return fixed


def infer_channel_axis(images: list, requested: int | None) -> int | None:
    if requested is not None:
        return requested
    for image in images:
        if getattr(image, "ndim", 0) == 3 and image.shape[-1] <= 4 and image.shape[0] > 4:
            return -1
    return requested


def main() -> None:
    args = parse_args()
    if args.verbose:
        logger, _ = io.logger_setup()
    else:
        logger = logging.getLogger(__name__)
    image_filter = args.img_filter or None
    test_dir = None if args.test_dir is None else os.fspath(args.test_dir)
    device, _ = models.assign_device(use_torch=True, gpu=args.use_gpu, device=args.gpu_device)
    normalize = not args.no_norm

    images, labels, image_names, test_images, test_labels, image_names_test = io.load_train_test_data(
        os.fspath(args.dir),
        test_dir,
        image_filter,
        args.mask_filter,
        args.look_one_level_down,
    )
    images = ensure_channel_axis(images)
    test_images = ensure_channel_axis(test_images)
    channel_axis = infer_channel_axis(images + test_images, args.channel_axis)
    model = models.CellposeModel(device=device, pretrained_model=args.pretrained_model)
    restore_dtype = install_dtype_guard(model.net, args.keep_bfloat16)
    try:
        cpmodel_path = train.train_seg(
            model.net,
            images,
            labels,
            train_files=image_names,
            test_data=test_images,
            test_labels=test_labels,
            test_files=image_names_test,
            compute_flows=False,
            load_files=True,
            normalize=normalize,
            channel_axis=channel_axis,
            learning_rate=args.learning_rate,
            weight_decay=args.weight_decay,
            n_epochs=args.n_epochs,
            batch_size=args.train_batch_size,
            min_train_masks=args.min_train_masks,
            nimg_per_epoch=None if args.nimg_per_epoch <= 0 else args.nimg_per_epoch,
            nimg_test_per_epoch=None if args.nimg_test_per_epoch <= 0 else args.nimg_test_per_epoch,
            save_path=os.path.realpath(args.dir),
            save_every=args.save_every,
            save_each=args.save_each,
            model_name=args.model_name_out,
            bsize=args.bsize,
        )[0]
    finally:
        restore_dtype()
    logger.info(">>>> model trained and saved to %s", cpmodel_path)


if __name__ == "__main__":
    main()
