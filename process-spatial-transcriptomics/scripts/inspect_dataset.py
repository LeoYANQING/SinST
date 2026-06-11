#!/usr/bin/env python3
"""Inspect an h5ad file without modifying it."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad

from pipeline_utils import audit_adata, json_safe, load_config, write_json


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit an AnnData spatial transcriptomics dataset."
    )
    parser.add_argument("input", type=Path, help="Input .h5ad file")
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument(
        "--profile",
        choices=("generic", "axolotl-stereoseq", "visium"),
        default=None,
    )
    parser.add_argument("--output", type=Path, default=None, help="Optional JSON output")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config, args.profile)
    adata = ad.read_h5ad(args.input, backed="r")
    try:
        audit = audit_adata(adata, config, args.input)
    finally:
        if adata.file is not None:
            adata.file.close()

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        write_json(audit, args.output)
    print(json.dumps(json_safe(audit), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
