#!/usr/bin/env python3
"""
Export a male-cns:v0.9 classification table as CSV.

Output: MCNS/MCNS_Data/male-cns-v0.9-classification.csv
Output columns:
    root_id, super_class, side

The super_class values are normalized as follows:
    ascending_neuron  -> ascending
    descending_neuron -> descending
    cb_intrinsic      -> central
    ol_intrinsic      -> optic

Requirements:
    - neuprint-python (pip install neuprint-python ; https://github.com/connectome-neuprint/neuprint-python)
    - a neuprint auth token: set NEUPRINT_APPLICATION_CREDENTIALS, or place the token in
      MCNS/Data_Processing/.neuprint_token.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from neuprint import Client


DEFAULT_SERVER = "neuprint.janelia.org"
DEFAULT_DATASET = "male-cns:v0.9"
BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR.parent / "MCNS_Data"          # MCNS/MCNS_Data
DEFAULT_OUTPUT = DATA_DIR / "male-cns-v0.9-classification.csv"
DEFAULT_TOKEN_FILE = BASE_DIR / ".neuprint_token"

SUPERCLASS_RENAMES = {
    "ascending_neuron": "ascending",
    "descending_neuron": "descending",
    "cb_intrinsic": "central",
    "ol_intrinsic": "optic",
}

SIDE_RENAMES = {
    "R": "right",
    "L": "left",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export male-cns-v0.9-classification.csv from male-cns:v0.9."
    )
    parser.add_argument(
        "--server",
        default=DEFAULT_SERVER,
        help=f"neuprint server hostname (default: {DEFAULT_SERVER})",
    )
    parser.add_argument(
        "--dataset",
        default=DEFAULT_DATASET,
        help=f"neuprint dataset name (default: {DEFAULT_DATASET})",
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_OUTPUT),
        help=f"output CSV path (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--token-file",
        default=str(DEFAULT_TOKEN_FILE),
        help=f"token file path (default: {DEFAULT_TOKEN_FILE})",
    )
    return parser.parse_args()


def load_token(token_file: str) -> str:
    token = os.environ.get("NEUPRINT_APPLICATION_CREDENTIALS")
    if token:
        return token

    token_path = Path(token_file)
    if token_path.exists():
        token = token_path.read_text(encoding="utf-8").strip()
        if token:
            return token

    raise RuntimeError(
        "No neuprint token found. Set NEUPRINT_APPLICATION_CREDENTIALS "
        f"or place the token in {token_path}."
    )


def build_client(server: str, dataset: str, token_file: str) -> Client:
    token = load_token(token_file)
    return Client(server, dataset=dataset, token=token)


def fetch_classification_table(client: Client):
    q = """
        MATCH (n:Neuron)
        RETURN n.bodyId AS root_id,
               n.superclass AS super_class,
               n.somaSide AS side
        ORDER BY root_id
    """
    df = client.fetch_custom(q)
    df["super_class"] = df["super_class"].replace(SUPERCLASS_RENAMES)
    df["side"] = df["side"].replace(SIDE_RENAMES)
    return df


def export_classification_table(client: Client, output_path: str) -> Path:
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    df = fetch_classification_table(client)
    df.to_csv(output_file, index=False)
    return output_file


def main() -> int:
    args = parse_args()
    client = build_client(args.server, args.dataset, args.token_file)
    output_file = export_classification_table(client, args.output)
    print(f"Wrote {output_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
