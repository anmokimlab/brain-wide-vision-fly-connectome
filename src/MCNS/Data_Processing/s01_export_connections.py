#!/usr/bin/env python3
"""
Export all neuron-to-neuron connections from a neuprint dataset as a CSV file.

Output: MCNS/MCNS_Data/male-cns-v0.9-all-connections.csv
Output columns:
    pre_root_id, post_root_id, neuropil, syn_count, nt_type

By default this script uses the dataset:
    male-cns:v0.9

Requirements:
    - neuprint-python (pip install neuprint-python ; https://github.com/connectome-neuprint/neuprint-python)
    - a neuprint auth token: set NEUPRINT_APPLICATION_CREDENTIALS, or place the token in
      MCNS/Data_Processing/.neuprint_token (get it from https://neuprint.janelia.org -> Account).

Notes:
    - pre_root_id/post_root_id are filled from neuprint's bodyId values.
    - neuropil is filled from the ROI name.
    - syn_count is the per-ROI connection weight.
    - nt_type is filled from the presynaptic neuron's predictedNt.
    - By default only primary ROIs are exported, plus "NotPrimary".
      If you want overlapping non-primary ROIs too, pass --include-nonprimary.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

import ujson
from tqdm import tqdm

from neuprint import Client


DEFAULT_SERVER = "neuprint.janelia.org"
DEFAULT_DATASET = "male-cns:v0.9"
BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR.parent / "MCNS_Data"          # MCNS/MCNS_Data
DEFAULT_OUTPUT = DATA_DIR / "male-cns-v0.9-all-connections.csv"
DEFAULT_TOKEN_FILE = BASE_DIR / ".neuprint_token"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export all neuprint connections with ROI and predicted neurotransmitter."
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
        "--batch-size",
        type=int,
        default=100,
        help="batch size for neuprint adjacency queries (default: 100)",
    )
    parser.add_argument(
        "--include-nonprimary",
        action="store_true",
        help="include overlapping non-primary ROIs instead of primary ROIs only",
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


def fetch_neuron_nt_map(client: Client) -> dict[int, str | None]:
    q = """
        MATCH (n:Neuron)
        RETURN n.bodyId as bodyId, n.predictedNt as predictedNt
        ORDER BY n.bodyId
    """
    neurons_df = client.fetch_custom(q)

    if "predictedNt" not in neurons_df.columns:
        raise RuntimeError(
            "The dataset did not return a 'predictedNt' column. "
            "Please confirm that this property exists in the selected dataset."
        )

    return dict(zip(neurons_df["bodyId"], neurons_df["predictedNt"]))


def iter_connection_rows(client: Client, source_bodies: list[int]):
    q = f"""
        MATCH (n:Neuron)-[e:ConnectsTo]->(m:Neuron)
        WHERE n.bodyId in {source_bodies}
        RETURN n.bodyId as bodyId_pre,
               m.bodyId as bodyId_post,
               e.weight as weight,
               e.roiInfo as roiInfo
        ORDER BY n.bodyId, m.bodyId
    """
    return client.fetch_custom(q)


def write_batch_rows(
    writer: csv.writer,
    connections_df,
    nt_map: dict[int, str | None],
    primary_rois: set[str],
    include_nonprimary: bool,
) -> int:
    rows_written = 0

    for row in connections_df.itertuples(index=False):
        roi_info = ujson.loads(row.roiInfo)
        nt_type = nt_map.get(row.bodyId_pre)

        if include_nonprimary:
            for roi, weights in roi_info.items():
                syn_count = int(weights.get("post", 0))
                if syn_count < 1:
                    continue
                writer.writerow(
                    [row.bodyId_pre, row.bodyId_post, roi, syn_count, nt_type]
                )
                rows_written += 1
            continue

        primary_total = 0
        for roi, weights in roi_info.items():
            if roi not in primary_rois:
                continue

            syn_count = int(weights.get("post", 0))
            if syn_count < 1:
                continue

            primary_total += syn_count
            writer.writerow(
                [row.bodyId_pre, row.bodyId_post, roi, syn_count, nt_type]
            )
            rows_written += 1

        not_primary_count = int(row.weight) - primary_total
        if not_primary_count > 0:
            writer.writerow(
                [
                    row.bodyId_pre,
                    row.bodyId_post,
                    "NotPrimary",
                    not_primary_count,
                    nt_type,
                ]
            )
            rows_written += 1

    return rows_written


def export_connections(
    client: Client,
    output_path: str,
    batch_size: int,
    include_nonprimary: bool,
) -> int:
    nt_map = fetch_neuron_nt_map(client)
    source_bodies = sorted(nt_map.keys())
    primary_rois = set(client.primary_rois)

    output_file = Path(output_path)
    temp_output = output_file.with_suffix(output_file.suffix + ".part")
    temp_output.parent.mkdir(parents=True, exist_ok=True)

    total_rows = 0
    with temp_output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["pre_root_id", "post_root_id", "neuropil", "syn_count", "nt_type"]
        )

        for start in tqdm(range(0, len(source_bodies), batch_size)):
            batch_bodies = source_bodies[start:start + batch_size]
            connections_df = iter_connection_rows(client, batch_bodies)
            total_rows += write_batch_rows(
                writer=writer,
                connections_df=connections_df,
                nt_map=nt_map,
                primary_rois=primary_rois,
                include_nonprimary=include_nonprimary,
            )

    temp_output.replace(output_file)
    return total_rows


def main() -> int:
    args = parse_args()

    try:
        client = build_client(args.server, args.dataset, args.token_file)
        row_count = export_connections(
            client=client,
            output_path=args.output,
            batch_size=args.batch_size,
            include_nonprimary=args.include_nonprimary,
        )
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote {row_count:,} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
