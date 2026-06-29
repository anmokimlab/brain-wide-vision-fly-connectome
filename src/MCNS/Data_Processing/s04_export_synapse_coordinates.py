#!/usr/bin/env python3
"""
Export all neuron-to-neuron synapse connections from a neuprint dataset as a CSV file.

Output: MCNS/MCNS_Data/male-cns-v0.9-synapse-coordinates.csv
Output columns:
    pre_x, pre_y, pre_z, pre_root_id, post_root_id, neuropil

One row is written per synapse connection, not per neuron-neuron pair.

Requirements:
    - neuprint-python (pip install neuprint-python ; https://github.com/connectome-neuprint/neuprint-python).
      This script also uses neuprint's internal query helpers
      (neuprint.queries.synapses._fetch_synapse_connections,
       neuprint.queries.neuroncriteria.NeuronCriteria), so a recent neuprint-python is required.
    - a neuprint auth token: set NEUPRINT_APPLICATION_CREDENTIALS, or place the token in
      MCNS/Data_Processing/.neuprint_token.

The run is resumable: progress is tracked in a <output>.progress.json file and partial
rows are written to <output>.part, so an interrupted export can continue where it stopped.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path

from neuprint import Client, fetch_custom, SynapseCriteria as SC
from neuprint.queries.synapses import _fetch_synapse_connections
from neuprint.queries.neuroncriteria import NeuronCriteria as NC
from tqdm import tqdm


DEFAULT_SERVER = "neuprint.janelia.org"
DEFAULT_DATASET = "male-cns:v0.9"
BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR.parent / "MCNS_Data"          # MCNS/MCNS_Data
DEFAULT_OUTPUT = DATA_DIR / "male-cns-v0.9-synapse-coordinates.csv"
DEFAULT_TOKEN_FILE = BASE_DIR / ".neuprint_token"
RETRY_DELAY_SECONDS = 180
MAX_RETRIES = 3


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export all synapse-level neuron-to-neuron coordinates from a neuprint dataset."
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
    parser.add_argument(
        "--batch-size",
        type=int,
        default=25,
        help="number of presynaptic neurons to query per batch (default: 25)",
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


def fetch_neuron_ids(client: Client) -> list[int]:
    q = """
        MATCH (n:Neuron)
        RETURN n.bodyId as bodyId
        ORDER BY n.bodyId
    """
    df = fetch_custom(q, client=client)
    return df["bodyId"].astype(int).tolist()


def load_progress(progress_file: Path) -> dict:
    if not progress_file.exists():
        return {
            "next_start": 0,
            "rows_written": 0,
            "batch_size": None,
        }

    with progress_file.open("r", encoding="utf-8") as f:
        return json.load(f)


def save_progress(progress_file: Path, next_start: int, rows_written: int, batch_size: int) -> None:
    progress = {
        "next_start": next_start,
        "rows_written": rows_written,
        "batch_size": batch_size,
        "updated_at": int(time.time()),
    }
    progress_file.write_text(json.dumps(progress, indent=2), encoding="utf-8")


def fetch_batch_with_retries(
    client: Client,
    batch_ids: list[int],
    synapse_criteria: SC,
) -> object:
    last_error = None

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            src_crit = NC(bodyId=batch_ids, label="Neuron", client=client)
            tgt_crit = NC(label="Neuron", client=client)
            return _fetch_synapse_connections(
                source_criteria=src_crit,
                target_criteria=tgt_crit,
                synapse_criteria=synapse_criteria,
                min_total_weight=1,
                nt=None,
                client=client,
            )
        except Exception as exc:
            last_error = exc
            if attempt == MAX_RETRIES:
                break

            print(
                f"Batch starting with bodyId {batch_ids[0]} failed on attempt {attempt}/{MAX_RETRIES}. "
                f"Retrying in {RETRY_DELAY_SECONDS} seconds...",
                file=sys.stderr,
                flush=True,
            )
            time.sleep(RETRY_DELAY_SECONDS)

    raise last_error


def export_synapse_coordinates(client: Client, output_path: str, batch_size: int) -> int:
    neuron_ids = fetch_neuron_ids(client)
    synapse_criteria = SC(primary_only=True, client=client)

    output_file = Path(output_path)
    temp_output = output_file.with_suffix(output_file.suffix + ".part")
    progress_file = output_file.with_suffix(output_file.suffix + ".progress.json")
    temp_output.parent.mkdir(parents=True, exist_ok=True)

    progress = load_progress(progress_file)
    next_start = int(progress.get("next_start", 0))
    total_rows = int(progress.get("rows_written", 0))
    previous_batch_size = progress.get("batch_size")

    if previous_batch_size not in (None, batch_size):
        raise RuntimeError(
            f"Existing progress file uses batch_size={previous_batch_size}, "
            f"but current run uses batch_size={batch_size}. "
            "Please keep the same batch size while resuming, or remove the .part and .progress.json files to restart."
        )

    if next_start > 0 and not temp_output.exists():
        raise RuntimeError(
            f"Found progress file {progress_file}, but partial output file {temp_output} is missing."
        )

    file_mode = "a" if temp_output.exists() and next_start > 0 else "w"
    with temp_output.open(file_mode, newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        if file_mode == "w":
            writer.writerow(
                ["pre_x", "pre_y", "pre_z", "pre_root_id", "post_root_id", "neuropil"]
            )
            save_progress(progress_file, next_start=0, rows_written=0, batch_size=batch_size)

        for start in tqdm(range(next_start, len(neuron_ids), batch_size), initial=next_start // batch_size):
            batch_ids = neuron_ids[start:start + batch_size]
            syn_df = fetch_batch_with_retries(
                client=client,
                batch_ids=batch_ids,
                synapse_criteria=synapse_criteria,
            )

            if len(syn_df) == 0:
                save_progress(
                    progress_file,
                    next_start=min(start + batch_size, len(neuron_ids)),
                    rows_written=total_rows,
                    batch_size=batch_size,
                )
                continue

            batch_rows = syn_df[
                ["x_pre", "y_pre", "z_pre", "bodyId_pre", "bodyId_post", "roi_post"]
            ].itertuples(index=False, name=None)

            for row in batch_rows:
                writer.writerow(row)
                total_rows += 1

            f.flush()
            save_progress(
                progress_file,
                next_start=min(start + batch_size, len(neuron_ids)),
                rows_written=total_rows,
                batch_size=batch_size,
            )

    temp_output.replace(output_file)
    progress_file.unlink(missing_ok=True)
    return total_rows


def main() -> int:
    args = parse_args()

    try:
        client = build_client(args.server, args.dataset, args.token_file)
        row_count = export_synapse_coordinates(
            client=client,
            output_path=args.output,
            batch_size=args.batch_size,
        )
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote {row_count:,} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
