#!/usr/bin/env python3
"""
Export optic-lobe ROI meshes from a neuprint dataset as per-ROI vertices/faces CSV files.

This is the MCNS analogue of the FAFB Data_Processing/s01_fetch_neuropil_neuron.ipynb mesh
fetch. The meshes are read by MCNS/Data_Processing/s11_FBP_output_layer.m.

For each ROI, this script writes:

    MCNS/Processed_Data/roi_mesh_csv/<ROI>/vertices.csv   (columns: x, y, z)
    MCNS/Processed_Data/roi_mesh_csv/<ROI>/faces.csv      (columns: v1, v2, v3)

Vertices are in raw neuprint units (8 nm voxels); face indices are zero-based.

Requirements:
    - neuprint-python (pip install neuprint-python ; https://github.com/connectome-neuprint/neuprint-python).
    - a neuprint auth token: set NEUPRINT_APPLICATION_CREDENTIALS, or place the token in
      MCNS/Data_Processing/.neuprint_token.
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

from neuprint import Client


DEFAULT_SERVER = "neuprint.janelia.org"
DEFAULT_DATASET = "male-cns:v0.9"
BASE_DIR = Path(__file__).resolve().parent
PROCESSED_DIR = BASE_DIR.parent / "Processed_Data"          # MCNS/Processed_Data
DEFAULT_OUTPUT_DIR = PROCESSED_DIR / "roi_mesh_csv"
DEFAULT_TOKEN_FILE = BASE_DIR / ".neuprint_token"
DEFAULT_ROIS = [
    "LOP(R)",
    "LOP(L)",
    "LO(R)",
    "LO(L)",
    "ME(R)",
    "ME(L)",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export selected male-cns:v0.9 ROI meshes as vertices/faces CSV files."
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
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help=f"output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--token-file",
        default=str(DEFAULT_TOKEN_FILE),
        help=f"token file path (default: {DEFAULT_TOKEN_FILE})",
    )
    parser.add_argument(
        "--rois",
        nargs="+",
        default=DEFAULT_ROIS,
        help="ROI names to export",
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


def _parse_face_index(token: str) -> int:
    vertex_token = token.split("/")[0]
    return int(vertex_token) - 1


def parse_obj_mesh(obj_bytes: bytes) -> tuple[list[list[float]], list[list[int]]]:
    vertices: list[list[float]] = []
    faces: list[list[int]] = []

    for raw_line in obj_bytes.decode("utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        if line.startswith("v "):
            parts = line.split()
            vertices.append([float(parts[1]), float(parts[2]), float(parts[3])])
            continue

        if line.startswith("f "):
            parts = line.split()[1:]
            indices = [_parse_face_index(part) for part in parts]

            if len(indices) < 3:
                continue

            if len(indices) == 3:
                faces.append(indices)
                continue

            # Triangulate polygon faces with a simple fan.
            for i in range(1, len(indices) - 1):
                faces.append([indices[0], indices[i], indices[i + 1]])

    return vertices, faces


def write_vertices_csv(path: Path, vertices: list[list[float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["x", "y", "z"])
        writer.writerows(vertices)


def write_faces_csv(path: Path, faces: list[list[int]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["v1", "v2", "v3"])
        writer.writerows(faces)


def export_roi_mesh_csvs(client: Client, rois: list[str], output_dir: str) -> list[Path]:
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)

    written_paths: list[Path] = []
    for roi in rois:
        obj_bytes = client.fetch_roi_mesh(roi)
        vertices, faces = parse_obj_mesh(obj_bytes)

        roi_dir = output_root / roi
        roi_dir.mkdir(parents=True, exist_ok=True)

        vertices_path = roi_dir / "vertices.csv"
        faces_path = roi_dir / "faces.csv"

        write_vertices_csv(vertices_path, vertices)
        write_faces_csv(faces_path, faces)

        written_paths.extend([vertices_path, faces_path])
        print(
            f"Wrote {roi}: {len(vertices)} vertices, {len(faces)} faces -> {roi_dir}"
        )

    return written_paths


def main() -> int:
    args = parse_args()
    client = build_client(args.server, args.dataset, args.token_file)
    export_roi_mesh_csvs(client, args.rois, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
