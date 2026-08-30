#!/usr/bin/env python3
"""Write the sample IDs present in a BIOM HDF5 table.

Step 02 uses this file to sample only community-cluster samples that can be
extracted from the BIOM abundance table in step 03.
"""

import argparse
import os
import tempfile

import h5py


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract /sample/ids from a BIOM HDF5 table."
    )
    parser.add_argument("--biom", required=True, help="Input BIOM HDF5 table")
    parser.add_argument("--output", required=True, help="Output text file")
    return parser.parse_args()


def main():
    args = parse_args()

    if not os.path.isfile(args.biom):
        raise SystemExit(f"No existe la tabla BIOM de entrada: {args.biom}")

    output_dir = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(output_dir, exist_ok=True)

    ids = []
    with h5py.File(args.biom, "r") as handle:
        if "sample" not in handle or "ids" not in handle["sample"]:
            raise SystemExit(f"El BIOM no contiene el dataset /sample/ids: {args.biom}")

        for raw_id in handle["sample"]["ids"][()]:
            sample_id = raw_id.decode("utf-8") if isinstance(raw_id, bytes) else str(raw_id)
            sample_id = sample_id.strip()
            if sample_id:
                ids.append(sample_id)

    if not ids:
        raise SystemExit(f"No se encontraron IDs de muestra en el BIOM: {args.biom}")

    unique_ids = set(ids)
    if len(unique_ids) != len(ids):
        raise SystemExit(
            "El BIOM contiene IDs de muestra duplicados; no se puede muestrear "
            f"de forma segura: {args.biom}"
        )

    fd, tmp_path = tempfile.mkstemp(
        prefix=".biom_sample_ids.",
        suffix=".tmp",
        dir=output_dir,
        text=True,
    )
    try:
        with os.fdopen(fd, "w") as out:
            for sample_id in ids:
                out.write(sample_id + "\n")
        os.replace(tmp_path, args.output)
    except Exception:
        try:
            os.unlink(tmp_path)
        except FileNotFoundError:
            pass
        raise

    print(f"BIOM sample IDs written: {args.output}")
    print(f"BIOM sample IDs count: {len(ids)}")


if __name__ == "__main__":
    main()
