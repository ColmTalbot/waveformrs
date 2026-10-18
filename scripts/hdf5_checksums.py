#!/usr/bin/env python
"""
Dump the shape, element type and a checksum of every dataset of an HDF5 file
with h5py, for testing the pure-Rust reader in ``src/hdf5.rs``
(``tests/nrsur.rs``).

The checksum of a dataset with values ``v_i`` (row-major, as 64-bit
patterns: ``float64`` bits or ``int64`` two's complement) is

    sum_i ((bits(v_i) ^ (i * K1)) * K2)   (mod 2**64)

with ``K1 = 0x9E3779B97F4A7C15`` and ``K2 = 0xBF58476D1CE4E5B9``.

Usage
-----
    python scripts/hdf5_checksums.py lal_data/NRSur7dq4_v1.0.h5 \
        tests/data/nrsur7dq4_hdf5_checksums.json
"""
import argparse
import hashlib
import json

import h5py
import numpy as np

K1 = np.uint64(0x9E3779B97F4A7C15)
K2 = np.uint64(0xBF58476D1CE4E5B9)


def checksum(values):
    if values.dtype == np.float64:
        bits = values.astype("<f8").view(np.uint64).ravel()
    elif values.dtype.kind == "i":
        bits = values.astype("<i8").view(np.uint64).ravel()
    else:
        raise ValueError(f"unsupported dtype {values.dtype}")
    index = np.arange(bits.size, dtype=np.uint64)
    with np.errstate(over="ignore"):
        return int(((bits ^ (index * K1)) * K2).sum(dtype=np.uint64))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path")
    parser.add_argument("output")
    args = parser.parse_args()
    datasets = dict()

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            values = obj[()]
            datasets[name] = dict(
                shape=list(obj.shape),
                dtype=values.dtype.str,
                checksum=f"{checksum(np.asarray(values)):016x}",
            )

    with h5py.File(args.path, "r") as ff:
        ff.visititems(visit)
    with open(args.path, "rb") as ff:
        md5 = hashlib.md5(ff.read()).hexdigest()
    with open(args.output, "w") as ff:
        json.dump(dict(md5=md5, h5py_version=h5py.__version__, datasets=datasets), ff)
    print(f"Wrote {len(datasets)} datasets to {args.output}")


if __name__ == "__main__":
    main()
