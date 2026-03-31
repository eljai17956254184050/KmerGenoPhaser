#!/usr/bin/env python3
"""
window_to_spectral_features_v2.py  ─  KmerGenoPhaser v1.1
==========================================================
Chromosome-level spectral feature extraction using CORRECTED complex encoding.

Fixes the encoding bug in the original window_to_spectral_features.py where
np.frombuffer(..., dtype=np.int8) was used, mapping bases to their raw ASCII
values (A=65, C=67, G=71, T=84).  This script replaces that with the
biologically meaningful vari-code complex encoding:

  A =  1+1j   (purine   + amino-group)
  G =  1-1j   (purine   + keto-group)
  C = -1+1j   (pyrimidine + amino-group)
  T = -1-1j   (pyrimidine + keto-group)
  N / other → 0+0j

The orthogonal axes preserve two independent biological distinctions in the
complex plane, leading to FFT peaks that reflect genuine compositional and
structural periodicity rather than ASCII numerical artefacts.

Use case
--------
Chromosome-level ancestry inference (low-resolution mode in unsupervised
pipeline).  Each sliding window of --window_size bp becomes one feature
vector.  Output is compatible with train_adaptive_unsupervised.py.

Output
------
Pickle dict:
  { window_id: np.ndarray(shape=(fft_size,), dtype=float64), ... }

  window_id format:  {SeqID}_{window_index:05d}
  e.g.  Chr1G_00001, Chr1G_00002, …

Usage
-----
python window_to_spectral_features_v2.py \\
    --input_fasta   target.fasta \\
    --output_pickle spectral_features.pkl \\
    --window_size   10000 \\
    --fft_size      1024

Author : Gengrui Zhu
Version: 1.1 (2026-03-31)  — Fixed ASCII→complex encoding
"""

import argparse
import os
import pickle
import sys

import numpy as np
from Bio import SeqIO
from scipy.fftpack import fft


# ──────────────────────────────────────────────────────────────────────────────
# Complex (vari-code) encoding  — identical to extract_block_features_fft.py
# ──────────────────────────────────────────────────────────────────────────────

_COMPLEX_MAP: dict = {
    'A':  1 + 1j,
    'C': -1 + 1j,
    'G':  1 - 1j,
    'T': -1 - 1j,
}


def encode_complex(seq: str) -> np.ndarray:
    """Map each ACGT base to its complex vari-code value; unknowns → 0+0j."""
    return np.array(
        [_COMPLEX_MAP.get(b, 0 + 0j) for b in seq.upper()],
        dtype=np.complex128,
    )


def compute_fft_features(seq: str, fft_size: int = 1024) -> np.ndarray:
    """
    Complex-encode → FFT → magnitude spectrum.

    Sequence truncated to `fft_size` or zero-padded; output dim = fft_size.
    """
    z = encode_complex(seq)
    if len(z) >= fft_size:
        z = z[:fft_size]
    else:
        z = np.pad(z, (0, fft_size - len(z)))
    return np.abs(fft(z)).astype(np.float64)


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            'Chromosome-level spectral feature extraction with complex encoding  '
            '(KmerGenoPhaser v1.1)'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Complex encoding:  A=1+j  C=-1+j  G=1-j  T=-1-j
OUTPUT INPUT_DIM = FFT_SIZE  (default 1024)
""",
    )
    parser.add_argument('--input_fasta',   required=True,
                        help='Target genome FASTA file')
    parser.add_argument('--output_pickle', required=True,
                        help='Output pickle path')
    parser.add_argument('--window_size',   type=int, default=10000,
                        help='Sliding window size in bp (default: 10000)')
    parser.add_argument('--fft_size',      type=int, default=1024,
                        help='FFT points / output spectrum dimension (default: 1024)')
    parser.add_argument('--target_chroms', nargs='+', default=None,
                        help='Only process these sequence IDs (default: all)')
    args = parser.parse_args()

    print("=" * 68)
    print("  window_to_spectral_features_v2.py  (KmerGenoPhaser v1.1)")
    print("=" * 68)
    print(f"  Window size: {args.window_size:,} bp")
    print(f"  FFT size   : {args.fft_size}")
    print(f"  INPUT_DIM  : {args.fft_size}")
    print(f"  Encoding   : complex vari-code  A=1+j C=-1+j G=1-j T=-1-j")
    print("=" * 68)

    target_set = set(args.target_chroms) if args.target_chroms else None

    window_data: dict = {}
    total_windows = 0
    skipped_seqs  = 0

    for record in SeqIO.parse(args.input_fasta, 'fasta'):
        if target_set and record.id not in target_set:
            continue

        seq_str  = str(record.seq)
        seq_len  = len(seq_str)
        n_windows = seq_len // args.window_size

        if n_windows == 0:
            print(f"  [SKIP] {record.id}: length {seq_len:,} < window_size "
                  f"{args.window_size:,}")
            skipped_seqs += 1
            continue

        print(f"  {record.id}: {seq_len:,} bp → {n_windows} window(s)")

        for i in range(n_windows):
            start  = i * args.window_size
            subseq = seq_str[start: start + args.window_size]

            # skip windows with no ACGT content
            valid = sum(1 for b in subseq.upper() if b in 'ACGT')
            if valid == 0:
                continue

            feat = compute_fft_features(subseq, args.fft_size)

            window_id = f"{record.id}_{i + 1:05d}"
            window_data[window_id] = feat
            total_windows += 1

    # ── save ─────────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(os.path.abspath(args.output_pickle)), exist_ok=True)
    with open(args.output_pickle, 'wb') as fh:
        pickle.dump(window_data, fh)

    print()
    print("=" * 68)
    print(f"  Total windows saved: {total_windows}")
    print(f"  Skipped sequences  : {skipped_seqs}")
    print(f"  Feature dim        : {args.fft_size}")
    print(f"  OUTPUT_PICKLE      : {args.output_pickle}")
    print()
    print(f">>> Set INPUT_DIM={args.fft_size} in your training script / conf <<<")
    print("=" * 68)


if __name__ == '__main__':
    main()
