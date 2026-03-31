# Changelog — KmerGenoPhaser

All notable changes to this project will be documented here.

---

## [v1.1] — 2026-03-31

### Summary

This release introduces **selectable feature encoding** for the `unsupervised`
module, fixing a latent encoding bug and significantly extending the feature
space available to the autoencoder.  All other modules (`supervised`, `snpml`,
`karyotype`) are unchanged.

### New Files

| File | Description |
|---|---|
| `lib/unsupervised/extract_block_features_fft.py` | New unified block-level feature extractor; replaces `extract_block_features.py` as the default in the pipeline |
| `lib/unsupervised/window_to_spectral_features_v2.py` | Fixed chromosome-window spectral extractor (see Bug Fixes below) |

### Bug Fixes

**`window_to_spectral_features.py` — Wrong FFT encoding (ASCII vs. complex)**

The original script encoded DNA bases using raw ASCII integer values via
`np.frombuffer(seq.encode(), dtype=np.int8)`, mapping A→65, C→67, G→71,
T→84.  These values carry no biological meaning and produce FFT spectra that
reflect ASCII table distances rather than sequence composition.

The corrected encoding (`window_to_spectral_features_v2.py` and the new
`extract_block_features_fft.py`) uses the vari-code complex mapping:

```
A =  1+1j   (purine   + amino-group)
G =  1-1j   (purine   + keto-group)
C = -1+1j   (pyrimidine + amino-group)
T = -1-1j   (pyrimidine + keto-group)
N / other → 0+0j
```

The real axis encodes purine/pyrimidine identity and the imaginary axis
encodes amino/keto identity, keeping both distinctions orthogonal in the
complex plane.  Empirical testing on sugarcane Chr1 data showed a reduction
in subgenome-inconsistent blocks from 94 → 48 after switching to this
encoding (Silhouette score maintained at 0.97).

### New Features

#### `--encoding` — selectable feature extraction strategy

`KmerGenoPhaser unsupervised` and `extract_block_features_fft.py` now accept
a `--encoding` argument:

| Value | Description | `INPUT_DIM` (k=1..5, fft=1024) |
|---|---|---|
| `kmer` | K-mer frequency only (original v1.0 behaviour) | 1364 |
| `fft` | Complex FFT magnitude spectrum only | 1024 |
| `concat` | K-mer + complex FFT concatenated **[new default]** | 2388 |

#### `--feature_mode` — block vs. genome resolution

| Value | Extractor called | Steps enabled |
|---|---|---|
| `block` | `extract_block_features_fft.py` | All (1-5) |
| `genome` | `window_to_spectral_features_v2.py` | 1-2 only |

#### `--fft_size` — configurable FFT window

Default 1024.  Controls both the number of FFT points applied to each
sequence and the resulting feature dimension for the FFT component.

#### `INPUT_DIM` is now auto-computed

`KmerGenoPhaser_unsupervised.sh` no longer requires `INPUT_DIM` to be
manually set in `conf/kmergenophaser.conf`.  It is computed at runtime from
`--encoding`, `--min_kmer`, `--max_kmer`, and `--fft_size` and passed
directly to `train_adaptive_unsupervised.py`.

### Changes to `conf/kmergenophaser.conf`

Three new defaults (all overridable via CLI):

```bash
FEATURE_MODE=block          # block | genome
ENCODING=concat             # kmer  | fft | concat
FFT_SIZE=1024
GENOME_WINDOW_SIZE=10000    # window size for genome mode
# INPUT_DIM is now computed at runtime — no need to edit this
```

### Upgrade Guide

If you are running the v1.0 pipeline and want identical behaviour:

```bash
KmerGenoPhaser unsupervised \
    ... \
    --encoding kmer       # restore k-mer only extraction
```

If you want the new default (concat):

1. Delete any cached `.pkl` and `*_block_distances.tsv` files from previous
   runs — the pipeline detects existing files and may skip re-extraction/training.
2. Run with `--encoding concat` (or omit `--encoding`, it is the new default).
3. `INPUT_DIM` will be reported in the log; you do **not** need to update
   `conf/kmergenophaser.conf`.

---

## [v1.0] — 2025-03  (initial public release)

First public release of KmerGenoPhaser.

### Modules

- **supervised** — K-mer specificity scoring + genome mapping
- **snpml** — SNP + Maximum-Likelihood ancestry block calling
- **unsupervised** — Autoencoder-based ancestry block discovery
- **karyotype** — Idiogram-style visualization
- **build-inputs** — Metadata file generator

### Unsupervised module feature extraction (v1.0)

- Block-level: `extract_block_features.py` — k-mer frequency only (k=1..5,
  dim=1364)
- Genome-level: `window_to_spectral_features.py` — ASCII FFT (has encoding
  bug fixed in v1.1)
- `INPUT_DIM` must be set manually in `conf/kmergenophaser.conf`

---

## Citation

> [Paper citation placeholder]

## License

This project is licensed under a **Non-Commercial Research License**.  
Free to use for academic and research purposes only. Commercial use is strictly prohibited.  
See the [LICENSE](./LICENSE) file for details.
