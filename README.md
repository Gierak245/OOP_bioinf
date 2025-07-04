## Overview

A lightweight, extensible bioinformatics toolkit for parsing multi‑record FASTA/FASTQ files and computing k‑mer counts with an efficient on‑disk caching mechanism. Includes a CLI for end‑to‑end processing,
cache archiving, and pruning of old entries.

## Features

### `Sequence` **class** with methods:

* from_fasta(text) / from_fastq(text) for single‑record parsing
* parse_fasta_file(filepath) / parse_fastq_file(filepath) for streaming multi‑record parsing
* expensive_kmer() decorated with @disk_cache to cache and load k‑mer computations (named "expensive" to reflect its computational cost)

### `disk_cache` decorator:

* Caches function outputs to disk using SHA‑256 keys
* Automatically creates cache directory
* Loads from cache on hit, computes & saves on miss

### `CLI (seq_toolkit.py)`:
* -i/--input: FASTA or FASTQ file
*-o/--output: cache directory and location for per sentence and aggregated results (all_kmers.pkl)
* --clear_cache: archive & prune cache entries older than 7 days

## Installation

1. Clone repository:
```
git clone https://github.com/Gierak245/OOP_bioinf
cd OOP_bioinf
```
*Note: Uses **only Python stdlib**; no external packages required.*
## Usage

```
python seqtools.py \
  -i path/to/input.fa \
  -o path/to/cache_dir \
  --clear_cache True
```
- Output directory will contain `all_kmers.pkl` and cached `.pkl` per sequence.
- On subsequent runs, cached results are reused automatically.
