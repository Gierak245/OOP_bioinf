# Sequence Toolkit

A modular bioinformatics toolkit for parsing FASTA/FASTQ files, computing k‑mer counts with on‑disk caching, and running arbitrarily extensible “analyzer” plugins (e.g. GC content, motif search).

---

## Features

### `Sequence` **class** with methods:

* from_fasta(text) / from_fastq(text) for single‑record parsing
* parse_fasta_file(filepath) / parse_fastq_file(filepath) for streaming multi‑record parsing
* expensive_kmer() decorated with @disk_cache to cache and load k‑mer computations (named "expensive" to reflect its computational cost)

### `disk_cache` decorator:

* Caches function outputs to disk using SHA‑256 keys
* Automatically creates cache directory
* Loads from cache on hit, computes & saves on miss

### **`AnalyzerMeta` metaclass**  
* Auto‑registers every `Analyzer` subclass in a central `Analyzer.registry`  

### **`Analyzer` (abstract base class)**  
* Defines the `run(self, seq: Sequence) -> Any` interface  

### **`Built‑in analyzers`**  
* `GCContentAnalyzer`  
* `MotifSearchAnalyzer(motif)`  
* (Easily add more by subclassing `Analyzer`!)  

### **`PluginManager`**  
* Discovers & instantiates plugins by name  
* Supports user‑specified order and dynamic additions at runtime  
* Runs all analyzers on one or many sequences, returning nested result dicts

---
### `CLI (seq_toolkit.py)`:
* -i/--input: FASTA or FASTQ file
*-o/--output: cache directory and location for per sentence and aggregated results (all_kmers.pkl)
* --clear_cache: archive & prune cache entries older than 7 days

## Installation

```
git clone https://github.com/Gierak245/OOP_bioinf
cd OOP_bioinf
```
*Note: Uses **only Python stdlib**; no external packages required.*

## Usage

### **As a library**
```
from seq_toolkit import Sequence
from Analyzer import PluginManager, GCContentAnalyzer, MotifSearchAnalyzer

# parse sequences
records = Sequence.parse_fasta_file("example.fa")

# configure and run analyzers
pm = PluginManager({
    "GCContentAnalyzer": {},
    "MotifSearchAnalyzer": {"motif": "ATG"}
})
pm.discover(["GCContentAnalyzer", "MotifSearchAnalyzer"])

for rec in records:
    print(rec.header, pm.run_instances(rec))
```
### CLI
```
python seqtools.py \
  -i path/to/input.fa \
  -o path/to/cache_dir \
  --clear_cache True
```
- Output directory will contain `all_kmers.pkl` and cached `.pkl` per sequence.
- On subsequent runs, cached results are reused automatically.

*Note: To test the function, example files in FASTA (ex_fasta.fa) and FASTQ (ex_fastq.fastq) format are provided.*
