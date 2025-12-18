# Per-Read Filtering Feature

## Overview

The fastq-filter tool has been extended to support **per-read threshold values** for paired-end and multi-read filtering. Previously, all reads in a tuple (e.g., R1 and R2 in paired-end sequencing) were evaluated against a single threshold. Now you can specify different thresholds for each read position.

## What Changed

### Backward Compatibility
- **Single threshold** behavior is unchanged - existing code and CLI commands work exactly as before
- All existing tests pass without modification

### New Feature: Per-Read Thresholds
- Filters now accept either:
  - A single value (original behavior)
  - A tuple/list of values (one per read position)

### Affected Filters
All four filter types now support per-read thresholds:

1. **`MinimumLengthFilter`** - Minimum read length
2. **`MaximumLengthFilter`** - Maximum read length  
3. **`AverageErrorRateFilter`** - Average error rate
4. **`MedianQualityFilter`** - Median quality score

## Usage

### Python API

```python
from fastq_filter import (
    MinimumLengthFilter,
    MedianQualityFilter,
    AverageErrorRateFilter
)
from dnaio import SequenceRecord

# Single threshold (original behavior)
single_filter = MinimumLengthFilter(50)
result = single_filter((r1, r2))  # Passes if ANY read >= 50bp

# Per-read thresholds (new feature)
per_read_filter = MinimumLengthFilter((80, 60))
result = per_read_filter((r1, r2))  # R1 must be ≥80bp AND R2 must be ≥60bp

# Works with quality filters too
quality_filter = MedianQualityFilter((35, 25))  # R1 ≥Q35, R2 ≥Q25
error_filter = AverageErrorRateFilter((0.001, 0.01))  # R1 ≤0.1%, R2 ≤1%
```

### Command Line Interface

Use **comma-separated values** to specify per-read thresholds:

```bash
# Single threshold (original)
fastq-filter -l 50 -o r1.fq -o r2.fq r1_in.fq r2_in.fq

# Per-read minimum length: R1 ≥ 80bp, R2 ≥ 60bp
fastq-filter -l 80,60 -o r1.fq -o r2.fq r1_in.fq r2_in.fq

# Per-read median quality: R1 ≥ Q30, R2 ≥ Q25
fastq-filter -Q 30,25 -o r1.fq -o r2.fq r1_in.fq r2_in.fq

# Combined filters with per-read thresholds
fastq-filter -l 80,60 -Q 30,25 -o r1.fq -o r2.fq r1_in.fq r2_in.fq

# Works with 3+ files (e.g., R1, R2, Index)
fastq-filter -l 90,90,6 -Q 30,30,35 -o r1.fq -o r2.fq -o i1.fq \\
  r1_in.fq r2_in.fq i1_in.fq
```

## Filtering Semantics

### Original Behavior (Single Threshold)
- **Length filters**: 
  - `MinimumLengthFilter`: Passes if **ANY** read meets threshold
  - `MaximumLengthFilter`: Fails if **ANY** read exceeds threshold
- **Quality filters**: Aggregate all bases across all reads, compute single metric

### New Behavior (Per-Read Thresholds)
- **All filters**: Each read must pass its respective threshold
- Evaluation is short-circuited: stops at first failing read
- More stringent: **ALL** reads must pass their individual thresholds

## Implementation Details

### C Extension Changes
- [src/fastq_filter/_filtersmodule.c](src/fastq_filter/_filtersmodule.c):
  - Added threshold array storage to `FastqFilter` struct
  - Constructors now accept sequences and allocate per-read threshold arrays
  - Filter `__call__` methods check array length matches number of reads
  - Custom getters return tuples when per-read thresholds are set
  - Memory management: arrays freed in `dealloc`

### Python Layer Changes
- [src/fastq_filter/__init__.py](src/fastq_filter/__init__.py):
  - Added `parse_threshold()` helper for comma-separated CLI values
  - Updated argument parser help text
  - CLI arguments now parse as strings and convert to int/float tuples
  - Logging displays comma-separated thresholds when applicable

## Testing

### Test Files
1. **`test_per_read_filters.py`** - Comprehensive unit tests showing:
   - Backward compatibility with single thresholds
   - Per-read quality filtering (average & median)
   - Per-read length filtering (min & max)
   - Triple-read filtering (R1, R2, Index)
   - Statistics tracking

2. **`demo_per_read_simple.py`** - CLI demo with realistic examples:
   - Creates sample paired-end FASTQ files
   - Demonstrates single vs. per-read thresholds
   - Shows combined filters
   - Verifies expected pass/fail behavior

### Running Tests

```bash
# Unit tests (original test suite)
pytest tests/test_filters.py -v  # All 94 tests pass

# New feature tests
python test_per_read_filters.py  # All per-read tests pass

# CLI demo
python demo_per_read_simple.py  # Shows working examples
```

## Use Cases

### When to Use Per-Read Thresholds

1. **Asymmetric paired-end reads**
   - R1 is longer/higher quality than R2
   - Different QC standards for each read

2. **Index reads**
   - Separate thresholds for data reads vs. barcode/index reads
   - Typically shorter with different quality requirements

3. **Targeted sequencing**
   - Forward and reverse reads cover different genomic regions
   - May have different quality or length expectations

4. **Legacy compatibility**
   - Matching QC from other tools that apply per-read thresholds
   - Recreating specific filtering pipelines

### Example: Real-World Scenario

```bash
# Illumina NovaSeq paired-end 2×150 with dual index
# R1/R2: 150bp target, Q30 minimum
# I1/I2: 8bp barcodes, Q20 minimum

fastq-filter \\
  -l 140,140,6,6 \\        # Min lengths: R1≥140, R2≥140, I1≥6, I2≥6
  -Q 30,30,20,20 \\         # Median Q: R1≥30, R2≥30, I1≥20, I2≥20
  -o r1.fq -o r2.fq -o i1.fq -o i2.fq \\
  r1_in.fq r2_in.fq i1_in.fq i2_in.fq
```

## Performance

- Per-read filtering adds minimal overhead:
  - Array length check (one comparison)
  - Individual metric calculation per read (was already looping)
  - Early termination on first failure (faster than aggregate)
- Single threshold behavior unchanged (same performance)

## Migration Guide

### Existing Code
No changes needed - single thresholds work exactly as before:

```python
# This continues to work
filt = MinimumLengthFilter(50)
filt((r1, r2))
```

### Adopting Per-Read Thresholds

```python
# Before: Single threshold applied to all reads
old_filter = MinimumLengthFilter(50)
# Passed if R1≥50 OR R2≥50 (either one)

# After: Per-read thresholds
new_filter = MinimumLengthFilter((80, 50))
# Passes only if R1≥80 AND R2≥50 (both must pass)
```

**Note**: Behavior changes when switching from single to per-read:
- Length filters change from "any" to "all" logic
- Quality filters change from aggregate to per-read evaluation
- May filter more reads - review thresholds carefully

## Summary

✓ **Backward compatible** - all existing tests pass  
✓ **API flexible** - accepts single value or tuple  
✓ **CLI intuitive** - comma-separated values  
✓ **Well-tested** - comprehensive test suite  
✓ **Documented** - examples and use cases provided  

The per-read filtering feature enables more precise QC for multi-read sequencing data while maintaining full compatibility with existing workflows.
