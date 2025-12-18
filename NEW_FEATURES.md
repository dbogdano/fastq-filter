# New Features in fastq-filter

## Per-Read Thresholds

You can now specify different quality/length thresholds for each read in paired-end or multi-read data by providing comma-separated values.

### Examples:

**Quality filtering with different thresholds:**
```bash
# R1 must have Q≥30, R2 must have Q≥25
fastq-filter -Q 30,25 -i R1.fastq R2.fastq -o filtered_R1.fastq filtered_R2.fastq
```

**Length filtering with different thresholds:**
```bash
# R1 must be ≥95bp, R2 must be ≥75bp, R3 must be ≥8bp
fastq-filter -l 95,75,8 -i R1.fastq R2.fastq R3.fastq -o filtered_R1.fastq filtered_R2.fastq filtered_R3.fastq
```

## Selective Read Filtering

You can now filter based on specific reads while keeping all paired reads synchronized using the `--apply-filters-to` option.

This is useful when you only care about the quality of one read (e.g., R1) but want to discard both reads in a pair when that read fails quality checks.

### Examples:

**Filter only R2, discard pairs when R2 fails:**
```bash
# Only check R2 quality (must be ≥Q30)
# If R2 fails, discard both R1 and R2
fastq-filter -Q 30 --apply-filters-to 2 -i R1.fastq R2.fastq -o filtered_R1.fastq filtered_R2.fastq
```

**Filter only R1, keep pairs synchronized:**
```bash
# Only check R1 quality (must be ≥Q30)
# If R1 fails, discard both R1 and R2
fastq-filter -Q 30 --apply-filters-to 1 -i R1.fastq R2.fastq -o filtered_R1.fastq filtered_R2.fastq
```

**Filter multiple specific reads:**
```bash
# Only check R1 and R3 (skip R2 quality check)
fastq-filter -Q 30 --apply-filters-to 1,3 -i R1.fastq R2.fastq R3.fastq -o filtered_R1.fastq filtered_R2.fastq filtered_R3.fastq
```

### Combining Features:

You can combine per-read thresholds with selective filtering:

```bash
# Use different thresholds (R1=35, R2=25)
# But only check R2 (the second read)
# Result: only R2 needs to be ≥Q25, R1 is kept regardless of its quality
fastq-filter -Q 35,25 --apply-filters-to 2 -i R1.fastq R2.fastq -o filtered_R1.fastq filtered_R2.fastq
```

## Important Notes:

1. **Read indices are 1-based** in the CLI (R1=1, R2=2, etc.)
2. **Pairs stay synchronized**: If any checked read fails, the entire read tuple is discarded
3. **Backward compatible**: Existing single-threshold behavior is unchanged
4. **Works with all filters**: `-Q`, `-q`, `-l`, `-L` all support these features
5. **Multi-file support**: Works with 2, 3, or more input files

## Use Cases:

### Use Case 1: Index Read Quality Control
Filter based on index read (R3) quality while keeping R1 and R2:
```bash
fastq-filter -Q 30 --apply-filters-to 3 -i R1.fastq R2.fastq I1.fastq -o filtered_R1.fastq filtered_R2.fastq filtered_I1.fastq
```

### Use Case 2: Asymmetric Paired-End Quality
R1 typically has better quality than R2. Require high quality for R1, lower for R2:
```bash
fastq-filter -Q 35,25 -i R1.fastq R2.fastq -o filtered_R1.fastq filtered_R2.fastq
```

### Use Case 3: Variable Length Trimming
After adapter trimming, reads may have different lengths:
```bash
# R1 minimum 95bp, R2 minimum 75bp (R2 is shorter after trimming)
fastq-filter -l 95,75 -i R1_trimmed.fastq R2_trimmed.fastq -o filtered_R1.fastq filtered_R2.fastq
```
