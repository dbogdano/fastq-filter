#!/usr/bin/env python3
"""
Demo: Create sample FASTQ files and test per-read filtering from the command line.
"""

import tempfile
import os
from pathlib import Path

# Sample paired-end reads (sequence and quality must be same length)
# R1: All reads 86bp long with different qualities
FASTQ_R1 = """@read1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
@read3
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
+
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
"""

# R2: Reads of different lengths with different qualities
FASTQ_R2 = """@read1
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
+
5555555555555555555555555555555555555555555555
@read2
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
@read3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
+
333333333333333333333333333333333333333333333333333333333333333333333333
"""

def main():
    # Create temp directory
    tmpdir = tempfile.mkdtemp(prefix="fastq_filter_demo_")
    print(f"Created temp directory: {tmpdir}")
    
    # Write input files
    r1_in = Path(tmpdir) / "input_R1.fastq"
    r2_in = Path(tmpdir) / "input_R2.fastq"
    r1_out = Path(tmpdir) / "output_R1.fastq"
    r2_out = Path(tmpdir) / "output_R2.fastq"
    
    r1_in.write_text(FASTQ_R1)
    r2_in.write_text(FASTQ_R2)
    print(f"Created input files: {r1_in} and {r2_in}")
    
    # Show read statistics
    print("\n" + "="*70)
    print("INPUT READ STATISTICS:")
    print("="*70)
    print("Read 1 (R1):")
    print("  - read1: 86bp, Quality Q40 (I = ASCII 73, 73-33=40)")
    print("  - read2: 86bp, Quality Q39 (H = ASCII 72, 72-33=39)")
    print("  - read3: 86bp, Quality Q38 (G = ASCII 71, 71-33=38)")
    print("\nRead 2 (R2):")
    print("  - read1: 46bp, Quality Q20 (5 = ASCII 53, 53-33=20)")
    print("  - read2: 68bp, Quality Q39 (H = ASCII 72, 72-33=39)")
    print("  - read3: 72bp, Quality Q18 (3 = ASCII 51, 51-33=18)")
    
    # Get python executable
    import sys
    python_exe = sys.executable
    
    # Example 1: Single threshold (all reads must meet same criteria)
    print("\n" + "="*70)
    print("EXAMPLE 1: Single threshold - Minimum length 60bp")
    print("="*70)
    cmd = (
        f"{python_exe} -m fastq_filter "
        f"-l 60 "
        f"-o {r1_out} -o {r2_out} "
        f"{r1_in} {r2_in}"
    )
    print(f"Command: {cmd}")
    print("\nExpected: All 3 pairs pass (R1=86bp, R2=46-72bp)")
    print("          With single threshold, passes if ANY read >= 60")
    os.system(cmd)
    
    # Count output reads
    output_count = r1_out.read_text().count("@read")
    print(f"\n✓ Output reads: {output_count}/3")
    
    # Example 2: Per-read thresholds
    print("\n" + "="*70)
    print("EXAMPLE 2: Per-read thresholds - R1 ≥ 80bp, R2 ≥ 60bp")
    print("="*70)
    cmd = (
        f"{python_exe} -m fastq_filter "
        f"-l 80,60 "
        f"-o {r1_out} -o {r2_out} "
        f"{r1_in} {r2_in}"
    )
    print(f"Command: {cmd}")
    print("\nExpected: 2 pairs pass:")
    print("  - read1: ❌ R1=86✓ (≥80) but R2=46✗ (< 60)")
    print("  - read2: ✓ R1=86✓ (≥80) and R2=68✓ (≥60)")
    print("  - read3: ✓ R1=86✓ (≥80) and R2=72✓ (≥60)")
    os.system(cmd)
    
    output_count = r1_out.read_text().count("@read")
    print(f"\n✓ Output reads: {output_count}/3")
    
    # Example 3: Per-read quality thresholds
    print("\n" + "="*70)
    print("EXAMPLE 3: Per-read quality - R1 median Q≥38, R2 median Q≥19")
    print("="*70)
    cmd = (
        f"{python_exe} -m fastq_filter "
        f"-Q 38,19 "
        f"-o {r1_out} -o {r2_out} "
        f"{r1_in} {r2_in}"
    )
    print(f"Command: {cmd}")
    print("\nExpected: 2 pairs pass:")
    print("  - read1: ✓ R1 Q40≥38 ✓ and R2 Q20≥19 ✓")
    print("  - read2: ✓ R1 Q39≥38 ✓ and R2 Q39≥19 ✓")
    print("  - read3: ❌ R1 Q38≥38 ✓ but R2 Q18<19 ✗")
    os.system(cmd)
    
    output_count = r1_out.read_text().count("@read")
    print(f"\n✓ Output reads: {output_count}/3")
    
    # Example 4: Combined filters with per-read thresholds
    print("\n" + "="*70)
    print("EXAMPLE 4: Combined - Length (80,60) AND Quality (38,19)")
    print("="*70)
    cmd = (
        f"{python_exe} -m fastq_filter "
        f"-l 80,60 -Q 38,19 "
        f"-o {r1_out} -o {r2_out} "
        f"{r1_in} {r2_in}"
    )
    print(f"Command: {cmd}")
    print("\nExpected: Only read2 passes both filters")
    os.system(cmd)
    
    output_count = r1_out.read_text().count("@read")
    print(f"\n✓ Output reads: {output_count}/3")
    
    print("\n" + "="*70)
    print("DEMO COMPLETE")
    print("="*70)
    print(f"\nTemp files are in: {tmpdir}")
    print("You can inspect the output files manually if desired.")
    print(f"To clean up: rm -rf {tmpdir}")

if __name__ == "__main__":
    main()
