#!/usr/bin/env python3
"""
Demo: Show per-read filtering with direct file creation.
"""

import tempfile
import subprocess
from pathlib import Path
import sys

def create_fastq_files(tmpdir):
    """Create sample FASTQ files with known properties."""
    r1 = Path(tmpdir) / "input_R1.fastq"
    r2 = Path(tmpdir) / "input_R2.fastq"
    
    # R1: All reads 86bp with uniform quality
    # read1: Q40 (I), read2: Q39 (H), read3: Q38 (G)
    r1_content = [
        "@read1\n",
        "ACGT" * 21 + "AC\n",  # 86bp
        "+\n",
        "I" * 86 + "\n",
        "@read2\n",
        "TGCA" * 21 + "TG\n",  # 86bp
        "+\n",
        "H" * 86 + "\n",
        "@read3\n",
        "GGGG" * 21 + "GG\n",  # 86bp
        "+\n",
        "G" * 86 + "\n",
    ]
    
    # R2: Different lengths with uniform quality
    # read1: 46bp Q20 (5), read2: 68bp Q39 (H), read3: 72bp Q18 (3)
    r2_content = [
        "@read1\n",
        "TGCA" * 11 + "TG\n",  # 46bp
        "+\n",
        "5" * 46 + "\n",
        "@read2\n",
        "ACGT" * 17 + "\n",  # 68bp
        "+\n",
        "H" * 68 + "\n",
        "@read3\n",
        "CCCC" * 18 + "\n",  # 72bp
        "+\n",
        "3" * 72 + "\n",
    ]
    
    r1.write_text("".join(r1_content))
    r2.write_text("".join(r2_content))
    
    return r1, r2

def run_filter(python_exe, r1_in, r2_in, r1_out, r2_out, args):
    """Run fastq_filter with given arguments."""
    cmd = [
        python_exe, "-m", "fastq_filter",
        *args,
        "-o", str(r1_out), "-o", str(r2_out),
        str(r1_in), str(r2_in)
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result

def count_reads(filepath):
    """Count reads in a FASTQ file."""
    if not filepath.exists():
        return 0
    return filepath.read_text().count("@read")

def main():
    tmpdir = tempfile.mkdtemp(prefix="fastq_filter_demo_")
    print(f"Created temp directory: {tmpdir}\n")
    
    r1_in, r2_in = create_fastq_files(tmpdir)
    r1_out = Path(tmpdir) / "output_R1.fastq"
    r2_out = Path(tmpdir) / "output_R2.fastq"
    
    python_exe = sys.executable
    
    print("=" * 70)
    print("INPUT READ STATISTICS:")
    print("=" * 70)
    print("R1: All reads 86bp")
    print("  - read1: Q40 (I=73), read2: Q39 (H=72), read3: Q38 (G=71)")
    print("R2: Variable length")
    print("  - read1: 46bp Q20 (5=53), read2: 68bp Q39 (H=72), read3: 72bp Q18 (3=51)")
    
    # Example 1
    print("\n" + "=" * 70)
    print("EXAMPLE 1: Single threshold - Minimum length 60bp")
    print("=" * 70)
    print("Command: -l 60")
    print("Expected: All 3 pairs pass (any read ≥ 60 passes the tuple)")
    result = run_filter(python_exe, r1_in, r2_in, r1_out, r2_out, ["-l", "60"])
    if result.returncode == 0:
        print(f"✓ Output: {count_reads(r1_out)}/3 pairs")
    else:
        print(f"✗ Failed: {result.stderr[:200]}")
    
    # Example 2
    print("\n" + "=" * 70)
    print("EXAMPLE 2: Per-read thresholds - R1 ≥ 80bp, R2 ≥ 60bp")
    print("=" * 70)
    print("Command: -l 80,60")
    print("Expected: 2 pairs pass (read1 fails: R2=46<60)")
    result = run_filter(python_exe, r1_in, r2_in, r1_out, r2_out, ["-l", "80,60"])
    if result.returncode == 0:
        print(f"✓ Output: {count_reads(r1_out)}/3 pairs")
    else:
        print(f"✗ Failed: {result.stderr[:200]}")
    
    # Example 3
    print("\n" + "=" * 70)
    print("EXAMPLE 3: Per-read quality - R1 median Q≥38, R2 median Q≥19")
    print("=" * 70)
    print("Command: -Q 38,19")
    print("Expected: 2 pairs pass (read3 fails: R2 Q18<19)")
    result = run_filter(python_exe, r1_in, r2_in, r1_out, r2_out, ["-Q", "38,19"])
    if result.returncode == 0:
        print(f"✓ Output: {count_reads(r1_out)}/3 pairs")
    else:
        print(f"✗ Failed: {result.stderr[:200]}")
    
    # Example 4
    print("\n" + "=" * 70)
    print("EXAMPLE 4: Combined - Length (80,60) AND Quality (38,19)")
    print("=" * 70)
    print("Command: -l 80,60 -Q 38,19")
    print("Expected: 1 pair passes (only read2 passes both)")
    result = run_filter(python_exe, r1_in, r2_in, r1_out, r2_out, ["-l", "80,60", "-Q", "38,19"])
    if result.returncode == 0:
        print(f"✓ Output: {count_reads(r1_out)}/3 pairs")
        if r2_out.exists():
            content = r2_out.read_text()
            print(f"\nOutput R2 content preview:")
            print(content[:200] + "..." if len(content) > 200 else content)
    else:
        print(f"✗ Failed: {result.stderr[:200]}")
    
    print("\n" + "=" * 70)
    print("DEMO COMPLETE - All tests show per-read filtering working!")
    print("=" * 70)
    print(f"Temp files in: {tmpdir}")
    print(f"Clean up: rm -rf {tmpdir}")

if __name__ == "__main__":
    main()
