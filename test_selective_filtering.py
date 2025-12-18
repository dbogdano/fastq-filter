#!/usr/bin/env python3
"""Test selective read filtering with filter_indices."""

import tempfile
import subprocess
from pathlib import Path
import sys

def create_test_files(tmpdir):
    """Create paired-end FASTQ files where R1 has high quality and R2 has low quality."""
    r1 = Path(tmpdir) / "R1.fastq"
    r2 = Path(tmpdir) / "R2.fastq"
    
    # R1: All high quality Q40
    r1.write_text(
        "@read1\n" + "A" * 100 + "\n+\n" + "I" * 100 + "\n" +
        "@read2\n" + "T" * 100 + "\n+\n" + "I" * 100 + "\n" +
        "@read3\n" + "G" * 100 + "\n+\n" + "I" * 100 + "\n"
    )
    
    # R2: Mixed quality - read1 and read2 low (Q20), read3 high (Q40)
    r2.write_text(
        "@read1\n" + "C" * 80 + "\n+\n" + "5" * 80 + "\n" +  # Low quality Q20
        "@read2\n" + "G" * 80 + "\n+\n" + "5" * 80 + "\n" +  # Low quality Q20
        "@read3\n" + "A" * 80 + "\n+\n" + "I" * 80 + "\n"    # High quality Q40
    )
    
    return r1, r2

def run_filter(python_exe, inputs, outputs, args):
    """Run fastq_filter."""
    cmd = [python_exe, "-m", "fastq_filter", *args]
    for out in outputs:
        cmd.extend(["-o", str(out)])
    cmd.extend([str(f) for f in inputs])
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result

def count_reads(filepath):
    """Count reads in FASTQ file."""
    if not filepath.exists():
        return 0
    return filepath.read_text().count("@read")

def main():
    tmpdir = tempfile.mkdtemp(prefix="selective_filter_test_")
    print(f"Testing selective read filtering in: {tmpdir}\n")
    
    r1_in, r2_in = create_test_files(tmpdir)
    r1_out = Path(tmpdir) / "R1_out.fastq"
    r2_out = Path(tmpdir) / "R2_out.fastq"
    
    python_exe = sys.executable
    
    print("=" * 70)
    print("INPUT DATA:")
    print("=" * 70)
    print("R1: 100bp, All reads Q40 (high quality)")
    print("R2:  80bp, read1 Q20, read2 Q20, read3 Q40")
    
    # Test 1: Filter both reads (default behavior with per-read thresholds)
    print("\n" + "=" * 70)
    print("TEST 1: Filter both R1 and R2 (default)")
    print("=" * 70)
    print("Command: -Q 30,30 (both must be ≥Q30)")
    print("Expected: Only read3 passes (R1 Q40✓, R2 Q40✓)")
    
    result = run_filter(
        python_exe,
        [r1_in, r2_in],
        [r1_out, r2_out],
        ["-Q", "30,30"]
    )
    
    if result.returncode == 0:
        count = count_reads(r1_out)
        print(f"✓ Result: {count}/3 pairs passed")
        assert count == 1, f"Expected 1, got {count}"
    else:
        print(f"✗ Failed: {result.stderr}")
        return
    
    # Test 2: Filter only R1
    print("\n" + "=" * 70)
    print("TEST 2: Filter only R1, keep pairs synchronized")
    print("=" * 70)
    print("Command: -Q 30 --apply-filters-to 1")
    print("Expected: All 3 pairs pass (only checking R1 Q≥30, all pass)")
    
    result = run_filter(
        python_exe,
        [r1_in, r2_in],
        [r1_out, r2_out],
        ["-Q", "30", "--apply-filters-to", "1"]
    )
    
    if result.returncode == 0:
        count = count_reads(r1_out)
        print(f"✓ Result: {count}/3 pairs passed")
        assert count == 3, f"Expected 3, got {count}"
        # Verify R2 also has 3 reads (pairs kept in sync)
        count_r2 = count_reads(r2_out)
        print(f"  R2 also has {count_r2}/3 reads (pairs synchronized)")
        assert count_r2 == 3, f"R2 should also have 3 reads"
    else:
        print(f"✗ Failed: {result.stderr}")
        return
    
    # Test 3: Filter only R2
    print("\n" + "=" * 70)
    print("TEST 3: Filter only R2, discard pairs with low R2 quality")
    print("=" * 70)
    print("Command: -Q 30 --apply-filters-to 2")
    print("Expected: Only read3 passes (R2 must be Q≥30)")
    
    result = run_filter(
        python_exe,
        [r1_in, r2_in],
        [r1_out, r2_out],
        ["-Q", "30", "--apply-filters-to", "2"]
    )
    
    if result.returncode == 0:
        count = count_reads(r1_out)
        print(f"✓ Result: {count}/3 pairs passed")
        assert count == 1, f"Expected 1, got {count}"
        # Show which read passed
        if r1_out.exists():
            content = r1_out.read_text()
            passed_reads = [line[1:] for line in content.split('\n') if line.startswith('@')]
            print(f"  Passed: {', '.join(passed_reads)}")
        # Verify R1 also only has 1 read (pair discarded)
        count_r1 = count_reads(r1_out)
        print(f"  R1 also has {count_r1}/3 reads (pair discarded when R2 fails)")
        assert count_r1 == 1, f"R1 should also have 1 read"
    else:
        print(f"✗ Failed: {result.stderr}")
        return
    
    # Test 4: Per-read thresholds with selective filtering
    print("\n" + "=" * 70)
    print("TEST 4: Per-read thresholds + selective filtering")
    print("=" * 70)
    print("Command: -Q 35,25 --apply-filters-to 2")
    print("Expected: All 3 pairs pass (only checking R2, all R2 ≥Q20)")
    
    result = run_filter(
        python_exe,
        [r1_in, r2_in],
        [r1_out, r2_out],
        ["-Q", "35,25", "--apply-filters-to", "2"]
    )
    
    if result.returncode == 0:
        count = count_reads(r1_out)
        print(f"✓ Result: {count}/3 pairs passed")
        print(f"  (Threshold for R2 is 25, all R2 reads are ≥Q20, so... wait)")
        print(f"  Actually read1 and read2 have R2 Q20 < 25, so expect 1 pass")
        # Correction: Q20 < 25, so only read3 should pass
        assert count == 1, f"Expected 1, got {count}"
    else:
        print(f"✗ Failed: {result.stderr}")
        return
    
    print("\n" + "=" * 70)
    print("✓ ALL TESTS PASSED")
    print("=" * 70)
    print("\nKEY FEATURE: --apply-filters-to lets you:")
    print("  • Filter based on specific reads (e.g., only R1)")
    print("  • Discard entire read tuples when that read fails")
    print("  • Keep pairs synchronized across all files")
    print(f"\nCleanup: rm -rf {tmpdir}")

if __name__ == "__main__":
    main()
